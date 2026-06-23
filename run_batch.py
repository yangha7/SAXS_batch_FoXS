#!/usr/bin/env python3
"""
Batch SAXS profile computation for DNA oligomers using FoXS.

Generates all possible DNA oligomer sequences of a given length,
builds 3D structures via DNA_Builder, and computes theoretical
SAXS profiles using FoXS (IMP).

Usage:
    python run_batch.py --length 4 --form B
    python run_batch.py --length 6 --form A B Z
    python run_batch.py --sequences ATCGATCG GCGCGCGC --form B
    python run_batch.py --pdb-dir /path/to/pdbs
"""

import argparse
import csv
import itertools
import os
import shutil
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

# Path to DNA_Builder (sibling directory)
DNA_BUILDER_DIR = Path(__file__).resolve().parent.parent / "DNA_Builder"

# Path to FoXS binary in the conda environment (mutable config)
_CONFIG = {
    "foxs_bin": Path.home() / "miniconda3" / "envs" / "foxs_env" / "bin" / "foxs",
}

# Default output directory
DEFAULT_OUTPUT_DIR = Path(__file__).resolve().parent / "output"

BASES = ["A", "T", "G", "C"]


def _set_foxs_bin(path: Path):
    """Update the FoXS binary path."""
    _CONFIG["foxs_bin"] = path


# ---------------------------------------------------------------------------
# Sequence generation
# ---------------------------------------------------------------------------

def generate_all_sequences(length: int) -> list[str]:
    """Generate all possible DNA sequences of the given length."""
    return ["".join(combo) for combo in itertools.product(BASES, repeat=length)]


def count_sequences(length: int) -> int:
    """Return 4^length."""
    return 4 ** length


# ---------------------------------------------------------------------------
# PDB generation via DNA_Builder
# ---------------------------------------------------------------------------

def build_pdb(sequence: str, form: str, output_path: Path) -> tuple[bool, str]:
    """
    Build a PDB file for the given sequence and form using DNA_Builder.

    Returns (success: bool, message: str).
    """
    try:
        result = subprocess.run(
            [
                sys.executable, "-m", "dna_builder",
                sequence,
                "--form", form,
                "--output", str(output_path),
            ],
            capture_output=True,
            text=True,
            timeout=60,
            cwd=str(DNA_BUILDER_DIR),
        )
        if result.returncode != 0:
            return False, f"DNA_Builder error: {result.stderr.strip()}"
        return True, result.stdout.strip()
    except subprocess.TimeoutExpired:
        return False, "DNA_Builder timed out"
    except Exception as e:
        return False, f"DNA_Builder exception: {e}"


# ---------------------------------------------------------------------------
# FoXS SAXS computation
# ---------------------------------------------------------------------------

def run_foxs(
    pdb_path: Path,
    output_dir: Path,
    max_q: float = 0.5,
    num_points: int = 500,
    hydrogens: bool = False,
) -> tuple[bool, str, Path | None]:
    """
    Run FoXS on a PDB file to compute the SAXS profile.

    FoXS writes the output .dat file in the current working directory
    with the name <pdb_filename>.dat.

    Returns (success, message, dat_path).
    """
    foxs_bin = _CONFIG["foxs_bin"]
    if not foxs_bin.exists():
        return False, f"FoXS not found at {foxs_bin}", None

    cmd = [
        str(foxs_bin),
        "-s", str(num_points),
        "-q", str(max_q),
        str(pdb_path),
    ]
    if hydrogens:
        cmd.insert(-1, "-h")

    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=120,
            cwd=str(output_dir),
        )

        # FoXS creates <filename>.dat in the cwd
        dat_name = pdb_path.name + ".dat"
        dat_path = output_dir / dat_name

        if dat_path.exists():
            return True, result.stdout.strip() + result.stderr.strip(), dat_path
        else:
            # Sometimes FoXS puts it next to the PDB
            alt_dat = pdb_path.parent / dat_name
            if alt_dat.exists():
                final = output_dir / dat_name
                shutil.move(str(alt_dat), str(final))
                return True, result.stdout.strip() + result.stderr.strip(), final
            return False, f"FoXS ran but no .dat file found. stderr: {result.stderr.strip()}", None

    except subprocess.TimeoutExpired:
        return False, "FoXS timed out", None
    except Exception as e:
        return False, f"FoXS exception: {e}", None


# ---------------------------------------------------------------------------
# Single job: build PDB + run FoXS
# ---------------------------------------------------------------------------

def process_sequence(
    sequence: str,
    form: str,
    pdb_dir: Path,
    saxs_dir: Path,
    max_q: float,
    num_points: int,
    hydrogens: bool,
) -> dict:
    """
    Process a single sequence: build PDB, run FoXS, return result dict.
    """
    form_upper = form.upper()
    pdb_name = f"{form_upper}_{sequence}.pdb"
    pdb_path = pdb_dir / pdb_name

    result = {
        "sequence": sequence,
        "form": form_upper,
        "pdb_file": pdb_name,
        "saxs_file": "",
        "status": "pending",
        "message": "",
    }

    # Step 1: Build PDB
    if not pdb_path.exists():
        ok, msg = build_pdb(sequence, form_upper, pdb_path)
        if not ok:
            result["status"] = "pdb_error"
            result["message"] = msg
            return result
    else:
        msg = "PDB already exists, skipping build"

    # Step 2: Run FoXS
    dat_expected = saxs_dir / f"{pdb_name}.dat"
    if dat_expected.exists():
        result["saxs_file"] = dat_expected.name
        result["status"] = "skipped"
        result["message"] = "SAXS profile already exists"
        return result

    ok, msg, dat_path = run_foxs(pdb_path, saxs_dir, max_q, num_points, hydrogens)
    if ok and dat_path:
        result["saxs_file"] = dat_path.name
        result["status"] = "success"
        result["message"] = "OK"
    else:
        result["status"] = "foxs_error"
        result["message"] = msg

    return result


# ---------------------------------------------------------------------------
# Batch processing for existing PDB directory
# ---------------------------------------------------------------------------

def process_pdb_file(
    pdb_path: Path,
    saxs_dir: Path,
    max_q: float,
    num_points: int,
    hydrogens: bool,
) -> dict:
    """Process a single existing PDB file through FoXS."""
    result = {
        "sequence": pdb_path.stem,
        "form": "",
        "pdb_file": pdb_path.name,
        "saxs_file": "",
        "status": "pending",
        "message": "",
    }

    dat_expected = saxs_dir / f"{pdb_path.name}.dat"
    if dat_expected.exists():
        result["saxs_file"] = dat_expected.name
        result["status"] = "skipped"
        result["message"] = "SAXS profile already exists"
        return result

    ok, msg, dat_path = run_foxs(pdb_path, saxs_dir, max_q, num_points, hydrogens)
    if ok and dat_path:
        result["saxs_file"] = dat_path.name
        result["status"] = "success"
        result["message"] = "OK"
    else:
        result["status"] = "foxs_error"
        result["message"] = msg

    return result


# ---------------------------------------------------------------------------
# Summary CSV
# ---------------------------------------------------------------------------

def write_summary(results: list[dict], summary_path: Path):
    """Write a CSV summary of all results."""
    fieldnames = ["sequence", "form", "pdb_file", "saxs_file", "status", "message"]
    with open(summary_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(results)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Batch SAXS profile computation for DNA oligomers using FoXS.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  python run_batch.py --length 4 --form B\n"
            "  python run_batch.py --length 6 --form A B Z\n"
            "  python run_batch.py --sequences ATCGATCG GCGCGCGC --form B\n"
            "  python run_batch.py --pdb-dir /path/to/pdbs\n"
        ),
    )

    # Sequence source (mutually exclusive)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument(
        "--length", "-l",
        type=int,
        help="Generate all possible sequences of this length (4^N sequences).",
    )
    source.add_argument(
        "--sequences", "-s",
        nargs="+",
        type=str,
        help="Specific DNA sequences to process.",
    )
    source.add_argument(
        "--pdb-dir",
        type=str,
        help="Directory containing existing PDB files to process.",
    )

    # DNA form
    parser.add_argument(
        "--form", "-f",
        nargs="+",
        type=str,
        choices=["A", "B", "Z"],
        default=["B"],
        help="DNA form(s) to generate (default: B).",
    )

    # FoXS parameters
    parser.add_argument(
        "--max-q",
        type=float,
        default=0.5,
        help="Maximum q value in Å⁻¹ (default: 0.5).",
    )
    parser.add_argument(
        "--num-points",
        type=int,
        default=500,
        help="Number of points in the SAXS profile (default: 500).",
    )
    parser.add_argument(
        "--hydrogens",
        action="store_true",
        help="Explicitly consider hydrogens in PDB files.",
    )

    # Execution
    parser.add_argument(
        "--workers", "-w",
        type=int,
        default=None,
        help="Number of parallel workers (default: CPU count).",
    )
    parser.add_argument(
        "--output-dir", "-o",
        type=str,
        default=None,
        help=f"Output directory (default: {DEFAULT_OUTPUT_DIR}).",
    )

    # FoXS binary override
    parser.add_argument(
        "--foxs-bin",
        type=str,
        default=None,
        help=f"Path to FoXS binary (default: {_CONFIG['foxs_bin']}).",
    )

    args = parser.parse_args()

    # Override FoXS binary path if specified
    foxs_bin = Path(args.foxs_bin) if args.foxs_bin else _CONFIG["foxs_bin"]
    _set_foxs_bin(foxs_bin)

    # Validate FoXS exists
    if not foxs_bin.exists():
        print(f"ERROR: FoXS binary not found at {foxs_bin}")
        print("Install IMP via conda: conda install -c conda-forge imp")
        print("Or specify path with --foxs-bin")
        sys.exit(1)

    output_dir = Path(args.output_dir) if args.output_dir else DEFAULT_OUTPUT_DIR
    workers = args.workers or os.cpu_count() or 4

    # -----------------------------------------------------------------------
    # Mode 1: Process existing PDB directory
    # -----------------------------------------------------------------------
    if args.pdb_dir:
        pdb_dir = Path(args.pdb_dir)
        if not pdb_dir.is_dir():
            print(f"ERROR: {pdb_dir} is not a directory")
            sys.exit(1)

        pdb_files = sorted(pdb_dir.glob("*.pdb"))
        if not pdb_files:
            print(f"No PDB files found in {pdb_dir}")
            sys.exit(1)

        saxs_dir = output_dir / "saxs"
        saxs_dir.mkdir(parents=True, exist_ok=True)

        print(f"Processing {len(pdb_files)} PDB files from {pdb_dir}")
        print(f"SAXS output: {saxs_dir}")
        print(f"FoXS params: max_q={args.max_q}, points={args.num_points}")
        print(f"Workers: {workers}")
        print()

        results = []
        t0 = time.time()

        with ProcessPoolExecutor(max_workers=workers) as executor:
            futures = {
                executor.submit(
                    process_pdb_file,
                    pdb_path, saxs_dir,
                    args.max_q, args.num_points, args.hydrogens,
                ): pdb_path
                for pdb_path in pdb_files
            }

            for i, future in enumerate(as_completed(futures), 1):
                pdb_path = futures[future]
                try:
                    result = future.result()
                    results.append(result)
                    status_icon = "✓" if result["status"] in ("success", "skipped") else "✗"
                    print(f"  [{i}/{len(pdb_files)}] {status_icon} {pdb_path.name}: {result['status']}")
                except Exception as e:
                    print(f"  [{i}/{len(pdb_files)}] ✗ {pdb_path.name}: {e}")
                    results.append({
                        "sequence": pdb_path.stem, "form": "",
                        "pdb_file": pdb_path.name, "saxs_file": "",
                        "status": "error", "message": str(e),
                    })

        elapsed = time.time() - t0
        summary_path = output_dir / "summary.csv"
        write_summary(results, summary_path)
        _print_summary(results, elapsed, summary_path)
        return

    # -----------------------------------------------------------------------
    # Mode 2: Generate sequences and process
    # -----------------------------------------------------------------------
    if args.length:
        sequences = generate_all_sequences(args.length)
        print(f"Generated {len(sequences)} sequences of length {args.length}")
    else:
        sequences = [s.upper() for s in args.sequences]
        # Validate sequences
        for seq in sequences:
            if not all(c in "ATGC" for c in seq):
                print(f"ERROR: Invalid sequence '{seq}'. Only A, T, G, C allowed.")
                sys.exit(1)
        print(f"Processing {len(sequences)} specified sequences")

    forms = [f.upper() for f in args.form]
    total_jobs = len(sequences) * len(forms)

    print(f"DNA forms: {', '.join(forms)}")
    print(f"Total jobs: {total_jobs}")
    print(f"FoXS params: max_q={args.max_q}, points={args.num_points}")
    print(f"Workers: {workers}")
    print()

    all_results = []
    t0 = time.time()

    for form in forms:
        form_dir = output_dir / f"{form}_form"
        pdb_dir = form_dir / "pdb"
        saxs_dir = form_dir / "saxs"
        pdb_dir.mkdir(parents=True, exist_ok=True)
        saxs_dir.mkdir(parents=True, exist_ok=True)

        print(f"--- {form}-form DNA ---")
        results = []

        with ProcessPoolExecutor(max_workers=workers) as executor:
            futures = {
                executor.submit(
                    process_sequence,
                    seq, form, pdb_dir, saxs_dir,
                    args.max_q, args.num_points, args.hydrogens,
                ): seq
                for seq in sequences
            }

            for i, future in enumerate(as_completed(futures), 1):
                seq = futures[future]
                try:
                    result = future.result()
                    results.append(result)
                    status_icon = "✓" if result["status"] in ("success", "skipped") else "✗"
                    print(f"  [{i}/{len(sequences)}] {status_icon} {form}_{seq}: {result['status']}")
                except Exception as e:
                    print(f"  [{i}/{len(sequences)}] ✗ {form}_{seq}: {e}")
                    results.append({
                        "sequence": seq, "form": form,
                        "pdb_file": f"{form}_{seq}.pdb", "saxs_file": "",
                        "status": "error", "message": str(e),
                    })

        # Write per-form summary
        summary_path = form_dir / "summary.csv"
        write_summary(results, summary_path)
        all_results.extend(results)
        print()

    elapsed = time.time() - t0

    # Write global summary
    global_summary = output_dir / "summary_all.csv"
    write_summary(all_results, global_summary)
    _print_summary(all_results, elapsed, global_summary)


def _print_summary(results: list[dict], elapsed: float, summary_path: Path):
    """Print a summary of the batch run."""
    n_success = sum(1 for r in results if r["status"] == "success")
    n_skipped = sum(1 for r in results if r["status"] == "skipped")
    n_error = sum(1 for r in results if r["status"] not in ("success", "skipped"))

    print("=" * 60)
    print(f"Batch complete in {elapsed:.1f}s")
    print(f"  Success:  {n_success}")
    print(f"  Skipped:  {n_skipped}")
    print(f"  Errors:   {n_error}")
    print(f"  Summary:  {summary_path}")
    print("=" * 60)

    if n_error > 0:
        print("\nFailed jobs:")
        for r in results:
            if r["status"] not in ("success", "skipped"):
                print(f"  {r['sequence']} ({r['form']}): {r['message']}")


if __name__ == "__main__":
    main()
