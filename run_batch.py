#!/usr/bin/env python3
"""
Batch SAXS profile computation using FoXS.

Computes theoretical SAXS profiles for all PDB/CIF files in a directory.

Usage:
    python run_batch.py --pdb-dir /path/to/pdbs
    python run_batch.py --pdb-dir /path/to/pdbs --plot
    python run_batch.py --pdb-dir /path/to/pdbs --plot --max-traces 10
    python run_batch.py --pdb-dir /path/to/pdbs --max-q 0.5 --num-points 500
    python run_batch.py --pdb-dir /path/to/pdbs --workers 8
"""

import argparse
import csv
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

# Path to FoXS binary in the conda environment (mutable config)
_CONFIG = {
    "foxs_bin": Path.home() / "miniconda3" / "envs" / "foxs_env" / "bin" / "foxs",
}

# Default output directory
DEFAULT_OUTPUT_DIR = Path(__file__).resolve().parent / "output"


def _set_foxs_bin(path: Path):
    """Update the FoXS binary path."""
    _CONFIG["foxs_bin"] = path


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
# Single job
# ---------------------------------------------------------------------------

def process_pdb_file(
    pdb_path: Path,
    saxs_dir: Path,
    max_q: float,
    num_points: int,
    hydrogens: bool,
) -> dict:
    """Process a single PDB file through FoXS."""
    result = {
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
    fieldnames = ["pdb_file", "saxs_file", "status", "message"]
    with open(summary_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(results)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Batch SAXS profile computation using FoXS.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  python run_batch.py --pdb-dir /path/to/pdbs\n"
            "  python run_batch.py --pdb-dir /path/to/pdbs --plot\n"
            "  python run_batch.py --pdb-dir /path/to/pdbs --plot --max-traces 10\n"
            "  python run_batch.py --pdb-dir /path/to/pdbs --max-q 0.5 --num-points 500\n"
            "  python run_batch.py --pdb-dir /path/to/pdbs --workers 8\n"
        ),
    )

    parser.add_argument(
        "--pdb-dir",
        type=str,
        required=True,
        help="Directory containing PDB files to process.",
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

    # Plotting
    parser.add_argument(
        "--plot",
        action="store_true",
        help="Generate an overlay plot of all SAXS profiles after the run.",
    )
    parser.add_argument(
        "--max-traces",
        type=int,
        default=50,
        help="Maximum number of profiles to include in the plot (default: 50).",
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

    pdb_dir = Path(args.pdb_dir)
    if not pdb_dir.is_dir():
        print(f"ERROR: {pdb_dir} is not a directory")
        sys.exit(1)

    pdb_files = sorted([*pdb_dir.glob("*.pdb"), *pdb_dir.glob("*.cif")])
    if not pdb_files:
        print(f"No PDB or CIF files found in {pdb_dir}")
        sys.exit(1)

    saxs_dir = output_dir / "saxs"
    saxs_dir.mkdir(parents=True, exist_ok=True)

    print(f"Processing {len(pdb_files)} structures from {pdb_dir}")
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
                    "pdb_file": pdb_path.name, "saxs_file": "",
                    "status": "error", "message": str(e),
                })

    elapsed = time.time() - t0
    summary_path = output_dir / "summary.csv"
    write_summary(results, summary_path)
    _print_summary(results, elapsed, summary_path)

    if args.plot:
        from plot_profiles import plot_overlay, collect_dat_files
        dat_files = collect_dat_files([str(saxs_dir)])
        if dat_files:
            plot_path = output_dir / "saxs_overlay.png"
            title = Path(args.pdb_dir).name + " — SAXS profiles"
            print()
            plot_overlay(dat_files, plot_path, title=title, max_traces=args.max_traces)
        else:
            print("No .dat files to plot.")


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
                print(f"  {r['pdb_file']}: {r['message']}")


if __name__ == "__main__":
    main()
