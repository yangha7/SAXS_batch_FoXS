#!/usr/bin/env python3
"""
Plot SAXS profiles computed by run_batch.py.

Reads .dat files produced by FoXS and generates overlay plots,
Kratky plots, and Guinier plots for analysis.

Usage:
    python plot_profiles.py output/B_form/saxs/
    python plot_profiles.py output/B_form/saxs/ --kratky
    python plot_profiles.py output/B_form/saxs/ --guinier
    python plot_profiles.py output/B_form/saxs/ --select "ATCG*" --log-scale
    python plot_profiles.py file1.dat file2.dat --compare
"""

import argparse
import fnmatch
import sys
from pathlib import Path

import numpy as np

try:
    import matplotlib
    matplotlib.use("Agg")  # Non-interactive backend
    import matplotlib.pyplot as plt
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def load_foxs_dat(filepath: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Load a FoXS .dat file.

    Returns (q, intensity, error) arrays.
    """
    q_vals, i_vals, e_vals = [], [], []
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) >= 3:
                try:
                    q_vals.append(float(parts[0]))
                    i_vals.append(float(parts[1]))
                    e_vals.append(float(parts[2]))
                except ValueError:
                    continue
    return np.array(q_vals), np.array(i_vals), np.array(e_vals)


def collect_dat_files(paths: list[str], select: str | None = None) -> list[Path]:
    """
    Collect .dat files from the given paths (files or directories).
    Optionally filter by a glob pattern on the filename.
    """
    dat_files = []
    for p in paths:
        path = Path(p)
        if path.is_dir():
            dat_files.extend(sorted(path.glob("*.dat")))
        elif path.is_file() and path.suffix == ".dat":
            dat_files.append(path)
        else:
            print(f"Warning: skipping {p} (not a .dat file or directory)")

    if select:
        dat_files = [f for f in dat_files if fnmatch.fnmatch(f.stem, select)]

    return dat_files


# ---------------------------------------------------------------------------
# Plotting functions
# ---------------------------------------------------------------------------

def plot_overlay(
    dat_files: list[Path],
    output_path: Path,
    log_scale: bool = True,
    title: str = "SAXS Profiles",
    max_traces: int = 50,
):
    """Plot I(q) vs q for all profiles overlaid."""
    fig, ax = plt.subplots(figsize=(10, 7))

    if len(dat_files) > max_traces:
        print(f"  Plotting first {max_traces} of {len(dat_files)} profiles "
              f"(use --max-traces to change)")
        dat_files = dat_files[:max_traces]

    for dat_file in dat_files:
        q, intensity, _ = load_foxs_dat(dat_file)
        label = dat_file.stem.replace(".pdb", "")
        ax.plot(q, intensity, linewidth=0.8, alpha=0.7, label=label)

    ax.set_xlabel("q (Å⁻¹)", fontsize=12)
    ax.set_ylabel("I(q)", fontsize=12)
    ax.set_title(title, fontsize=14)

    if log_scale:
        ax.set_yscale("log")

    # Only show legend if few enough traces
    if len(dat_files) <= 20:
        ax.legend(fontsize=7, ncol=2, loc="upper right")

    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(output_path, dpi=150)
    plt.close(fig)
    print(f"  Saved: {output_path}")


def plot_kratky(
    dat_files: list[Path],
    output_path: Path,
    title: str = "Kratky Plot",
    max_traces: int = 50,
):
    """Plot q²·I(q) vs q (Kratky plot) for folding assessment."""
    fig, ax = plt.subplots(figsize=(10, 7))

    if len(dat_files) > max_traces:
        dat_files = dat_files[:max_traces]

    for dat_file in dat_files:
        q, intensity, _ = load_foxs_dat(dat_file)
        label = dat_file.stem.replace(".pdb", "")
        ax.plot(q, q**2 * intensity, linewidth=0.8, alpha=0.7, label=label)

    ax.set_xlabel("q (Å⁻¹)", fontsize=12)
    ax.set_ylabel("q² · I(q)", fontsize=12)
    ax.set_title(title, fontsize=14)

    if len(dat_files) <= 20:
        ax.legend(fontsize=7, ncol=2, loc="upper right")

    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(output_path, dpi=150)
    plt.close(fig)
    print(f"  Saved: {output_path}")


def plot_guinier(
    dat_files: list[Path],
    output_path: Path,
    q_max_guinier: float = 0.1,
    title: str = "Guinier Plot",
    max_traces: int = 50,
):
    """Plot ln(I(q)) vs q² (Guinier plot) for Rg estimation."""
    fig, ax = plt.subplots(figsize=(10, 7))

    if len(dat_files) > max_traces:
        dat_files = dat_files[:max_traces]

    for dat_file in dat_files:
        q, intensity, _ = load_foxs_dat(dat_file)
        mask = (q <= q_max_guinier) & (intensity > 0)
        q_g = q[mask]
        i_g = intensity[mask]
        label = dat_file.stem.replace(".pdb", "")
        ax.plot(q_g**2, np.log(i_g), linewidth=0.8, alpha=0.7, label=label)

    ax.set_xlabel("q² (Å⁻²)", fontsize=12)
    ax.set_ylabel("ln(I(q))", fontsize=12)
    ax.set_title(title, fontsize=14)

    if len(dat_files) <= 20:
        ax.legend(fontsize=7, ncol=2, loc="upper right")

    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(output_path, dpi=150)
    plt.close(fig)
    print(f"  Saved: {output_path}")


def plot_comparison(
    dat_files: list[Path],
    output_path: Path,
    title: str = "SAXS Profile Comparison",
):
    """
    Detailed comparison of 2-5 profiles with residuals.
    """
    if len(dat_files) < 2:
        print("Need at least 2 files for comparison")
        return
    if len(dat_files) > 5:
        print("Comparison limited to 5 profiles; using first 5")
        dat_files = dat_files[:5]

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 9), height_ratios=[3, 1],
                                    sharex=True)

    # Load reference (first file)
    q_ref, i_ref, _ = load_foxs_dat(dat_files[0])
    ref_label = dat_files[0].stem.replace(".pdb", "")
    ax1.plot(q_ref, i_ref, linewidth=1.5, label=f"{ref_label} (ref)")

    for dat_file in dat_files[1:]:
        q, intensity, _ = load_foxs_dat(dat_file)
        label = dat_file.stem.replace(".pdb", "")
        ax1.plot(q, intensity, linewidth=1.0, alpha=0.8, label=label)

        # Compute ratio to reference (interpolate if needed)
        if len(q) == len(q_ref) and np.allclose(q, q_ref):
            ratio = intensity / i_ref
        else:
            ratio = np.interp(q_ref, q, intensity) / i_ref
            q = q_ref
        ax2.plot(q, ratio, linewidth=1.0, alpha=0.8, label=label)

    ax1.set_ylabel("I(q)", fontsize=12)
    ax1.set_yscale("log")
    ax1.set_title(title, fontsize=14)
    ax1.legend(fontsize=9)
    ax1.grid(True, alpha=0.3)

    ax2.set_xlabel("q (Å⁻¹)", fontsize=12)
    ax2.set_ylabel("I(q) / I_ref(q)", fontsize=12)
    ax2.axhline(y=1.0, color="gray", linestyle="--", alpha=0.5)
    ax2.grid(True, alpha=0.3)

    fig.tight_layout()
    fig.savefig(output_path, dpi=150)
    plt.close(fig)
    print(f"  Saved: {output_path}")


# ---------------------------------------------------------------------------
# Rg extraction
# ---------------------------------------------------------------------------

def estimate_rg(dat_file: Path, q_max: float = 0.05) -> float | None:
    """
    Estimate radius of gyration (Rg) from Guinier analysis.

    Fits ln(I) = ln(I0) - (Rg²/3)·q² in the low-q region.
    """
    q, intensity, _ = load_foxs_dat(dat_file)
    mask = (q > 0) & (q <= q_max) & (intensity > 0)
    q_g = q[mask]
    i_g = intensity[mask]

    if len(q_g) < 5:
        return None

    # Linear fit: ln(I) = a + b·q²  where b = -Rg²/3
    coeffs = np.polyfit(q_g**2, np.log(i_g), 1)
    b = coeffs[0]

    if b >= 0:
        return None  # Non-physical

    rg = np.sqrt(-3.0 * b)
    return rg


def print_rg_table(dat_files: list[Path], q_max: float = 0.05):
    """Print a table of estimated Rg values."""
    print(f"\n{'Sequence':<30s} {'Rg (Å)':>10s}")
    print("-" * 42)
    for dat_file in dat_files:
        name = dat_file.stem.replace(".pdb", "")
        rg = estimate_rg(dat_file, q_max)
        if rg is not None:
            print(f"{name:<30s} {rg:>10.2f}")
        else:
            print(f"{name:<30s} {'N/A':>10s}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Plot and analyze SAXS profiles from FoXS .dat files.",
    )
    parser.add_argument(
        "paths",
        nargs="+",
        help="Paths to .dat files or directories containing .dat files.",
    )
    parser.add_argument(
        "--select",
        type=str,
        default=None,
        help="Glob pattern to filter filenames (e.g., 'B_ATCG*').",
    )
    parser.add_argument(
        "--output", "-o",
        type=str,
        default=None,
        help="Output image filename (default: auto-generated).",
    )
    parser.add_argument(
        "--log-scale",
        action="store_true",
        default=True,
        help="Use log scale for I(q) axis (default: True).",
    )
    parser.add_argument(
        "--linear-scale",
        action="store_true",
        help="Use linear scale for I(q) axis.",
    )
    parser.add_argument(
        "--kratky",
        action="store_true",
        help="Generate Kratky plot (q²·I(q) vs q).",
    )
    parser.add_argument(
        "--guinier",
        action="store_true",
        help="Generate Guinier plot (ln(I) vs q²).",
    )
    parser.add_argument(
        "--compare",
        action="store_true",
        help="Detailed comparison of 2-5 profiles with residuals.",
    )
    parser.add_argument(
        "--rg",
        action="store_true",
        help="Estimate and print Rg values from Guinier analysis.",
    )
    parser.add_argument(
        "--max-traces",
        type=int,
        default=50,
        help="Maximum number of traces to plot (default: 50).",
    )
    parser.add_argument(
        "--title",
        type=str,
        default=None,
        help="Plot title.",
    )

    args = parser.parse_args()

    if not HAS_MATPLOTLIB and not args.rg:
        print("ERROR: matplotlib is required for plotting.")
        print("Install it: pip install matplotlib")
        sys.exit(1)

    dat_files = collect_dat_files(args.paths, args.select)
    if not dat_files:
        print("No .dat files found.")
        sys.exit(1)

    print(f"Found {len(dat_files)} SAXS profiles")

    log_scale = not args.linear_scale

    # Determine output directory
    first_path = Path(args.paths[0])
    if first_path.is_dir():
        out_dir = first_path.parent
    else:
        out_dir = first_path.parent

    # Rg table
    if args.rg:
        print_rg_table(dat_files)
        if not (args.kratky or args.guinier or args.compare):
            return

    # Plots
    if args.compare:
        out_path = Path(args.output) if args.output else out_dir / "comparison.png"
        title = args.title or "SAXS Profile Comparison"
        plot_comparison(dat_files, out_path, title=title)

    elif args.kratky:
        out_path = Path(args.output) if args.output else out_dir / "kratky.png"
        title = args.title or "Kratky Plot"
        plot_kratky(dat_files, out_path, title=title, max_traces=args.max_traces)

    elif args.guinier:
        out_path = Path(args.output) if args.output else out_dir / "guinier.png"
        title = args.title or "Guinier Plot"
        plot_guinier(dat_files, out_path, title=title, max_traces=args.max_traces)

    else:
        out_path = Path(args.output) if args.output else out_dir / "saxs_overlay.png"
        title = args.title or "SAXS Profiles"
        plot_overlay(dat_files, out_path, log_scale=log_scale,
                     title=title, max_traces=args.max_traces)


if __name__ == "__main__":
    main()
