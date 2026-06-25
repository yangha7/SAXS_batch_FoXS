# SAXS_Batch

Batch computation of Small-Angle X-ray Scattering (SAXS) profiles for PDB
structures using [FoXS](https://modbase.compbio.ucsf.edu/foxs/) (Fast SAXS
Profile Computation with Debye Formula).

## Overview

This tool runs FoXS on all PDB files in a directory in parallel and aggregates
the results into a summary CSV.

## Requirements

- Python 3.11+
- IMP/FoXS installed via conda:
  ```bash
  conda activate foxs_env
  ```

## Usage

```bash
# Compute SAXS profiles for all PDB files in a directory
python run_batch.py --pdb-dir /path/to/pdb/files

# Control FoXS parameters
python run_batch.py --pdb-dir /path/to/pdb/files --max-q 0.5 --num-points 500

# Parallel execution (default: number of CPU cores)
python run_batch.py --pdb-dir /path/to/pdb/files --workers 8

# Plot all profiles overlaid
python plot_profiles.py output/
```

## Output Structure

```
output/
├── saxs/
│   ├── structure1.pdb.dat
│   ├── structure2.pdb.dat
│   └── ...
└── summary.csv
```

## FoXS Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--max-q` | 0.5 | Maximum q value (Å⁻¹) |
| `--num-points` | 500 | Number of points in the profile |
| `--hydrogens` | false | Explicitly consider hydrogens |

## References

- Schneidman-Duhovny D, Hammel M, Sali A. *FoXS: a web server for rapid
  computation and fitting of SAXS profiles.* Nucleic Acids Research, 2010.
- IMP: Integrative Modeling Platform. https://integrativemodeling.org/
