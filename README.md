# SAXS_Batch

Batch computation of Small-Angle X-ray Scattering (SAXS) profiles for DNA
oligomers using [FoXS](https://modbase.compbio.ucsf.edu/foxs/) (Fast SAXS
Profile Computation with Debye Formula).

## Overview

This tool generates all possible DNA oligomer sequences of a given length,
builds their 3D structures using
[DNA_Builder](../DNA_Builder), computes theoretical SAXS profiles via FoXS,
and aggregates the results for analysis.

## Requirements

- Python 3.11+
- [DNA_Builder](../DNA_Builder) (sibling directory)
- IMP/FoXS installed via conda:
  ```bash
  conda activate foxs_env
  ```

## Usage

```bash
# Compute SAXS profiles for all 4-mer B-DNA sequences
python run_batch.py --length 4 --form B

# Compute for all 6-mers in A, B, and Z forms
python run_batch.py --length 6 --form A B Z

# Compute for specific sequences
python run_batch.py --sequences ATCGATCG GCGCGCGC --form B

# Use existing PDB files from a directory
python run_batch.py --pdb-dir /path/to/pdb/files

# Control FoXS parameters
python run_batch.py --length 4 --form B --max-q 0.5 --num-points 500

# Parallel execution (default: number of CPU cores)
python run_batch.py --length 4 --form B --workers 8

# Plot all profiles overlaid
python plot_profiles.py output/B_form/
```

## Output Structure

```
output/
├── B_form/
│   ├── pdb/
│   │   ├── B_AAAA.pdb
│   │   ├── B_AAAC.pdb
│   │   └── ...
│   ├── saxs/
│   │   ├── B_AAAA.pdb.dat
│   │   ├── B_AAAC.pdb.dat
│   │   └── ...
│   └── summary.csv
├── A_form/
│   └── ...
└── Z_form/
    └── ...
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
