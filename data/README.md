# Input data

The input data are **not** stored in this repository. Three public GEO datasets and three
third-party reference resources are used, plus a small set of derived tables.

See **[../docs/DATA.md](../docs/DATA.md)** for the full list, the expected directory layout, and how
to obtain each item.

Quick summary of what to place here:

```
data/
├── raw/
│   ├── GSE32472/                  # microarray series matrix + derived expression matrices
│   │   └── GSE25286_mouse_validation/   # mouse hyperoxia validation set (GSE25286)
│   ├── GSE220135/                 # RNA-seq counts/expression, day 14 and day 28
│   ├── E-MTAB-14507/              # optional additional neonatal cohort
│   └── GSE32472.xlsx              # severity grades at 28 days
└── derived/                       # training/validation tables extracted from the above
```
