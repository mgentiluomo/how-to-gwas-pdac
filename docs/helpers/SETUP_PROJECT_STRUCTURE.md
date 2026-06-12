# Setting Up Your PDAC GWAS Project Structure

This guide helps you set up a flexible, portable project folder that works on **any system** without hard-coded paths.

---

## Step 1: Create Your Project Folder

Choose any location on your computer and create a folder named `gwas_tutorial`:

**Windows (bash in WSL2):**
```bash
mkdir -p ~/gwas_tutorial
cd ~/gwas_tutorial
```

**macOS/Linux (bash):**
```bash
mkdir -p ~/gwas_tutorial
cd ~/gwas_tutorial
```

> **Note:** You can place this folder anywhere. The structure will work from any location because we use relative paths.

---

## Step 2: Create the Directory Structure

```bash
# Create all directories
mkdir -p \
    scripts/dev \
    scripts/01B_genotyping_qc \
    scripts/02_population_stratification \
    scripts/03_imputation \
    demo_data \
    tools/bin \
    data_processed \
    results/{qc,pop_structure,imputation,association,finemapping,meta_analysis}

echo "✓ Directory structure created"
tree -L 2  # (or: find . -type d -not -path '*/\.*' | sort)
```

---

## Final Structure

After setup, your project looks like this:

```
gwas_tutorial/
│
├── scripts/
│   ├── dev/                         # Utility scripts
│   │   ├── download_demo_data.sh
│   │   ├── tools_setup.sh
│   │   ├── test.sh
│   │   └── init_project.sh
│   │
│   ├── 01B_genotyping_qc/           # QC pipeline scripts
│   │   ├── 01_initial_qc_stats.sh
│   │   ├── 02_sample_callrate.sh
│   │   ├── ... (9 QC scripts total)
│   │   └── 09_qc_summary.sh
│   │
│   ├── 02_population_stratification/
│   ├── 03_imputation/
│   └── ...other sections...
│
├── demo_data/                       # Downloaded demo dataset files
│   ├── pdac_demo.bed
│   ├── pdac_demo.bim
│   ├── pdac_demo.fam
│   ├── phenotype.txt
│   ├── covariates.txt
│   ├── survival.txt
│   └── sample_ancestry.tsv
│
├── tools/                           # Local tool installation
│   └── bin/                         # Executables (plink, plink2, etc)
│       ├── plink
│       └── plink2
│
├── data_processed/                  # Processed data outputs
│   └── (will be populated during QC)
│
└── results/                         # Analysis results organized by workflow
    ├── qc/
    ├── pop_structure/
    ├── imputation/
    ├── association/
    ├── finemapping/
    └── meta_analysis/
```

---

## Step 3: Get the Scripts

Follow the main guide at `getting_started.qmd` which provides automated and manual options for downloading scripts in the correct structure.

The **automated approach** uses:
```bash
bash scripts/dev/init_project.sh    # Creates folder structure
bash scripts/dev/download_demo_data.sh  # Downloads data + verifies
bash scripts/dev/tools_setup.sh     # Installs PLINK2 and PLINK1.9
bash scripts/dev/test.sh            # Tests everything
```

The **manual approach** is detailed in `getting_started.qmd` Step 2-6.

---

## Portable Structure Design

This structure is **fully portable** because:

1. **No hard-coded paths** — All paths are relative (e.g., `./scripts/01B_genotyping_qc/`)
2. **Local tools** — Tools installed to `./tools/bin/`, not system-wide
3. **Works anywhere** — Folder can be on Desktop, Documents, USB drive, etc.
4. **Easy backups** — Entire folder is self-contained

---

## Directory Descriptions

| Folder | Purpose |
|--------|---------|
| `scripts/dev/` | Utility scripts (download, setup, test) |
| `scripts/01B_genotyping_qc/` | QC pipeline scripts |
| `scripts/02_population_stratification/` | Population structure analysis |
| `scripts/03_imputation/` | Genotype imputation |
| `demo_data/` | Demo dataset (7 files) |
| `tools/bin/` | PLINK, PLINK2, and other tools |
| `data_processed/` | Intermediate outputs |
| `results/` | Final analysis results |

---

## Quick Verification

To check your structure is correct:

```bash
# Verify folders exist
ls -la scripts/dev/
ls -la scripts/01B_genotyping_qc/
ls -la demo_data/
ls -la tools/bin/

# Verify demo data downloaded (after Step 4)
ls -lh demo_data/ | wc -l  # Should show 7 files + header

# Verify tools installed (after Step 5)
./tools/bin/plink2 --version
./tools/bin/plink --version
```

---

## Next Steps

Once your structure is set up, follow the main guide to:
1. Download scripts (Step 3 in `getting_started.qmd`)
2. Download demo data (Step 4)
3. Install tools (Step 5)
4. Test everything (Step 6)
5. Run your first QC pipeline!
