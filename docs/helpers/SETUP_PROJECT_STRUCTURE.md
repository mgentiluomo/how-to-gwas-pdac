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
# Create base directories
mkdir -p \
    scripts \
    scripts/dev \
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
│   │   ├── init_project.sh
│   │   ├── section_manifest.txt
│   │   ├── script_manifest.txt
│   │   └── tool_manifest.tsv
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
│   ├── bin/                         # Command links used by the tutorial
│   │   ├── plink
│   │   ├── plink2
│   │   ├── metal
│   │   └── regenie
│   ├── plink1.9/
│   ├── plink2/
│   ├── metal/
│   ├── micromamba/
│   ├── micromamba-root/
│   └── regenie/
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

The current script-copy step is manifest-based. It copies every section that has a `scripts/` folder, so future sections can be added without changing the setup test:

```bash
git clone https://github.com/mgentiluomo/how-to-gwas-pdac.git

mkdir -p scripts/dev
find how-to-gwas-pdac/scripts/dev -maxdepth 1 -type f -exec cp {} scripts/dev/ \;

: > scripts/dev/section_manifest.txt
: > scripts/dev/script_manifest.txt

for section_dir in how-to-gwas-pdac/sections/*/; do
  section=$(basename "$section_dir")
  if [ -d "$section_dir/scripts" ]; then
    mkdir -p "scripts/$section"
    echo "scripts/$section" >> scripts/dev/section_manifest.txt
    find "$section_dir/scripts" -maxdepth 1 -type f -exec cp {} "scripts/$section/" \;
    find "scripts/$section" -maxdepth 1 -type f | sort >> scripts/dev/script_manifest.txt
  fi
done

find scripts/dev -maxdepth 1 -type f -name "*.sh" | sort >> scripts/dev/script_manifest.txt
sort -u scripts/dev/section_manifest.txt -o scripts/dev/section_manifest.txt
sort -u scripts/dev/script_manifest.txt -o scripts/dev/script_manifest.txt
rm -r how-to-gwas-pdac
```

The **automated approach** uses:
```bash
bash scripts/dev/init_project.sh    # Creates folder structure
bash scripts/dev/download_demo_data.sh  # Downloads data + verifies
bash scripts/dev/tools_setup.sh     # Installs tools and writes tool_manifest.tsv
bash scripts/dev/test.sh            # Tests everything
```

The **manual approach** is detailed in `getting_started.qmd` Step 2-6.

---

## Portable Structure Design

This structure is **fully portable** because:

1. **No hard-coded paths** — All paths are relative (e.g., `./scripts/01B_genotyping_qc/`)
2. **Manifest-based checks** — `section_manifest.txt`, `script_manifest.txt`, and `tool_manifest.tsv` tell the setup test what should exist
3. **Local tools** — Tools are installed or linked under `./tools/`, not system-wide
4. **Works anywhere** — Folder can be on Desktop, Documents, USB drive, etc.
5. **Easy backups** — Entire folder is self-contained, except for large downloaded files you may choose to regenerate

---

## Directory Descriptions

| Folder | Purpose |
|--------|---------|
| `scripts/dev/` | Utility scripts (download, setup, test) |
| `scripts/01B_genotyping_qc/` | QC pipeline scripts |
| `scripts/<section>/` | Section scripts copied from `sections/<section>/scripts/` |
| `demo_data/` | Demo dataset (7 files) |
| `tools/bin/` | Command links for PLINK1.9, PLINK2, METAL, REGENIE, and future tools |
| `tools/micromamba/` | Project-local micromamba executable |
| `tools/micromamba-root/` | Project-local environments such as `regenie_env` |
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

# Verify manifests
cat scripts/dev/section_manifest.txt
cat scripts/dev/script_manifest.txt
cat scripts/dev/tool_manifest.tsv

# Verify demo data downloaded (after Step 4)
ls -lh demo_data/ | wc -l  # Should show 7 files + header

# Verify tools installed (after Step 5)
bash scripts/dev/test.sh
```

---

## Next Steps

Once your structure is set up, follow the main guide to:
1. Download scripts (Step 3 in `getting_started.qmd`)
2. Download demo data (Step 4)
3. Install tools (Step 5)
4. Test everything (Step 6)
5. Run your first QC pipeline!
