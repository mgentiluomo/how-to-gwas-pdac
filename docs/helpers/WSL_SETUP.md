---
title: "WSL Setup for PDAC GWAS Tutorial"
author: "Contributors"
date: 2026-06-12
---

# Running the PDAC GWAS Tutorial on Windows (via WSL)

**For Windows users:** Using Windows Subsystem for Linux (WSL) is **strongly recommended** for this tutorial. The QC pipeline uses PLINK2, R scripts, and bash—tools that work best in a Linux environment.

## Why WSL?

| Aspect | Windows CMD/PowerShell | WSL (Linux) |
|--------|----------------------|-------------|
| **PLINK2 availability** | Requires manual Windows binary download | Native Linux package (`plink2`) |
| **Bash scripts** | PowerShell emulation (fragile) | Full bash compatibility |
| **File paths** | Windows-style (backslashes, drive letters) | Unix-style (forward slashes) |
| **Performance** | Slower for file I/O intensive tasks | Native Linux speed |
| **Scripting** | Complex escaping, quoting issues | Standard shell behavior |

**Recommendation:** Install WSL2 with Ubuntu, then follow this guide.

---

## Prerequisites

### 1. Enable WSL2 on Windows

Open **PowerShell as Administrator** and run:

```powershell
wsl --install
# or, if you already have WSL1:
wsl --set-default-version 2
```

Restart your computer, then verify:

```powershell
wsl --list --verbose
```

You should see Ubuntu (or your chosen distro) with VERSION = 2.

### 2. Install Ubuntu in WSL

```powershell
wsl --install -d Ubuntu
```

Or use the Microsoft Store to download Ubuntu 22.04 LTS or later.

### 3. Launch WSL Terminal

From PowerShell or Windows Terminal:

```powershell
wsl
```

Or click **Ubuntu** in Windows Terminal tabs.

---

## Installation Instructions (in WSL)

Once you're in the WSL terminal, run these commands **once** to set up the environment:

### Step 1: Update Package Manager

```bash
sudo apt-get update -y
sudo apt-get upgrade -y
```

*Note:* The `sudo` prompt may ask for your password (the one you set during WSL Ubuntu installation).

### Step 2: Install PLINK2

```bash
sudo apt-get install -y plink2
```

Verify:

```bash
plink2 --version
```

### Step 3: Install R

```bash
sudo apt-get install -y r-base r-base-dev
```

Verify:

```bash
R --version
```

### Step 4: Install Utilities

```bash
sudo apt-get install -y curl wget git
```

### Step 5: Clone or Navigate to the Repository

If you haven't cloned the repo yet:

```bash
cd /tmp
git clone https://github.com/<OWNER>/how-to-gwas-pdac.git
cd how-to-gwas-pdac
```

Or, if you cloned on Windows, access it via `/mnt/s` (for S: drive):

```bash
cd /mnt/s/Github/how-to-gwas-pdac
```

---

## Running the QC Pipeline in WSL

### Using Your gwas_tutorial Project

Once your `gwas_tutorial` folder is set up (see `getting_started.qmd`), run scripts from your local copy:

```bash
# Navigate to your project
cd ~/gwas_tutorial

# Run utility scripts from scripts/dev/
bash scripts/dev/download_demo_data.sh    # Download + verify data
bash scripts/dev/tools_setup.sh           # Install PLINK2 and PLINK1.9
bash scripts/dev/test.sh                  # Test everything

# Run QC pipeline from scripts/01B_genotyping_qc/
bash scripts/01B_genotyping_qc/01_initial_qc_stats.sh
bash scripts/01B_genotyping_qc/02_sample_callrate.sh
bash scripts/01B_genotyping_qc/03_sex_check.sh
bash scripts/01B_genotyping_qc/04_heterozygosity.sh
bash scripts/01B_genotyping_qc/05_variant_callrate.sh
bash scripts/01B_genotyping_qc/06_hardy_weinberg.sh
bash scripts/01B_genotyping_qc/07_relatedness.sh
bash scripts/01B_genotyping_qc/08_maf_filter.sh
bash scripts/01B_genotyping_qc/09_qc_summary.sh
```

Each script prints the next command to run at the end.

---

## Accessing Your Files

### Windows Path → WSL Path Mapping

| Windows | WSL |
|---------|-----|
| `C:\Users\...` | `/home/username/...` |
| `C:\` drive | `/mnt/c/` |
| `S:\` drive | `/mnt/s/` |
| `D:\` drive | `/mnt/d/` |

### Example: Accessing S: Drive Files

Your repo at `S:\Github\how-to-gwas-pdac` is accessible in WSL as:

```bash
cd /mnt/s/Github/how-to-gwas-pdac
ls -la
```

### Editing Files

You have two options:

1. **Edit in Windows (VS Code)** → Run in WSL
   - Open files in VS Code from Windows
   - Make changes
   - Run scripts in WSL terminal
   - Both environments see the same files (via `/mnt/s/`)

2. **Clone inside WSL** → Edit and run in WSL
   - Clone the repo inside WSL (`/home/...` or `/tmp/...`)
   - Edit and run entirely in WSL
   - Push commits from WSL git
   - **Faster**, but requires separate workflow

---

## Troubleshooting

### Issue: "plink2 not found"

**Solution:** Check installation:

```bash
which plink2
plink2 --version
```

If not found, reinstall:

```bash
sudo apt-get remove plink2
sudo apt-get install -y plink2
```

### Issue: "Permission denied" on script execution

**Solution:** Make the script executable:

```bash
chmod +x sections/01B_genotyping_qc/scripts/*.sh
bash sections/01B_genotyping_qc/scripts/01_initial_qc_stats.sh
```

### Issue: "File not found" when accessing Windows files

**Solution:** Use `/mnt/` prefix for Windows drives:

```bash
# ✗ Wrong
cd S:\Github\how-to-gwas-pdac

# ✓ Correct
cd /mnt/s/Github/how-to-gwas-pdac
```

### Issue: Slow performance on /mnt/s/ (Windows drive)

**Reason:** WSL filesystem performance is slower when accessing Windows drives.

**Solution:** Clone the repo inside WSL:

```bash
git clone https://github.com/<OWNER>/how-to-gwas-pdac.git ~/how-to-gwas-pdac
cd ~/how-to-gwas-pdac
```

---

## Typical Workflow

### With Your Portable gwas_tutorial Folder

This is the **recommended approach** that works anywhere:

```bash
# 1. Create your project folder (on Windows or WSL, doesn't matter)
mkdir -p ~/gwas_tutorial
cd ~/gwas_tutorial

# 2. In WSL terminal, clone the repo temporarily
git clone https://github.com/mgentiluomo/how-to-gwas-pdac.git

# 3. Copy scripts to your local project
for section_dir in how-to-gwas-pdac/sections/*/; do
  section=$(basename "$section_dir")
  mkdir -p "scripts/$section"
  cp "$section_dir/scripts"/* "scripts/$section/"
done
cp how-to-gwas-pdac/scripts/dev/* scripts/dev/
rm -rf how-to-gwas-pdac

# 4. Run utility scripts
bash scripts/dev/download_demo_data.sh
bash scripts/dev/tools_setup.sh
bash scripts/dev/test.sh

# 5. Run your QC pipeline
bash scripts/01B_genotyping_qc/01_initial_qc_stats.sh
# ... continue with other steps
```

This approach is **fully portable** — your `gwas_tutorial` folder can be:
- On your home directory (`~/gwas_tutorial`)
- On a USB drive (`/media/usb/gwas_tutorial`)
- On a network share
- Anywhere with read/write access

---

## Quick Reference: WSL Commands

```bash
# Check WSL version and distro
wsl --list --verbose

# Switch default distro
wsl --set-default <distro-name>

# Enter WSL
wsl

# Exit WSL
exit

# Run single WSL command from PowerShell
wsl bash -c "cd /mnt/s/Github/how-to-gwas-pdac && ls"

# Update WSL
wsl --update
```

---

## Next Steps

1. Follow the **Installation Instructions** above
2. Download the demo dataset: `bash demo_dataset/download_data.sh`
3. Run the first QC script: `bash sections/01B_genotyping_qc/scripts/01_initial_qc_stats.sh`
4. Check output and troubleshoot any errors

Good luck! 🚀
