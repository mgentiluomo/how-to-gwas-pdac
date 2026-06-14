---
title: "WSL Setup for PDAC GWAS Tutorial"
author: "Contributors"
date: 2026-06-12
---

# Running the PDAC GWAS Tutorial on Windows (via WSL)

**For Windows users:** Using Windows Subsystem for Linux (WSL) is **strongly recommended** for this tutorial. The pipeline uses bash scripts and Linux/macOS command-line tools, which work most consistently inside WSL.

## Why WSL?

| Aspect | Windows CMD/PowerShell | WSL (Linux) |
|--------|----------------------|-------------|
| **Tool setup** | Requires separate Windows binaries and PATH edits | Tutorial scripts install tools into the project folder |
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

The short video below shows how to open Ubuntu from Windows Terminal:

<video controls width="100%">
  <source src="illustration/wsl_ubuntu.mp4" type="video/mp4">
  Your browser does not support the video tag.
</video>

---

## Basic WSL Preparation

Once you're in the WSL terminal, run these commands **once** to prepare Ubuntu. These commands install only the basic utilities needed to download and run the tutorial setup scripts. PLINK2, PLINK1.9, METAL, micromamba, and REGENIE are installed later by `bash scripts/dev/tools_setup.sh`.

### Step 1: Update Package Manager

```bash
sudo apt-get update -y
sudo apt-get upgrade -y
```

*Note:* The `sudo` prompt may ask for your password (the one you set during WSL Ubuntu installation).

### Step 2: Install Basic Utilities

```bash
sudo apt-get install -y curl wget git unzip tar bzip2
```

Install R, which is used for QC plots and checked by the setup test:

```bash
sudo apt-get install -y r-base r-base-dev
```

Verify the basic utilities:

```bash
git --version
curl --version
wget --version
```

Verify R too:

```bash
R --version
```


---

## Accessing Your Files

### Windows Path → WSL Path Mapping

| Windows | WSL |
|---------|-----|
| `C:\Users\...` | `/home/username/...` |
| `C:\` drive | `/mnt/c/` |
| `S:\` drive | `/mnt/s/` |
| `D:\` drive | `/mnt/d/` |

### Example: Accessing S: Drive and Example folder

Your repo at `S:\Example` is accessible in WSL as:

```bash
cd /mnt/s/Example
ls -la
```
You can return back to your HOME directory with:

```bash
cd $HOME
ls -la
```


## Troubleshooting

### Issue: WSL cannot connect to the internet on a managed network

**Reason:** Some universities, hospitals, and institutes use a network proxy. Windows may already know this proxy, but WSL may not use the Windows setting automatically.

Ask your local IT team for the proxy address and port. It usually looks like:

```bash
http://proxy.example.org:3128/
```

Replace the example proxy below with your institution's real proxy:

```bash
PROXY="http://proxy.example.org:3128/"
NO_PROXY="localhost,127.0.0.1"
```

Set the proxy for the current WSL terminal:

```bash
export http_proxy="$PROXY"
export https_proxy="$PROXY"
export ftp_proxy="$PROXY"
export HTTP_PROXY="$PROXY"
export HTTPS_PROXY="$PROXY"
export FTP_PROXY="$PROXY"
export no_proxy="$NO_PROXY"
export NO_PROXY="$NO_PROXY"
```

Test the connection:

```bash
curl -I https://github.com
wget --spider https://github.com
```

If this works, make it permanent for future WSL terminals:

```bash
cat >> "$HOME/.bashrc" <<'EOF'

# Proxy for WSL network access.
# Replace this example with your institution's proxy.
export http_proxy="http://proxy.example.org:3128/"
export https_proxy="http://proxy.example.org:3128/"
export ftp_proxy="http://proxy.example.org:3128/"
export HTTP_PROXY="http://proxy.example.org:3128/"
export HTTPS_PROXY="http://proxy.example.org:3128/"
export FTP_PROXY="http://proxy.example.org:3128/"
export no_proxy="localhost,127.0.0.1"
export NO_PROXY="localhost,127.0.0.1"
EOF

source "$HOME/.bashrc"
```

`apt-get` may need its own proxy setting because it runs with `sudo`:

```bash
sudo tee /etc/apt/apt.conf.d/95proxy >/dev/null <<EOF
Acquire::http::Proxy "$PROXY";
Acquire::https::Proxy "$PROXY";
Acquire::ftp::Proxy "$PROXY";
EOF

sudo apt-get update
```

Git can also be configured separately:

```bash
git config --global http.proxy "$PROXY"
git config --global https.proxy "$PROXY"
git ls-remote https://github.com/git/git.git HEAD
```

If your proxy requires a username or password, ask your IT team for the recommended WSL setup. Avoid putting passwords into shared scripts or tutorial files.

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
git clone https://github.com/mgentiluomo/how-to-gwas-pdac.git ~/how-to-gwas-pdac
cd ~/how-to-gwas-pdac
```

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

1. Follow [Before you start: setup and your first QC pipeline](../../getting_started.qmd)
2. Return to this page if you need help with WSL paths or basic WSL commands.

Good luck!
