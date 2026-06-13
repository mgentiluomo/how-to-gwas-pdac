#!/usr/bin/env bash

################################################################################
# Tools Setup Script for PDAC GWAS Tutorial
#
# Downloads and installs tutorial-managed tools, then writes a tool manifest
# Works on Linux (WSL/native) and macOS
#
# Usage: bash scripts/dev/tools_setup.sh
#
# Creates:
#   - tools/bin/plink2 (symlink to PLINK2 binary)
#   - tools/bin/plink  (symlink to PLINK1.9 binary)
#   - tools/bin/metal  (symlink to METAL binary)
#   - tools/bin/regenie (wrapper for the micromamba REGENIE environment)
#   - scripts/dev/tool_manifest.tsv (tools checked by test.sh)
#   - Updates PATH to include tools/bin/
#
################################################################################

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
cd "$PROJECT_ROOT"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Create tools directory structure (please add more tools here as needed)
mkdir -p tools/bin
mkdir -p tools/plink2
mkdir -p tools/plink1.9
mkdir -p tools/metal
mkdir -p tools/regenie
mkdir -p tools/micromamba
mkdir -p tools/micromamba-root

TOOL_MANIFEST="scripts/dev/tool_manifest.tsv"
TOOL_MANIFEST_ROWS=()
REGENIE_ENV_NAME="${REGENIE_ENV_NAME:-regenie_env}"
MICROMAMBA_BIN="$PROJECT_ROOT/tools/micromamba/bin/micromamba"
MICROMAMBA_ROOT="$PROJECT_ROOT/tools/micromamba-root"

echo ""
echo "╔════════════════════════════════════════════════════════════════╗"
echo "║   PDAC GWAS Tutorial — Download and Install Tools              ║"
echo "╚════════════════════════════════════════════════════════════════╝"
echo ""

# ============================================================================
# Check for required utilities
# ============================================================================
check_requirements() {
    echo "Checking for required utilities..."
    echo ""

    local MISSING=0

    # Check for a download utility
    if command -v wget &> /dev/null; then
        echo -e "  ${GREEN}✓${NC} wget found"
    elif command -v curl &> /dev/null; then
        echo -e "  ${GREEN}✓${NC} curl found"
    else
        echo -e "  ${YELLOW}⚠${NC}  wget/curl not found. Attempting to install wget..."

        # Detect package manager and install
        if command -v apt-get &> /dev/null; then
            if sudo apt-get update -qq && sudo apt-get install -y wget -qq 2>/dev/null; then
                echo -e "  ${GREEN}✓${NC} wget installed via apt-get"
            else
                echo -e "  ${RED}✗${NC} Failed to install wget via apt-get (may need sudo access)"
                MISSING=1
            fi
        elif command -v yum &> /dev/null; then
            if sudo yum install -y wget -q 2>/dev/null; then
                echo -e "  ${GREEN}✓${NC} wget installed via yum"
            else
                echo -e "  ${RED}✗${NC} Failed to install wget via yum (may need sudo access)"
                MISSING=1
            fi
        elif command -v brew &> /dev/null; then
            if brew install wget -q 2>/dev/null; then
                echo -e "  ${GREEN}✓${NC} wget installed via brew"
            else
                echo -e "  ${RED}✗${NC} Failed to install wget via brew"
                MISSING=1
            fi
        else
            echo -e "  ${RED}✗${NC} No package manager found. Cannot install wget automatically."
            MISSING=1
        fi
    fi

    # Check for unzip
    if ! command -v unzip &> /dev/null; then
        echo -e "  ${YELLOW}⚠${NC}  unzip not found. Attempting to install..."

        # Detect package manager and install
        if command -v apt-get &> /dev/null; then
            if sudo apt-get install -y unzip -qq 2>/dev/null; then
                echo -e "  ${GREEN}✓${NC} unzip installed via apt-get"
            else
                echo -e "  ${RED}✗${NC} Failed to install unzip via apt-get (may need sudo access)"
                MISSING=1
            fi
        elif command -v yum &> /dev/null; then
            if sudo yum install -y unzip -q 2>/dev/null; then
                echo -e "  ${GREEN}✓${NC} unzip installed via yum"
            else
                echo -e "  ${RED}✗${NC} Failed to install unzip via yum (may need sudo access)"
                MISSING=1
            fi
        elif command -v brew &> /dev/null; then
            if brew install unzip -q 2>/dev/null; then
                echo -e "  ${GREEN}✓${NC} unzip installed via brew"
            else
                echo -e "  ${RED}✗${NC} Failed to install unzip via brew"
                MISSING=1
            fi
        else
            echo -e "  ${RED}✗${NC} No package manager found. Cannot install unzip automatically."
            MISSING=1
        fi
    else
        echo -e "  ${GREEN}✓${NC} unzip found"
    fi

    # Check for tar
    if ! command -v tar &> /dev/null; then
        echo -e "  ${YELLOW}⚠${NC}  tar not found. Attempting to install..."

        if command -v apt-get &> /dev/null; then
            if sudo apt-get install -y tar -qq 2>/dev/null; then
                echo -e "  ${GREEN}✓${NC} tar installed via apt-get"
            else
                echo -e "  ${RED}✗${NC} Failed to install tar via apt-get (may need sudo access)"
                MISSING=1
            fi
        elif command -v yum &> /dev/null; then
            if sudo yum install -y tar -q 2>/dev/null; then
                echo -e "  ${GREEN}✓${NC} tar installed via yum"
            else
                echo -e "  ${RED}✗${NC} Failed to install tar via yum (may need sudo access)"
                MISSING=1
            fi
        elif command -v brew &> /dev/null; then
            if brew install gnu-tar -q 2>/dev/null; then
                echo -e "  ${GREEN}✓${NC} tar installed via brew"
            else
                echo -e "  ${RED}✗${NC} Failed to install tar via brew"
                MISSING=1
            fi
        else
            echo -e "  ${RED}✗${NC} No package manager found. Cannot install tar automatically."
            MISSING=1
        fi
    else
        echo -e "  ${GREEN}✓${NC} tar found"
    fi

    # Check for bzip2 support used by micromamba archives
    if ! command -v bzip2 &> /dev/null; then
        echo -e "  ${YELLOW}⚠${NC}  bzip2 not found. Attempting to install..."

        if command -v apt-get &> /dev/null; then
            if sudo apt-get install -y bzip2 -qq 2>/dev/null; then
                echo -e "  ${GREEN}✓${NC} bzip2 installed via apt-get"
            else
                echo -e "  ${RED}✗${NC} Failed to install bzip2 via apt-get (may need sudo access)"
                MISSING=1
            fi
        elif command -v yum &> /dev/null; then
            if sudo yum install -y bzip2 -q 2>/dev/null; then
                echo -e "  ${GREEN}✓${NC} bzip2 installed via yum"
            else
                echo -e "  ${RED}✗${NC} Failed to install bzip2 via yum (may need sudo access)"
                MISSING=1
            fi
        elif command -v brew &> /dev/null; then
            if brew install bzip2 -q 2>/dev/null; then
                echo -e "  ${GREEN}✓${NC} bzip2 installed via brew"
            else
                echo -e "  ${RED}✗${NC} Failed to install bzip2 via brew"
                MISSING=1
            fi
        else
            echo -e "  ${RED}✗${NC} No package manager found. Cannot install bzip2 automatically."
            MISSING=1
        fi
    else
        echo -e "  ${GREEN}✓${NC} bzip2 found"
    fi

    if [[ $MISSING -eq 1 ]]; then
        echo ""
        echo -e "${RED}Required utilities are missing.${NC}"
        echo ""
        echo "If you don't have sudo access, please ask your system administrator to install:"
        echo "  • wget or curl (download tool)"
        echo "  • unzip (archive extraction)"
        echo "  • tar (archive extraction)"
        echo "  • bzip2 (micromamba archive extraction)"
        echo ""
        echo "Manual installation commands:"
        echo "  Ubuntu/Debian: sudo apt-get install -y wget unzip tar bzip2"
        echo "  CentOS/RHEL:   sudo yum install -y wget unzip tar bzip2"
        echo "  macOS:         brew install wget unzip bzip2"
        echo ""
        echo "Alternatively, you can download PLINK2, PLINK1.9, and METAL manually and place them in:"
        echo "  • tools/plink2/plink2"
        echo "  • tools/plink1.9/plink"
        echo "  • tools/metal/metal"
        echo ""
        return 1
    fi

    echo -e "  ${GREEN}✓${NC} All utilities available"
    echo ""
}

download_file() {
    local URL=$1
    local FILENAME=$2

    if command -v wget &> /dev/null; then
        wget -q -O "$FILENAME" "$URL"
    elif command -v curl &> /dev/null; then
        curl -fsSL -o "$FILENAME" "$URL"
    else
        return 127
    fi
}

register_tool() {
    local COMMAND=$1
    local LABEL=$2
    local REQUIRED=${3:-required}
    local VERSION_ARGS=${4:---version}
    local INSTALL_HINT=${5:-"Install $LABEL, then run Step 5 again"}

    TOOL_MANIFEST_ROWS+=("${COMMAND}"$'\t'"${LABEL}"$'\t'"${REQUIRED}"$'\t'"${VERSION_ARGS}"$'\t'"${INSTALL_HINT}")
}

register_default_tools() {
    TOOL_MANIFEST_ROWS=()

    # Add future tools here after their installer/check step is added.
    # Format: command, display name, required|optional, version arguments, help text.
    register_tool "plink2" "PLINK2" "required" "--version" "Run bash scripts/dev/tools_setup.sh"
    register_tool "plink" "PLINK1.9" "required" "--version" "Run bash scripts/dev/tools_setup.sh"
    register_tool "metal" "METAL" "required" "-" "Run bash scripts/dev/tools_setup.sh"
    register_tool "regenie" "REGENIE" "required" "--version" "Run bash scripts/dev/tools_setup.sh"
    register_tool "R" "R" "required" "--version" "Install R, then run bash scripts/dev/tools_setup.sh again"
}

write_tool_manifest() {
    mkdir -p "$(dirname "$TOOL_MANIFEST")"

    {
        echo "# command<TAB>label<TAB>required<TAB>version_args<TAB>install_hint"
        printf '%s\n' "${TOOL_MANIFEST_ROWS[@]}"
    } > "$TOOL_MANIFEST"

    echo -e "  ${GREEN}✓${NC} Tool manifest written: $TOOL_MANIFEST"
}

verify_tool_manifest() {
    echo ""
    echo -e "${BLUE}Step 8: Verifying tools from manifest${NC}"
    echo ""

    local TOOLS_OK=true
    local COMMAND
    local LABEL
    local REQUIRED
    local VERSION_ARGS
    local INSTALL_HINT

    if [ ! -f "$TOOL_MANIFEST" ]; then
        echo -e "  ${RED}✗${NC} $TOOL_MANIFEST not found"
        return 1
    fi

    while IFS=$'\t' read -r COMMAND LABEL REQUIRED VERSION_ARGS INSTALL_HINT; do
        COMMAND=${COMMAND%$'\r'}
        [ -z "$COMMAND" ] && continue
        [[ "$COMMAND" == \#* ]] && continue

        LABEL=${LABEL:-$COMMAND}
        REQUIRED=${REQUIRED:-required}
        VERSION_ARGS=${VERSION_ARGS:-}
        INSTALL_HINT=${INSTALL_HINT:-"Install $LABEL, then run Step 5 again"}

        if command -v "$COMMAND" &> /dev/null; then
            local VERSION_OUTPUT=""
            if [ -n "$VERSION_ARGS" ] && [ "$VERSION_ARGS" != "-" ]; then
                local VERSION_PARTS=()
                read -r -a VERSION_PARTS <<< "$VERSION_ARGS"
                VERSION_OUTPUT="$("$COMMAND" "${VERSION_PARTS[@]}" 2>&1 | head -1 || true)"
            fi

            if [ -n "$VERSION_OUTPUT" ]; then
                echo -e "  ${GREEN}✓${NC} $LABEL: $VERSION_OUTPUT"
            else
                echo -e "  ${GREEN}✓${NC} $LABEL"
            fi
        elif [ "$REQUIRED" = "optional" ]; then
            echo -e "  ${YELLOW}!${NC} $LABEL not found (optional)"
        else
            echo -e "  ${RED}✗${NC} $LABEL not found"
            echo "    $INSTALL_HINT"
            TOOLS_OK=false
        fi
    done < "$TOOL_MANIFEST"

    if [ "$TOOLS_OK" = true ]; then
        return 0
    fi

    return 1
}

# ============================================================================
# Detect OS and Architecture
# ============================================================================
detect_system() {
    local OS
    local ARCH

    OS=$(uname -s)
    ARCH=$(uname -m)

    if [[ "$OS" == "Linux" ]]; then
        echo "Linux"
        if [[ "$ARCH" == "x86_64" ]]; then
            # Check for AVX2 support
            if grep -q avx2 /proc/cpuinfo 2>/dev/null; then
                echo "x86_64_avx2"
            else
                echo "x86_64"
            fi
        elif [[ "$ARCH" == "i686" ]]; then
            echo "i686"
        else
            echo "unknown"
        fi
    elif [[ "$OS" == "Darwin" ]]; then
        echo "macOS"
        if [[ "$ARCH" == "arm64" ]]; then
            echo "arm64"
        elif [[ "$ARCH" == "x86_64" ]]; then
            # Check for AVX2 support (older Intel Macs)
            if sysctl -a | grep -q avx2; then
                echo "avx2"
            else
                echo "x86_64"
            fi
        else
            echo "unknown"
        fi
    else
        echo "unknown"
        echo "unknown"
    fi
}

# ============================================================================
# Download and extract PLINK2
# ============================================================================
install_plink2() {
    local OS=$1
    local ARCH=$2

    echo ""
    echo -e "${BLUE}Step 1: Installing PLINK2${NC}"
    echo "  OS: $OS, Architecture: $ARCH"
    echo ""

    local URL
    local FILENAME

    if [[ "$OS" == "Linux" ]]; then
        if [[ "$ARCH" == "x86_64_avx2" ]]; then
            FILENAME="plink2_linux_avx2_20260504.zip"
            URL="https://s3.amazonaws.com/plink2-assets/alpha7/plink2_linux_avx2_20260504.zip"
        elif [[ "$ARCH" == "x86_64" ]]; then
            FILENAME="plink2_linux_x86_64_20260504.zip"
            URL="https://s3.amazonaws.com/plink2-assets/alpha7/plink2_linux_x86_64_20260504.zip"
        elif [[ "$ARCH" == "i686" ]]; then
            FILENAME="plink2_linux_i686_20260504.zip"
            URL="https://s3.amazonaws.com/plink2-assets/alpha7/plink2_linux_i686_20260504.zip"
        else
            echo -e "${RED}✗ Unsupported Linux architecture: $ARCH${NC}"
            return 1
        fi
    elif [[ "$OS" == "macOS" ]]; then
        if [[ "$ARCH" == "arm64" ]]; then
            FILENAME="plink2_mac_arm64_20260504.zip"
            URL="https://s3.amazonaws.com/plink2-assets/alpha7/plink2_mac_arm64_20260504.zip"
        elif [[ "$ARCH" == "avx2" ]]; then
            FILENAME="plink2_mac_avx2_20260504.zip"
            URL="https://s3.amazonaws.com/plink2-assets/alpha7/plink2_mac_avx2_20260504.zip"
        elif [[ "$ARCH" == "x86_64" ]]; then
            FILENAME="plink2_mac_20260504.zip"
            URL="https://s3.amazonaws.com/plink2-assets/alpha7/plink2_mac_20260504.zip"
        else
            echo -e "${RED}✗ Unsupported macOS architecture: $ARCH${NC}"
            return 1
        fi
    else
        echo -e "${RED}✗ Unsupported OS: $OS${NC}"
        return 1
    fi

    echo "  Downloading: $FILENAME"
    pushd tools/plink2 >/dev/null

    rm -f "$FILENAME"
    if download_file "$URL" "$FILENAME"; then
        echo -e "  ${GREEN}✓${NC} Downloaded successfully"

        # Extract
        echo "  Extracting..."
        unzip -oq "$FILENAME"
        rm "$FILENAME"

        # Check for plink2 binary
        if [[ -f "plink2" ]]; then
            chmod +x plink2
            echo -e "  ${GREEN}✓${NC} PLINK2 ready"
        else
            echo -e "  ${RED}✗${NC} plink2 binary not found after extraction"
            popd >/dev/null
            return 1
        fi
    else
        echo -e "  ${RED}✗${NC} Download failed. Check URL and internet connection."
        popd >/dev/null
        return 1
    fi

    popd >/dev/null
}

# ============================================================================
# Download and extract PLINK1.9
# ============================================================================
install_plink1_9() {
    local OS=$1
    local ARCH=$2

    echo ""
    echo -e "${BLUE}Step 2: Installing PLINK1.9${NC}"
    echo "  OS: $OS, Architecture: $ARCH"
    echo ""

    local URL
    local FILENAME

    if [[ "$OS" == "Linux" ]]; then
        if [[ "$ARCH" == *"x86_64"* ]]; then
            FILENAME="plink_linux_x86_64_20250819.zip"
            URL="https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20250819.zip"
        elif [[ "$ARCH" == "i686" ]]; then
            FILENAME="plink_linux_i686_20250819.zip"
            URL="https://s3.amazonaws.com/plink1-assets/plink_linux_i686_20250819.zip"
        else
            echo -e "${RED}✗ Unsupported Linux architecture for PLINK1.9: $ARCH${NC}"
            return 1
        fi
    elif [[ "$OS" == "macOS" ]]; then
        FILENAME="plink_mac_20250819.zip"
        URL="https://s3.amazonaws.com/plink1-assets/plink_mac_20250819.zip"
    else
        echo -e "${RED}✗ Unsupported OS: $OS${NC}"
        return 1
    fi

    echo "  Downloading: $FILENAME"
    pushd tools/plink1.9 >/dev/null

    rm -f "$FILENAME"
    if download_file "$URL" "$FILENAME"; then
        echo -e "  ${GREEN}✓${NC} Downloaded successfully"

        # Extract
        echo "  Extracting..."
        unzip -oq "$FILENAME"
        rm "$FILENAME"

        # Check for plink binary
        if [[ -f "plink" ]]; then
            chmod +x plink
            echo -e "  ${GREEN}✓${NC} PLINK1.9 ready"
        else
            echo -e "  ${RED}✗${NC} plink binary not found after extraction"
            popd >/dev/null
            return 1
        fi
    else
        echo -e "  ${RED}✗${NC} Download failed. Check URL and internet connection."
        popd >/dev/null
        return 1
    fi

    popd >/dev/null
}

# ============================================================================
# Download and extract METAL
# ============================================================================
install_metal() {
    local OS=$1

    echo ""
    echo -e "${BLUE}Step 3: Installing METAL${NC}"
    echo "  OS: $OS"
    echo ""

    local URL
    local FILENAME

    if [[ "$OS" == "Linux" ]]; then
        FILENAME="Linux-metal.tar.gz"
        URL="https://csg.sph.umich.edu/abecasis/Metal/download/Linux-metal.tar.gz"
    elif [[ "$OS" == "macOS" ]]; then
        FILENAME="Darwin-metal.tar.gz"
        URL="https://csg.sph.umich.edu/abecasis/Metal/download/Darwin-metal.tar.gz"
    else
        echo -e "${RED}✗ Unsupported OS for METAL: $OS${NC}"
        return 1
    fi

    echo "  Downloading: $FILENAME"
    pushd tools/metal >/dev/null

    find . -mindepth 1 -maxdepth 1 -exec rm -rf {} +
    if download_file "$URL" "$FILENAME"; then
        echo -e "  ${GREEN}✓${NC} Downloaded successfully"

        echo "  Extracting..."
        tar -xzf "$FILENAME"
        rm "$FILENAME"

        local METAL_BINARY
        METAL_BINARY="$(find . -type f | while IFS= read -r file; do
            case "$(basename "$file")" in
                metal|METAL)
                    printf '%s\n' "$file"
                    break
                    ;;
            esac
        done)"

        if [ -n "$METAL_BINARY" ]; then
            chmod +x "$METAL_BINARY"
            METAL_BINARY=${METAL_BINARY#./}
            echo "$METAL_BINARY" > .metal_binary_path
            echo -e "  ${GREEN}✓${NC} METAL ready"
        else
            echo -e "  ${RED}✗${NC} metal binary not found after extraction"
            popd >/dev/null
            return 1
        fi
    else
        echo -e "  ${RED}✗${NC} Download failed. Check URL and internet connection."
        popd >/dev/null
        return 1
    fi

    popd >/dev/null
}

# ============================================================================
# Download and install micromamba
# ============================================================================
detect_micromamba_platform() {
    local OS
    local ARCH

    OS="$(uname -s)"
    ARCH="$(uname -m)"

    if [[ "$OS" == "Linux" ]]; then
        if [[ "$ARCH" == "x86_64" ]]; then
            echo "linux-64"
        elif [[ "$ARCH" == "aarch64" || "$ARCH" == "arm64" ]]; then
            echo "linux-aarch64"
        else
            echo "unknown"
        fi
    elif [[ "$OS" == "Darwin" ]]; then
        if [[ "$ARCH" == "x86_64" ]]; then
            echo "osx-64"
        elif [[ "$ARCH" == "arm64" ]]; then
            echo "osx-arm64"
        else
            echo "unknown"
        fi
    else
        echo "unknown"
    fi
}

install_micromamba() {
    echo ""
    echo -e "${BLUE}Step 4: Installing micromamba${NC}"
    echo "  Location: tools/micromamba/"
    echo ""

    if [ -x "$MICROMAMBA_BIN" ]; then
        MICROMAMBA_VERSION="$("$MICROMAMBA_BIN" --version 2>/dev/null || true)"
        echo -e "  ${GREEN}✓${NC} micromamba already installed (${MICROMAMBA_VERSION:-version unknown})"
        return 0
    fi

    local PLATFORM
    PLATFORM="$(detect_micromamba_platform)"

    if [[ "$PLATFORM" == "unknown" ]]; then
        echo -e "  ${RED}✗${NC} Unsupported platform for micromamba: $(uname -s) $(uname -m)"
        return 1
    fi

    local URL="https://micro.mamba.pm/api/micromamba/${PLATFORM}/latest"
    local FILENAME="micromamba-${PLATFORM}.tar.bz2"

    echo "  Platform: $PLATFORM"
    echo "  Downloading: $FILENAME"

    pushd tools/micromamba >/dev/null

    rm -rf bin info "$FILENAME"
    if download_file "$URL" "$FILENAME"; then
        tar -xjf "$FILENAME" bin/micromamba
        rm "$FILENAME"
        chmod +x bin/micromamba
        echo -e "  ${GREEN}✓${NC} micromamba ready"
    else
        echo -e "  ${RED}✗${NC} micromamba download failed"
        popd >/dev/null
        return 1
    fi

    popd >/dev/null
}

# ============================================================================
# Install REGENIE with project-local micromamba
# ============================================================================
install_regenie() {
    echo ""
    echo -e "${BLUE}Step 5: Installing REGENIE with micromamba${NC}"
    echo "  Environment: $REGENIE_ENV_NAME"
    echo "  Environment root: tools/micromamba-root/"
    echo ""

    export MAMBA_ROOT_PREFIX="$MICROMAMBA_ROOT"

    if "$MICROMAMBA_BIN" env list | awk '{print $1}' | grep -Fxq "$REGENIE_ENV_NAME"; then
        echo "  Environment exists. Ensuring REGENIE is installed..."
        "$MICROMAMBA_BIN" install -y -n "$REGENIE_ENV_NAME" -c conda-forge -c bioconda regenie
    else
        echo "  Creating $REGENIE_ENV_NAME..."
        "$MICROMAMBA_BIN" create -y -n "$REGENIE_ENV_NAME" -c conda-forge -c bioconda regenie
    fi

    mkdir -p tools/regenie
    cat > tools/regenie/regenie <<EOF
#!/usr/bin/env bash
SCRIPT_DIR="\$(cd "\$(dirname "\${BASH_SOURCE[0]}")" && pwd)"
export MAMBA_ROOT_PREFIX="\$(cd "\$SCRIPT_DIR/../micromamba-root" && pwd)"
exec "\$SCRIPT_DIR/../micromamba/bin/micromamba" run -n "$REGENIE_ENV_NAME" regenie "\$@"
EOF
    chmod +x tools/regenie/regenie

    if tools/regenie/regenie --version >/dev/null 2>&1; then
        echo -e "  ${GREEN}✓${NC} REGENIE ready"
    else
        echo -e "  ${RED}✗${NC} REGENIE was installed but could not be run"
        return 1
    fi
}

# ============================================================================
# Create symbolic links
# ============================================================================
create_symlinks() {
    echo ""
    echo -e "${BLUE}Step 6: Creating symbolic links${NC}"

    # Remove old symlinks if they exist
    rm -f tools/bin/plink2 tools/bin/plink tools/bin/metal tools/bin/regenie

    # Create symlinks
    ln -s ../plink2/plink2 tools/bin/plink2
    ln -s ../plink1.9/plink tools/bin/plink
    if [ -f tools/metal/.metal_binary_path ]; then
        local METAL_BINARY
        METAL_BINARY="$(cat tools/metal/.metal_binary_path)"
        ln -s "../metal/$METAL_BINARY" tools/bin/metal
    else
        echo -e "  ${RED}✗${NC} METAL binary path not found"
        return 1
    fi
    ln -s ../regenie/regenie tools/bin/regenie

    echo -e "  ${GREEN}✓${NC} Symlinks created in tools/bin/"
}

# ============================================================================
# Add to PATH
# ============================================================================
add_to_path() {
    echo ""
    echo -e "${BLUE}Step 7: Configuring PATH${NC}"

    local TOOLS_BIN="$(cd tools/bin && pwd)"
    local PROFILE_FILE="$HOME/.bashrc"
    local PATH_LINE="export PATH=\"$TOOLS_BIN:\$PATH\""

    if [[ "$(uname -s)" == "Darwin" && "${SHELL:-}" == */zsh ]]; then
        PROFILE_FILE="$HOME/.zshrc"
    fi

    # Check if already in PATH
    if [[ ":$PATH:" == *":$TOOLS_BIN:"* ]]; then
        echo -e "  ${YELLOW}⚠${NC}  Already in PATH"
    else
        # Add to current session
        export PATH="$TOOLS_BIN:$PATH"
        echo -e "  ${GREEN}✓${NC} Added to current session PATH"

        # Add to shell profile for persistence
        if ! grep -Fq "$TOOLS_BIN" "$PROFILE_FILE" 2>/dev/null; then
            echo "$PATH_LINE" >> "$PROFILE_FILE"
            echo -e "  ${GREEN}✓${NC} Added to $PROFILE_FILE for next session"
        fi
    fi
}

# ============================================================================
# Main execution
# ============================================================================
main() {
    echo "Checking system requirements..."
    check_requirements || exit 1

    echo "Detecting system..."
    echo ""

    local OS
    local ARCH
    local DETECTED_SYSTEM
    DETECTED_SYSTEM="$(detect_system)"
    OS="$(printf '%s\n' "$DETECTED_SYSTEM" | sed -n '1p')"
    ARCH="$(printf '%s\n' "$DETECTED_SYSTEM" | sed -n '2p')"

    if [[ "$ARCH" == "unknown" ]]; then
        echo -e "${RED}✗ Could not detect your system architecture${NC}"
        echo "Detected: OS=$OS, Architecture=$ARCH"
        return 1
    fi

    echo "Detected: OS=$OS, Architecture=$ARCH"

    # Install tools
    install_plink2 "$OS" "$ARCH" || exit 1
    install_plink1_9 "$OS" "$ARCH" || exit 1
    install_metal "$OS" || exit 1
    install_micromamba || exit 1
    install_regenie || exit 1

    # Setup
    create_symlinks
    add_to_path
    register_default_tools
    write_tool_manifest
    verify_tool_manifest || exit 1

    echo ""
    echo "╔════════════════════════════════════════════════════════════════╗"
    echo "║                    Setup Complete!                            ║"
    echo "╚════════════════════════════════════════════════════════════════╝"
    echo ""
    echo "Required tools are listed in:"
    echo "  $TOOL_MANIFEST"
    echo ""
    echo "Next: Test your setup"
    echo "  bash scripts/dev/test.sh"
    echo ""
}

main "$@"
