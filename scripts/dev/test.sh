#!/usr/bin/env bash

# GWAS Tutorial Setup Test
# Tests folders, copied scripts, demo data, and required tools

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
cd "$PROJECT_ROOT"

if [ -d "tools/bin" ]; then
    export PATH="$(cd tools/bin && pwd):$PATH"
fi

echo ""
echo "=================================="
echo "GWAS Tutorial Setup Test"
echo "=================================="
echo ""

# Test folders
echo "1. Checking folder structure..."
FOLDERS=("scripts" "scripts/dev" "demo_data" "tools/bin" "data_processed" "results/qc")
FOLDER_OK=true
BASE_FOLDER_OK=true
SECTION_OK=true

for folder in "${FOLDERS[@]}"; do
    if [ -d "$folder" ]; then
        echo "   ✓ $folder"
    else
        echo "   ✗ $folder (missing)"
        BASE_FOLDER_OK=false
        FOLDER_OK=false
    fi
done

SECTION_MANIFEST="scripts/dev/section_manifest.txt"
if [ -f "$SECTION_MANIFEST" ]; then
    SECTION_COUNT=0
    while IFS= read -r folder; do
        folder=${folder%$'\r'}
        [ -z "$folder" ] && continue
        SECTION_COUNT=$((SECTION_COUNT + 1))
        if [ -d "$folder" ]; then
            echo "   ✓ $folder"
        else
            echo "   ✗ $folder (missing)"
            SECTION_OK=false
            FOLDER_OK=false
        fi
    done < "$SECTION_MANIFEST"
    if [ "$SECTION_COUNT" -eq 0 ]; then
        echo "   ✗ $SECTION_MANIFEST is empty"
        SECTION_OK=false
        FOLDER_OK=false
    fi
else
    echo "   ✗ $SECTION_MANIFEST (missing; run Step 3 again)"
    SECTION_OK=false
    FOLDER_OK=false
fi

# Test demo data files
echo ""
echo "2. Checking demo dataset files..."
DATA_FILES=("pdac_demo.bed" "pdac_demo.bim" "pdac_demo.fam" "phenotype.txt" "covariates.txt" "survival.txt" "sample_ancestry.tsv")
DATA_OK=true

for file in "${DATA_FILES[@]}"; do
    if [ -f "demo_data/$file" ]; then
        echo "   ✓ $file"
    else
        echo "   ✗ $file (missing)"
        DATA_OK=false
    fi
done

echo ""
echo "3. Checking copied scripts..."
SCRIPT_MANIFEST="scripts/dev/script_manifest.txt"
SCRIPT_OK=true

if [ -f "$SCRIPT_MANIFEST" ]; then
    SCRIPT_COUNT=0
    while IFS= read -r script; do
        script=${script%$'\r'}
        [ -z "$script" ] && continue
        SCRIPT_COUNT=$((SCRIPT_COUNT + 1))
        if [ -f "$script" ]; then
            echo "   ✓ $script"
        else
            echo "   ✗ $script (missing)"
            SCRIPT_OK=false
        fi
    done < "$SCRIPT_MANIFEST"
    if [ "$SCRIPT_COUNT" -eq 0 ]; then
        echo "   ✗ $SCRIPT_MANIFEST is empty"
        SCRIPT_OK=false
    fi
else
    echo "   ✗ $SCRIPT_MANIFEST (missing; run Step 3 again)"
    SCRIPT_OK=false
fi

# Test tools
echo ""
echo "4. Checking installed tools from manifest..."
TOOL_MANIFEST="scripts/dev/tool_manifest.tsv"
TOOLS_OK=true

if [ -f "$TOOL_MANIFEST" ]; then
    TOOL_COUNT=0
    while IFS=$'\t' read -r tool_command label required version_args install_hint check_scope; do
        tool_command=${tool_command%$'\r'}
        [ -z "$tool_command" ] && continue
        [[ "$tool_command" == \#* ]] && continue

        label=${label:-$tool_command}
        required=${required:-required}
        version_args=${version_args:-}
        install_hint=${install_hint:-"Install $label, then run Step 5 again"}
        install_hint=${install_hint%$'\r'}
        check_scope=${check_scope%$'\r'}
        case "$tool_command" in
            plink2|plink|metal|regenie)
                check_scope=${check_scope:-project}
                ;;
            *)
                check_scope=${check_scope:-path}
                ;;
        esac
        TOOL_COUNT=$((TOOL_COUNT + 1))

        TOOL_PATH=""
        if [ "$check_scope" = "project" ]; then
            TOOL_PATH="$PROJECT_ROOT/tools/bin/$tool_command"
            if [ ! -x "$TOOL_PATH" ]; then
                TOOL_PATH=""
            fi
        elif command -v "$tool_command" &> /dev/null; then
            TOOL_PATH="$(command -v "$tool_command")"
        fi

        if [ -n "$TOOL_PATH" ]; then
            VERSION_OUTPUT=""
            if [ -n "$version_args" ] && [ "$version_args" != "-" ]; then
                VERSION_PARTS=()
                read -r -a VERSION_PARTS <<< "$version_args"
                VERSION_OUTPUT="$("$TOOL_PATH" "${VERSION_PARTS[@]}" 2>&1 | head -1 || true)"
            fi

            DISPLAY_PATH="$TOOL_PATH"
            if [[ "$DISPLAY_PATH" == "$PROJECT_ROOT/"* ]]; then
                DISPLAY_PATH="${DISPLAY_PATH#"$PROJECT_ROOT"/}"
            fi

            if [ -n "$VERSION_OUTPUT" ]; then
                echo "   ✓ $label ($VERSION_OUTPUT) [$DISPLAY_PATH]"
            else
                echo "   ✓ $label [$DISPLAY_PATH]"
            fi
        elif [ "$required" = "optional" ]; then
            echo "   ! $label (optional; not found)"
        else
            echo "   ✗ $label (not found)"
            echo "     $install_hint"
            TOOLS_OK=false
        fi
    done < "$TOOL_MANIFEST"

    if [ "$TOOL_COUNT" -eq 0 ]; then
        echo "   ✗ $TOOL_MANIFEST is empty"
        TOOLS_OK=false
    fi
else
    echo "   ✗ $TOOL_MANIFEST (missing; run Step 5 again)"
    TOOLS_OK=false
fi

# Test wget
if command -v wget &> /dev/null; then
    echo "   ✓ wget"
fi

# Test git
if command -v git &> /dev/null; then
    echo "   ✓ git"
fi

# Summary
echo ""
echo "=================================="
echo "Summary"
echo "=================================="

if [ "$FOLDER_OK" = true ] && [ "$DATA_OK" = true ] && [ "$SCRIPT_OK" = true ] && [ "$TOOLS_OK" = true ]; then
    echo "✓ Initial test done! Folders, scripts, data, and required tools are working."
    echo ""
    echo "Next command to run the QC pipeline:"
    echo "  bash scripts/01B_genotyping_qc/01_initial_qc_stats.sh"
    echo ""
    exit 0
else
    echo "✗ Some issues found. Please check the output above and fix them."
    echo ""
    if [ "$BASE_FOLDER_OK" = false ]; then
        echo "  • Missing base folders: Run Step 2 again"
    fi
    if [ "$SECTION_OK" = false ]; then
        echo "  • Missing section folder list: Run Step 3 again"
    fi
    if [ "$DATA_OK" = false ]; then
        echo "  • Missing data files: Run Step 4 again"
    fi
    if [ "$SCRIPT_OK" = false ]; then
        echo "  • Missing scripts: Run Step 3 again"
    fi
    if [ "$TOOLS_OK" = false ]; then
        echo "  • Missing tools or tool manifest: Run Step 5 again"
    fi
    echo ""
    exit 1
fi
