#!/usr/bin/env bash
# just for my branch demo, not part of the tutorial setup
# delete after merge
set -euo pipefail

echo "==> Fetching latest Quarto release..."
QUARTO_VERSION=$(curl -sSL "https://api.github.com/repos/quarto-dev/quarto-cli/releases/latest" \
  | python3 -c "import sys, json; print(json.load(sys.stdin)['tag_name'].lstrip('v'))")
echo "    Quarto ${QUARTO_VERSION}"

echo "==> Downloading Quarto..."
curl -sSL \
  "https://github.com/quarto-dev/quarto-cli/releases/download/v${QUARTO_VERSION}/quarto-${QUARTO_VERSION}-linux-amd64.tar.gz" \
  -o /tmp/quarto.tar.gz

echo "==> Extracting Quarto..."
mkdir -p /tmp/quarto
tar -xzf /tmp/quarto.tar.gz -C /tmp/quarto --strip-components=1
export PATH="/tmp/quarto/bin:$PATH"

echo "==> Quarto $(quarto --version) ready"

echo "==> Rendering site..."
quarto render

echo "==> Done — output in _site/"