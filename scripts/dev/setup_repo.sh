#!/bin/bash
# Initialise the local repo with YOUR git identity and push it to GitHub.
#
# Prerequisites:
#   - git installed
#   - either the GitHub CLI `gh` (recommended), or an empty repo already
#     created on github.com to push into.
#
# Usage:
#   1. Edit the three variables below.
#   2. bash setup_repo.sh
set -euo pipefail

GIT_NAME="Manuel Gentiluomo"
GIT_EMAIL="manuel.gentiluomo@unipi.it"
# Owner can be your user (mgentiluomo) or a UniPi org, e.g. unipi-genetics
REPO_OWNER="mgentiluomo"
REPO_NAME="how-to-gwas-pdac"
VISIBILITY="public"        # or: private

git init -b main
git config user.name  "${GIT_NAME}"
git config user.email "${GIT_EMAIL}"
git add .
git commit -m "Initial scaffold: README, CONTRIBUTING, demo dataset pipeline, section folders"

if command -v gh >/dev/null 2>&1; then
    gh repo create "${REPO_OWNER}/${REPO_NAME}" --"${VISIBILITY}" \
        --source=. --remote=origin --push
    echo "Repository created and pushed: https://github.com/${REPO_OWNER}/${REPO_NAME}"
else
    echo "GitHub CLI not found."
    echo "1) Create an EMPTY repo at https://github.com/new"
    echo "   (name: ${REPO_NAME}; do NOT add README/LICENSE/.gitignore)."
    echo "2) Then run:"
    echo "   git remote add origin https://github.com/${REPO_OWNER}/${REPO_NAME}.git"
    echo "   git push -u origin main"
fi
