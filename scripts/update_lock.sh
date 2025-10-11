#!/usr/bin/env bash
set -euo pipefail

# Regenerate requirements.lock from requirements.txt for reproducible environments.
echo "[MolBridge] Regenerating requirements.lock from requirements.txt"
pip install --upgrade pip >/dev/null 2>&1 || true
pip install -r requirements.txt
pip freeze --exclude-editable > requirements.lock
echo "[MolBridge] Done. Review diff and commit if intentional."