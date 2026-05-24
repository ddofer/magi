#!/usr/bin/env bash
# Fail if git would track files larger than MAX_BYTES under clinvar_pipeline/.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
MAX_MB="${1:-50}"
MAX_BYTES=$((MAX_MB * 1024 * 1024))

cd "$ROOT"

if ! git rev-parse --git-dir >/dev/null 2>&1; then
  echo "check_large_files: no git repo at $ROOT (run from magi repo root or init git first)"
  echo "Scanning working tree instead ..."
  mapfile -t BIG < <(find . -type f ! -path './.git/*' ! -path './.venv/*' -size +"${MAX_MB}"M 2>/dev/null | sed 's|^\./||' | sort)
else
  mapfile -t BIG < <(
    git ls-files -z | while IFS= read -r -d '' f; do
      [[ -f "$f" ]] || continue
      sz=$(stat -c%s "$f" 2>/dev/null || stat -f%z "$f")
      if (( sz > MAX_BYTES )); then
        echo "$f ($(( sz / 1024 / 1024 )) MB)"
      fi
    done | sort
  )
fi

if ((${#BIG[@]})); then
  echo "ERROR: files larger than ${MAX_MB} MB would be tracked:"
  printf '  %s\n' "${BIG[@]}"
  echo "Add them to .gitignore or use Git LFS before pushing."
  exit 1
fi

echo "OK: no tracked files > ${MAX_MB} MB under $(basename "$ROOT")"
