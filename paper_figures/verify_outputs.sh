#!/usr/bin/env bash
# Verify every generated panel is genuinely print-ready.
#
#   bash verify_outputs.sh
#
# Checks each PDF for: (1) zero embedded rasters - i.e. real vector, not a
# screenshot in a PDF wrapper; (2) no Type 3 fonts - journal production systems
# routinely reject them; (3) fonts actually embedded. Also checks SVGs for
# embedded bitmaps. Exits non-zero if anything fails.
#
# Needs poppler-utils (pdfimages, pdffonts, pdfinfo):  apt-get install poppler-utils
set -uo pipefail
cd "$(dirname "$0")"

for tool in pdfimages pdffonts pdfinfo; do
  command -v $tool >/dev/null || { echo "MISSING: $tool (apt-get install poppler-utils)"; exit 2; }
done

fail=0
printf "%-56s %8s %-22s %s\n" "FILE" "RASTERS" "FONT TYPE" "VERDICT"
printf '%.0s-' {1..110}; echo

while IFS= read -r f; do
  # _audit_reference holds a deliberate browser capture; not a deliverable.
  case "$f" in *_audit_reference*) continue;; esac
  rasters=$(pdfimages -list "$f" 2>/dev/null | tail -n +3 | grep -c . || true)
  ftypes=$(pdffonts "$f" 2>/dev/null | tail -n +3 | awk 'NF{print $2" "$3}' | sort -u | tr '\n' ',' | sed 's/,$//')
  notemb=$(pdffonts "$f" 2>/dev/null | tail -n +3 | awk 'NF && $(NF-3)=="no"' | grep -c . || true)

  verdict="OK"
  [ "$rasters" != "0" ] && { verdict="FAIL: contains raster images"; fail=1; }
  case "$ftypes" in *"Type 3"*) verdict="FAIL: Type 3 font"; fail=1;; esac
  [ "$notemb" != "0" ] && { verdict="FAIL: font not embedded"; fail=1; }

  printf "%-56s %8s %-22s %s\n" "$(basename "$f")" "$rasters" "${ftypes:-none}" "$verdict"
done < <(find . -name '*.pdf' | sort)

echo
while IFS= read -r f; do
  n=$(grep -c 'data:image' "$f" 2>/dev/null || true)
  [ "$n" != "0" ] && { echo "FAIL: $(basename "$f") embeds $n bitmap(s)"; fail=1; }
done < <(find . -name '*.svg' | sort)

echo
if [ "$fail" -eq 0 ]; then
  echo "PASS - every PDF is pure vector with embedded non-Type-3 fonts; no SVG embeds bitmaps."
else
  echo "FAILURES ABOVE - do not submit these files."
fi
exit $fail
