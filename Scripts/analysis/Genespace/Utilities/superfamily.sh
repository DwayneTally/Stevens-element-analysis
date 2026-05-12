#!/usr/bin/env bash
set -euo pipefail

ROOT="${1:-$PWD}"                
OUT="${2:-genome_superfamily.tsv}"

tmp="$(mktemp)"
trap 'rm -f "$tmp"' EXIT

echo -e "genome\tsuperfamily" > "$OUT"

#Loop over genespace_* dirs
for d in "$ROOT"/genespace_run_*; do
  [[ -d "$d" ]] || continue

  base="$(basename "$d")"

  #Extract superfamily name
  sf="${base#genespace_run_}"
  sf="${sf%_dropIDs}"

  #Pull genome IDs from BED files
  if [[ -d "$d/bed" ]] && compgen -G "$d/bed/*.bed" > /dev/null; then
    for f in "$d"/bed/*.bed; do
      b="$(basename "$f")"
      g="${b%.bed}"
      echo -e "${g}\t${sf}" >> "$tmp"
    done
  else
    echo "[warn] No bed/*.bed found in: $d" >&2
  fi

done

#write final table
sort -u "$tmp" >> "$OUT"

echo "[done] Wrote $OUT"
echo "[info] Unique genomes: $(tail -n +2 "$OUT" | cut -f1 | sort -u | wc -l)"
echo "[info] Superfamilies:  $(tail -n +2 "$OUT" | cut -f2 | sort -u | wc -l)"

#Report if any genome appears in >1 superfamily
awk 'NR==1{next} {a[$1]= (a[$1]=="" ? $2 : a[$1]","$2)} END{for(g in a) if(a[g]~ /,/) print g"\t"a[g]}' "$OUT" \
  > superfamily_conflicts.tsv

if [[ -s superfamily_conflicts.tsv ]]; then
  echo "[warn] Some genomes appear in multiple superfamilies -> superfamily_conflicts.tsv"
else
  rm -f superfamily_conflicts.tsv
  echo "[ok] No superfamily conflicts detected"
fi

