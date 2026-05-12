#!/bin/bash
COMPARISON="Elements_vis_ALG/data/stevens.tsv"
OUTDIR="output"
mkdir -p "$OUTDIR"

get_busco_path () {
  BUSCO_SPECIES="$1"
  BUSCO_BASE="../busco_beetle_result"
  BUSCO_PAT_GCA="${BUSCO_BASE}/busco_${BUSCO_SPECIES}_GCA_*/run_*/full_table.tsv"
  BUSCO_PAT_GCF="${BUSCO_BASE}/busco_${BUSCO_SPECIES}_GCF_*/run_*/full_table.tsv"

  BUSCO_MATCH=$(compgen -G "$BUSCO_PAT_GCA" | head -n1 || true)
  if [[ -n "${BUSCO_MATCH:-}" && -f "$BUSCO_MATCH" ]]; then
    echo "$BUSCO_MATCH"
    return 0
  fi
  BUSCO_MATCH=$(compgen -G "$BUSCO_PAT_GCF" | head -n1 || true)
  if [[ -n "${BUSCO_MATCH:-}" && -f "$BUSCO_MATCH" ]]; then
    echo "$BUSCO_MATCH"
    return 0
  fi
  return 1
}

SPECIES_LIST=(
  Liopterus_haemorrhoidalis
  Carabus_problematicus
  Nebria_salina
  Nebria_brevicollis
  Clivina_fossor
  Agonum_fuliginosum
  Dascilus_cervinus
  Agrilus_canescens
  Hololepta_plana
  Stenus_bimaculatus
  Philonthus_cognatus
  Dorcus_parallelipipedus
  Dorcus_hopei
  Xestobium_rufovillosum
  Tribolium_confusum
  Diaperis_boleti
  Aethina_tumida
  Platypus_cylindrus
  Orchestes_rusci
  Ceutorhynchus_assimilis
  Arhopalus_rusticus
  Octodonta_nipae
  Cryptocephalus_moraei
  Cryptocephalus_primarius
  Gastrophysa_polygoni
  Longitarsus_flavicornis
  Neocrepidodera_transversa
  Hermaeophaga_mercurialis
  Agelastica_alni
  Ophraella_communa
  Lochmaea_crataegi
  Galerucella_nymphaeae
  Galeruca_laticollis
  Diorhabda_carinulata
  Diorhabda_elongata
  Diorhabda_carinata
  Diorhabda_sublineata
)

for SPECIES in "${SPECIES_LIST[@]}"; do
  echo "=== ${SPECIES} ==="
  BUSCO_PATH=""
  if BUSCO_PATH=$(get_busco_path "$SPECIES"); then
    echo "Using BUSCO: $BUSCO_PATH"
  else
    echo "BUSCO file not found for ${SPECIES} Skipping..."
    continue
  fi

  SHORT="$(echo "$SPECIES" | awk -F'_' '{printf "%s.%s", substr($1,1,1), $2}')"

  Rscript element_vis.R \
    -b "$BUSCO_PATH" \
    -c "$COMPARISON" \
    -s "$SHORT" \
    -o "${OUTDIR}/${SHORT}.pdf" \
    --stacked \
    --stackedOutPlot "${OUTDIR}/${SHORT}_stacked.pdf" \
    --threshold 0.1
done