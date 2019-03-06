#!/usr/bin/env bash

# Author: Ken Youens-Clark <kyclark@email.arizona.edu>

set -u

CWD=$(cd $(dirname "$0") && pwd)
cd "$CWD"

FASTA="190102_ANT28-5_all_eASV_seqs_derep.fasta"
DB_DIR="blast"

[[ -f "$FASTA.gz" ]] && gunzip "$FASTA.gz"

if [[ -f "$FASTA" ]]; then
    echo "Making BLAST DB from FASTA \"$FASTA\""
    makeblastdb -in "$FASTA" -out "blast" -input_type fasta -dbtype nucl
    gzip "$FASTA"
    echo "Done."
else
    echo "Missing FASTA \"$FASTA\""
fi
