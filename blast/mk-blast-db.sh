#!/usr/bin/env bash

FASTA="190102_ANT28-5_all_eASV_seqs_derep.fasta"
DB_NAME="blast"

[[ -f "$FASTA.gz" ]] && gunzip "$FASTA.gz"

if [[ -f "$FASTA" ]]; then
    echo "Making BLAST DB \"$DB_NAME\" from FASTA \"$FASTA\""
    makeblastdb -in "$FASTA" -out "$DB_NAME" -input_type fasta -dbtype nucl
    gzip "$FASTA"
    echo "Done."
else
    echo "Missing FASTA \"$FASTA\""
fi
