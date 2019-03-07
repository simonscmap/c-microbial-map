#!/usr/bin/env python3
"""
Author : Ken Youens-Clark <kyclark@email.arizona.edu>
Date   : 2019-02-01
Purpose: BLAST hits to CMAP visualization
"""

import argparse
import csv
import os
import sys
import pandas as pd
import sqlite3
from opedia import db
from subprocess import getstatusoutput


# --------------------------------------------------
def get_args():
    """get command-line arguments"""
    parser = argparse.ArgumentParser(
        description='BLAST hits to CMAP visualization',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument(
        '-q',
        '--query',
        help='Query file for BLAST',
        metavar='str',
        type=str,
        nargs='+',
        required=True)

    parser.add_argument(
        '-b',
        '--blast_db',
        help='BLAST db',
        metavar='str',
        type=str,
        required=True)

    parser.add_argument(
        '-c',
        '--centroids_db',
        help='Centroids/SQLite db',
        metavar='str',
        type=str,
        default=None)

    parser.add_argument(
        '-p',
        '--blast_program',
        help='BLAST program',
        metavar='str',
        type=str,
        default='blastn')

    parser.add_argument(
        '-i',
        '--perc_identity',
        help='BLAST percent identity',
        metavar='float',
        type=float,
        default=0.)

    parser.add_argument(
        '-Q',
        '--qcov_hsp_perc',
        help='BLAST percent query coverage per hsp',
        metavar='float',
        type=float,
        default=0.)

    parser.add_argument(
        '-o',
        '--outdir',
        help='Output directory',
        metavar='str',
        type=str,
        default='out')

    return parser.parse_args()


# --------------------------------------------------
def warn(msg):
    """Print a message to STDERR"""
    print(msg, file=sys.stderr)


# --------------------------------------------------
def die(msg='Something bad happened'):
    """warn() and exit with error"""
    warn(msg)
    sys.exit(1)


# --------------------------------------------------
def main():
    """Make a jazz noise here"""
    args = get_args()
    query = args.query
    blast_db = args.blast_db
    centroids_db = args.centroids_db
    blast_prg = args.blast_program
    out_dir = args.outdir

    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)

    blast_dir = os.path.dirname(os.path.abspath(blast_db))
    if not os.path.isdir(blast_dir):
        die('--blast_db dir "{}" is not a dir'.format(blast_dir))

    blast_db_name = os.path.basename(blast_db)
    if not [f for f in os.listdir(blast_dir) if f.startswith(blast_db_name)]:
        die('No BLAST "{}" files in "{}"'.format(blast_db_name, blast_dir))

    if centroids_db and not os.path.isfile(centroids_db):
        die('--centroids "{}" is not a file'.format(centroids_db))

    print('Running BLAST')
    hits = run_blast(
        blast_db=blast_db,
        blast_prg=blast_prg,
        perc_identity=args.perc_identity,
        qcov_hsp_perc=args.qcov_hsp_perc,
        query=query,
        out_dir=out_dir)

    print('Centroids query')
    centroids_file = cmap_query(
        blast_hits=hits, centroids_db=centroids_db, out_dir=out_dir)

    print('Plotting')
    plot(centroids_file=centroids_file, out_dir=out_dir)

    print('Done')


# --------------------------------------------------
def plot(centroids_file, out_dir):
    """Given CMAP location data, plot distribution"""

    cwd = os.path.abspath(os.path.dirname(sys.argv[0]))
    plot = os.path.join(cwd, 'plot.r')

    if os.path.isfile(plot):
        rv, out = getstatusoutput('{} {}'.format(plot, centroids_file))
        print(out)

        if rv == 0:
            return 1
        else:
            die('Error plotting ({}):\n{}\n'.format(rv, out))
    else:
        die('Cannot find "{}"'.format(plot))


# --------------------------------------------------
def cmap_query(blast_hits, centroids_db, out_dir):
    """Given BLAST hits, query CMAP/SQLite for location"""

    # TODO: fix spelling of "esv_tempreature" => "esv_temperature"
    # if the CMAP tblESV column is changed
    flds = [
        'lat', 'lon', 'depth', 'relative_abundance', 'esv_tempreature',
        'esv_salinity', 'cruise_name'
    ]
    qry = 'select {} from tblesv where centroid=?'.format(', '.join(flds))
    cursor = sqlite3.connect(
        centroids_db).cursor() if centroids_db else db.dbConnect().cursor()

    out_file = os.path.join(out_dir, 'oce-input.csv')
    out_fh = open(out_file, 'wt')
    out_fh.write(','.join([
        'latitude', 'longitude', 'depth', 'Relative_Abundance', 'temperature',
        'salinity', 'cruise_name'
    ]) + '\n')

    seen = set()
    with open(blast_hits) as csvfile:
        blast_flds = [
            'qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
            'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'
        ]
        reader = csv.DictReader(csvfile, fieldnames=blast_flds, delimiter='\t')

        for hit_num, row in enumerate(reader):
            seq_id = row['sseqid']
            print('{:5}: {}'.format(hit_num + 1, seq_id))

            if seq_id in seen:
                continue
            seen.add(seq_id)

            rows = cursor.execute(qry, (seq_id, )).fetchall()

            if rows:
                for row in rows:
                    out_fh.write(','.join(map(lambda x: str(x).strip(), row)) +
                                 '\n')
            else:
                warn('Found no match for centroid "{}"'.format(seq_id))

    out_fh.close()
    return out_file


# --------------------------------------------------
def run_blast(blast_db, blast_prg, query, perc_identity, qcov_hsp_perc,
              out_dir):
    """Given user query and params, run BLAST"""

    hits_file = os.path.join(out_dir, 'hits.tab')

    cmd = '{} -query {} -db {} -out {} -outfmt 6'.format(
        blast_prg, query, blast_db, hits_file)

    if perc_identity > 0.:
        cmd += ' -perc_identity {}'.format(perc_identity)

    if qcov_hsp_perc > 0.:
        cmd += ' -qcov_hsp_perc {}'.format(qcov_hsp_perc)

    (rv, blast_out) = getstatusoutput(cmd)

    if rv == 0:
        with open(hits_file) as fh:
            if len(fh.readlines()) > 0:
                return hits_file
            else:
                die('No hits from BLAST')
    else:
        die('Failed to run BLAST ({}):\n{}'.format(rv, blast_out))


# --------------------------------------------------
if __name__ == '__main__':
    main()
