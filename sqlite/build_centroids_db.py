#!/usr/bin/env python3
"""
Author : Ken Youens-Clark <kyclark@email.arizona.edu>
Date   : 2019-03-01
Purpose: Build Centroids db
"""

import argparse
import gzip
import os
import sys
import sqlite3
from opedia import db
from subprocess import getstatusoutput
from Bio import SeqIO


# --------------------------------------------------
def get_args():
    """get command-line arguments"""
    parser = argparse.ArgumentParser(
        description='Build centroids db',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument(
        'fasta', metavar='FILE', help='FASTA file of centroid sequences')

    parser.add_argument(
        '-d',
        '--dbname',
        help='SQLite DB name',
        metavar='DIR',
        type=str,
        default='centroids.db')

    parser.add_argument(
        '-o',
        '--outdir',
        help='Output directory',
        metavar='DIR',
        type=str,
        default='')

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
def import_centroid(sqlite_db, cmap_db, centroid_hash):
    """Import centroid"""

    select = "select * from tblesv where centroid='{}'"
    data = cmap_db.dbFetch(select.format(centroid_hash))
    if data.shape[0] == 0:
        print('Found nothing for "{}"'.format(centroid_hash))
        return 0

    cur = sqlite_db.cursor()

    # Remove all centroid data first
    cur.execute('delete from tblesv where centroid=?',
                (centroid_hash, ))

    flds = [
        'lat', 'lon', 'depth', 'relative_abundance', 'esv_temperature',
        'esv_salinity', 'cruise_name', 'size_frac_lower', 'size_frac_upper'
    ]

    insert = 'insert into tblesv (centroid, {}) values (?, {})'.format(
        ', '.join(flds), ', '.join(['?' for f in flds]))

    num_imported = 0
    for j, row in data.iterrows():
        vals = [centroid_hash] + list(
            map(lambda f: str(row[f]).rstrip(), flds))
        cur.execute(insert, vals)
        num_imported += 1

    sqlite_db.commit()

    return num_imported


# --------------------------------------------------
def sqlite_connect(db_name, cwd):
    if not os.path.isfile(db_name):
        print('Cannot find db "{}", creating...'.format(db_name))
        create = os.path.join(cwd, 'create.sql')
        if not os.path.isfile(create):
            die('Cannot find "{}"'.format(create))

        rv, out = getstatusoutput('sqlite3 {} < {}'.format(db_name, create))
        if rv != 0: die('Failed to create "{}"\n{}'.format(db_name, out))

    return sqlite3.connect(db_name)


# --------------------------------------------------
def main():
    """Make a jazz noise here"""
    args = get_args()
    fasta = args.fasta
    db_name = args.dbname

    if not os.path.isfile(fasta):
        die('"{}" is not a file'.format(fasta))

    cwd = os.path.abspath(os.path.dirname(sys.argv[0]))
    out_dir = args.outdir or cwd

    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)

    sqlite_db = sqlite_connect(db_name, cwd)
    total_imported = 0

    _, ext = os.path.splitext(os.path.basename(fasta))
    fh = gzip.open(fasta, 'rt') if ext == '.gz' else open(fasta)

    for i, rec in enumerate(SeqIO.parse(fh, 'fasta'), start=1):
        centroid_hash = rec.id
        num_imported = import_centroid(sqlite_db, db, centroid_hash)
        print('{:6}: {} ({})'.format(i, centroid_hash, num_imported))
        total_imported += num_imported

    print('Done, imported {} centroid values.'.format(total_imported))


# --------------------------------------------------
if __name__ == '__main__':
    main()
