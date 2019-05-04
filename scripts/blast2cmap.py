#!/usr/bin/env python3
"""
Author : Ken Youens-Clark <kyclark@email.arizona.edu>
Date   : 2019-02-01
Purpose: BLAST hits to CMAP visualization
"""

import argparse
import csv
import os
import shutil
import sqlite3
import subprocess
import sys
import tempfile
import pandas as pd
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
        choices=['blastn', 'tblastn'],
        default='blastn')

    parser.add_argument(
        '-i',
        '--perc_identity',
        help='BLAST percent identity',
        metavar='float',
        type=float,
        default=97.)

    parser.add_argument(
        '-Q',
        '--qcov_hsp_perc',
        help='BLAST percent query coverage per hsp',
        metavar='float',
        type=float,
        default=100.)

    parser.add_argument(
        '-o',
        '--outdir',
        help='Output directory',
        metavar='str',
        type=str,
        default='out')

    parser.add_argument(
        '-n',
        '--num_cpus',
        help='Num CPUs for parallel',
        metavar='int',
        type=int,
        default=8)

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
    frac_files = cmap_query(
        blast_hits=hits,
        centroids_db=centroids_db,
        out_dir=out_dir,
        perc_identity=args.perc_identity)

    plot(frac_files=frac_files, out_dir=out_dir, num_procs=args.num_cpus)

    print('Done')


# --------------------------------------------------
def plot(frac_files, out_dir, num_procs):
    """Given CMAP location data, plot distribution"""

    cwd = os.path.abspath(os.path.dirname(sys.argv[0]))
    plot_r = os.path.join(cwd, 'plot.r')

    if not os.path.isfile(plot_r):
        die('Missing "{}"'.format(plot_r))

    jobfile = tempfile.NamedTemporaryFile(delete=False, mode='wt')
    job_tmpl = "{prg} -f {file} -o {out_dir} -t '{title}' -l '{legend}'\n"
    for file in frac_files:
        base, _ = os.path.splitext(os.path.basename(file))
        frac_dir = os.path.join(out_dir, base)
        d = dict(map(lambda s: s.split('_', maxsplit=1), base.split('__')))
        title = 'eASV* Abundance For Cruise "{}" ({} μM size frac.)'.format(
            d['cruise'], d['frac'])
        legend = '*eASV ID="{}" is {}% similar to query "{}"'.format(
            d['asv'], d['pident'], d['qseqid'])
        jobfile.write(
            job_tmpl.format(
                prg=plot_r,
                file=file,
                out_dir=frac_dir,
                title=title,
                legend=legend))

    jobfile.close()

    run_job_file(jobfile=jobfile.name, msg='Plotting', num_procs=num_procs)

    return 1


# --------------------------------------------------
def cmap_query(blast_hits, centroids_db, out_dir, perc_identity):
    """Given BLAST hits, query CMAP/SQLite for location"""

    # TODO: fix spelling of "esv_tempreature" => "esv_temperature"
    # if the CMAP tblESV column is changed
    qry_flds = [
        'centroid', 'lat', 'lon', 'depth', 'relative_abundance',
        'esv_tempreature', 'esv_salinity', 'cruise_name', 'size_frac_lower',
        'size_frac_upper'
    ]
    qry = 'select {} from tblesv where centroid=?'.format(', '.join(qry_flds))
    cursor = sqlite3.connect(
        centroids_db).cursor() if centroids_db else db.dbConnect().cursor()

    out_flds = [
        'centroid', 'latitude', 'longitude', 'depth', 'relative_abundance',
        'temperature', 'salinity', 'cruise_name', 'size_frac_lower',
        'size_frac_upper', 'pident', 'qseqid'
    ]

    data_dir = os.path.join(out_dir, 'data')

    if not os.path.isdir(data_dir):
        os.makedirs(data_dir)

    out_file = os.path.join(data_dir, 'oce-input.csv')
    out_fh = open(out_file, 'wt')
    out_fh.write(','.join(out_flds) + '\n')

    seen = set()
    with open(blast_hits) as csvfile:
        blast_flds = [
            'qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
            'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'
        ]
        reader = csv.DictReader(csvfile, fieldnames=blast_flds, delimiter='\t')
        write_flds = qry_flds + ['pident', 'qseqid']

        for hit_num, blast_hit in enumerate(reader):
            seq_id = blast_hit['sseqid']
            print('{:5}: {}'.format(hit_num + 1, seq_id))

            if seq_id in seen:
                continue
            seen.add(seq_id)

            rows = cursor.execute(qry, (seq_id, )).fetchall()

            if rows:
                for row in rows:
                    d = dict(zip(qry_flds, row))
                    d['pident'] = blast_hit['pident']
                    d['qseqid'] = blast_hit['qseqid']

                    upper = d['size_frac_upper']
                    if not isinstance(upper, (int, float)):
                        d['size_frac_upper'] = d['size_frac_lower']

                    out_fh.write(','.join(
                        map(lambda x: str(d[x]).strip(), write_flds)) + '\n')
            else:
                warn('Found no match for centroid "{}"'.format(seq_id))

    out_fh.close()

    frac_files = []
    df = pd.read_csv(out_file)
    df['size'] = df['size_frac_lower'].astype(
        str) + '-' + df['size_frac_upper'].astype(str)

    for centroid in df['centroid'].unique():
        for cruise_name in df['cruise_name'].unique():
            for frac in df['size'].unique():
                frac_df = df[(df['centroid'] == centroid)
                             & (df['cruise_name'] == cruise_name) &
                             (df['size'] == frac)]

                if not frac_df.empty:
                    pident = frac_df['pident'].unique()[0]
                    qseqid = frac_df['qseqid'].unique()[0]
                    t = '__'.join([
                        'asv_{}', 'cruise_{}', 'qseqid_{}', 'pident_{:.02f}',
                        'frac_{}'
                    ]) + '.csv'
                    frac_out = os.path.join(
                        data_dir,
                        t.format(centroid, cruise_name, qseqid, pident, frac))
                    frac_df.to_csv(frac_out, index=False)
                    frac_files.append(frac_out)

    return frac_files


# --------------------------------------------------
def line_count(fname):
    """Count the number of lines in a file"""

    n = 0
    for _ in open(fname):
        n += 1

    return n


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
def run_job_file(jobfile, msg='Running job', num_procs=1):
    """Run a job file if there are jobs"""

    num_jobs = line_count(jobfile)
    warn('{} (# jobs = {} # cpus = {})'.format(msg, num_jobs, num_procs))

    if num_jobs > 0:
        try:
            if shutil.which('parallel'):
                try:
                    cmd = 'parallel --halt soon,fail=1 -P {} < {}'.format(
                        num_procs, jobfile)
                    subprocess.run(cmd, shell=True, check=True)
                except subprocess.CalledProcessError as err:
                    die('Error:\n{}\n{}\n'.format(err.stderr, err.stdout))
            else:
                for cmd in open(jobfile):
                    subprocess.run(cmd, shell=True, check=True)
        except Exception as e:
            die('Error: {}'.format(e))
        finally:
            os.remove(jobfile)

    return True


# --------------------------------------------------
if __name__ == '__main__':
    main()
