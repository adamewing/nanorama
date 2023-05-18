#!/usr/bin/env python

import os
import sys
import argparse
import subprocess
import psutil

from time import sleep, time_ns
from collections import defaultdict as dd
from operator import itemgetter
from uuid import uuid4

import pysam
import pandas as pd

import plotly.express as px

from bx.intervals.intersection import Intersecter, Interval

import pandas as pd
import numpy as np

import logging
FORMAT = '%(asctime)s %(message)s'
logging.basicConfig(format=FORMAT)
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def map_reads(fastq, ref, read_data, targets, threads=8):

    FNULL = open(os.devnull, 'w')

    mm2_cmd  = ['minimap2', '-t', str(threads), '-x', 'map-ont', '-a', ref, fastq]
    view_cmd = ['samtools', 'view', '-u', '-']

    aln  = subprocess.Popen(mm2_cmd, stdout=subprocess.PIPE, stderr=FNULL)
    view = subprocess.Popen(view_cmd, stdin=aln.stdout, stdout=subprocess.PIPE, stderr=FNULL)

    save = pysam.set_verbosity(0)
    bamstream = pysam.AlignmentFile(view.stdout, 'rb')
    pysam.set_verbosity(save)

    r = 0

    for read in bamstream:
        if not read.is_unmapped and not read.is_secondary and not read.is_supplementary:
            target = 'NA'
            targeted = 'N'

            if read.reference_name in targets:
                for t in targets[read.reference_name].find(read.reference_start, read.reference_end):
                    target = f'{read.reference_name}:{t.start}-{t.end}'
                    targeted = 'Y'

            read_data.loc[str(uuid4())] = [target, targeted, read.reference_name, read.reference_start, read.reference_end, len(read.seq), np.log10(len(read.seq))]
            r += 1

    logger.info(f'mapped {r} reads from {fastq}')
    return read_data


def check_fopen(fn):
    for proc in psutil.process_iter():
        try:
            flist = proc.open_files()
            if flist:
                for nt in flist:
                    if os.path.basename(nt.path) == fn:
                        return True
        except:
            pass

    return False


def enrich(df, genome_len, target_len):
    enrich_data = pd.DataFrame(columns=['min_length', 'on_tgt_cov', 'off_tgt_cov', 'enrichment'])

    for minlen in range(100,1225,25):
        df_len = df[df['length'] >= minlen]
        offtgt_len = genome_len-target_len
        offtgt_cov = df_len[df_len['targeted'] == 'N']['length'].sum()
        target_cov = df_len[df_len['targeted'] == 'Y']['length'].sum()
        enr = (target_cov/target_len)/(offtgt_cov/offtgt_len)

        enrich_data.loc[str(uuid4())] = [minlen, target_cov, offtgt_cov, enr]

    return enrich_data


def main(args):
    chrom_bins = dd(dict)
    chrom_lens = {}

    read_data = pd.DataFrame(columns=['target', 'targeted', 'chrom', 'start', 'end', 'length', 'log10_length'])

    if os.path.exists(args.plotdir):
        logger.warning(f'output directory already exists: {args.plotdir}')
        sys.exit(1)

    else:
        os.mkdir(args.plotdir)

        assert os.path.exists(args.plotdir)

        with open(args.fai) as fai:
            for line in fai:
                chrom, length = line.strip().split()[:2]
                length = int(length)
                chrom_lens[chrom] = length

        genome_len = sum(list(chrom_lens.values()))
        logger.info(f'genome length: {genome_len}')

        targets = dd(Intersecter)
        target_len = 0

        with open(args.targets) as bed:
            for line in bed:
                c = line.strip().split()
                chrom, start, end = c[:3]
                start = int(start)
                end = int(end)
                target_len += (end-start)
                targets[chrom].add_interval(Interval(start, end))

        logger.info(f'total target length: {target_len}')

        with open(args.plotdir+'/stats.txt', 'w') as stats_out:
            stats_out.write(f'command line: {" ".join(sys.argv)}\n')
            stats_out.write(f'genome length: {genome_len}\ntarget length: {target_len}\n')

        done_fastqs = {}

        while True:
            for outdir in args.outdir[0]:
                for fn in os.listdir(outdir):
                    if fn.endswith('.fastq') or fn.endswith('.fastq.gz'):
                        if fn not in done_fastqs:
                            if check_fopen(fn):
                                logger.info(f'file is open: {fn}, skip for now')
                                continue

                            logger.info(f'mapping {fn} and updating read data')

                            read_data = map_reads(outdir+'/'+fn, args.mmi, read_data, targets)
                            read_data.to_csv(args.plotdir+'/read_data.csv')
                            logger.info(f'saved current data to {args.plotdir}/read_data.csv')

                            enrich_data = enrich(read_data, genome_len, target_len)
                            enrich_data.to_csv(args.plotdir+'/enrichment_data.csv')
                            logger.info(f'saved current data to {args.plotdir}/enrichment_data.csv')

                            done_fastqs[fn] = True

            sleep(2.0)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='fastq live coverage monitor')
    parser.add_argument('-o', '--outdir', required=True, help='minknow output directory containing fastqs, can be called more than once', action='append', nargs='+')
    parser.add_argument('-m', '--mmi', required=True, help='minimap2 index (.mmi)')
    parser.add_argument('-f', '--fai', required=True, help='reference genome .fai (from samtools faidx)')
    parser.add_argument('-t', '--targets', required=True, help='target .bed')
    parser.add_argument('-d', '--plotdir', required=True, help='output directory')
    args = parser.parse_args()
    main(args)

