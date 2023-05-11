#!/usr/bin/env python

import os
import argparse

from dash import Dash, html, dcc
import plotly.express as px
import pandas as pd

import logging
FORMAT = '%(asctime)s %(message)s'
logging.basicConfig(format=FORMAT)
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def load_data(plotdir):
    assert os.path.exists(plotdir)

    logger.info(f'load {args.plotdir}/read_data.csv')
    read_data = pd.read_csv(args.plotdir+'/read_data.csv')

    logger.info(f'load {args.plotdir}/enrichment_data.csv')
    enrich_data = pd.read_csv(args.plotdir+'/enrichment_data.csv')

    return read_data, enrich_data


def load_metrics(plotdir):
    assert os.path.exists(plotdir)

    genome_len = 0
    target_len = 0

    with open(args.plotdir+'/stats.txt') as stats:
        for line in stats:
            c = line.strip().split()
            if c[0] == 'genome':
                genome_len = int(c[-1])

            if c[0] == 'target':
                target_len = int(c[-1])

    return genome_len, target_len


def make_plots(read_data, enrich_data, plotdir, sample_size):
    sample_size = int(sample_size)

    if read_data.shape[0] < sample_size:
        sample_size = read_data.shape[0]

    sampled_data = read_data.sample(sample_size, random_state=1)

    figs = {}

    figs['strip'] = px.strip(sampled_data, x="targeted", y="log10_length", hover_name="target", width=500, height=500)
    figs['strip'].update_traces(marker=dict(size=5, opacity=0.2))

    figs['violin'] = px.violin(sampled_data, x="targeted", y="log10_length", width=500, height=500)

    figs['line'] = px.line(enrich_data, x='min_length', y='enrichment', markers=True, width=1000, height=400)
    figs['line'].update_layout(plot_bgcolor='#ffffff')

    return figs


def run_dash(args):
    app = Dash(__name__)

    genome_len, target_len = load_metrics(args.plotdir)
    read_data, enrich_data = load_data(args.plotdir)

    tgt_pct = '%.2f' % (target_len/genome_len*100)

    figs = make_plots(read_data, enrich_data, args.plotdir, args.samplesize)

    app.layout = html.Div(children=[
        html.H3(id='title', children='Live Read Monitor for ONT output'),

        html.Div(id='metrics', children=f'genome size: {genome_len} | target size: {target_len} ({tgt_pct}%)'),

        dcc.Graph(
            id='strip',
            figure=figs['strip'],
            style={'display': 'inline-block'}
        ),

        dcc.Graph(
            id='violin',
            figure=figs['violin'],
            style={'display': 'inline-block'}
        ),

        dcc.Graph(
            id='line',
            figure=figs['line']
        )
    ])

    app.run_server(host="127.0.0.1", port="8050", debug=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='fastq live coverage monitor')
    parser.add_argument('-d', '--plotdir', required=True, help='output directory')
    parser.add_argument('-s', '--samplesize', default=10000, help='sample size for plotting (default = 10000)')
    args = parser.parse_args()
    run_dash(args)

    

