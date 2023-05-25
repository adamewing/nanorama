#!/usr/bin/env python

import os
import sys
import argparse
from collections import defaultdict as dd

import dash
from dash import html, dcc
import pandas as pd
import plotly.graph_objs as go
import plotly.express as px
from dash.dependencies import Input, Output, State

from bx.intervals.intersection import Intersecter, Interval

import logging
FORMAT = '%(asctime)s %(message)s'
logging.basicConfig(format=FORMAT)
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def load_data(plotdir):
    logger.info(f'load {plotdir}/read_data.csv')
    #sz = os.path.getsize(plotdir+'/read_data.csv')
    return pd.read_csv(plotdir+'/read_data.csv')


def new_data(plotdir, old_sz):
    if old_sz != os.path.getsize(plotdir+'/read_data.csv'):
        return True
    return False


def load_enrich(plotdir):
    logger.info(f'load {plotdir}/enrichment_data.csv')
    return pd.read_csv(plotdir+'/enrichment_data.csv')


def load_metrics(plotdir):
    assert os.path.exists(plotdir)

    genome_len = 0
    target_len = 0

    with open(plotdir+'/stats.txt') as stats:
        for line in stats:
            c = line.strip().split()
            if c[0] == 'genome':
                genome_len = int(c[-1])

            if c[0] == 'target':
                target_len = int(c[-1])

    return genome_len, target_len


def load_targets(bed_fn):
    targets = dd(Intersecter)

    with open(bed_fn) as bed:
        for line in bed:
            c = line.strip().split()
            chrom, start, end = c[:3]
            start = int(start)
            end = int(end)

            targets[chrom].add_interval(Interval(start, end))
    
    return targets


def dash_lrt(plotdir, target_bed):
    if not os.path.exists(plotdir):
        sys.exit(f'directory not found: {plotdir}')
    
    targets = load_targets(target_bed)

    # Define the initial browser view
    initial_chr = 'chr11'
    initial_start = 1
    initial_end = 20000000

    genome_len, target_len = load_metrics(plotdir)
    tgt_pct = '%.2f' % (target_len/genome_len*100)

    df = load_data(plotdir)

    # Filter the data for the initial view
    df_view = df[(df['chrom'] == initial_chr) & (df['start'] >= initial_start) & (df['end'] <= initial_end)]

    app = dash.Dash(__name__)

    app.layout = html.Div([
        html.H3(id='title', children='Live Read Monitor for ONT output'),
        html.Div(id='metrics', children=f'genome size: {genome_len} | target size: {target_len} ({tgt_pct}%)'),

        html.H4(id='curve-desc', children=f'Overall enrichment:'),
        
        dcc.Graph(id='curve-plot'),
        html.Div([
            html.Button('Update Enrichment', id='update-curve-button')
        ]),
        
        html.H4(id='genome-desc', children=f'Enrichment Browser:'),
        html.Div([
            dcc.Input(
                id='genome-coordinates-input',
                type='text',
                value=f'{initial_chr}:{initial_start}-{initial_end}'
            ),
            html.Button('Submit', id='submit-button', n_clicks=0)
        ]),
        dcc.Graph(id='scatter-plot'),
        dcc.Graph(id='violin-plot', style={'display': 'inline-block'}),
        dcc.Graph(id='strip-plot', style={'display': 'inline-block'})
    ])


    @app.callback(
        [Output('scatter-plot', 'figure'),
        Output('violin-plot', 'figure'),
        Output('strip-plot', 'figure')],

        [Input('submit-button', 'n_clicks'),
        Input('scatter-plot', 'relayoutData')],

        [State('genome-coordinates-input', 'value')]
    )
    def update_figure(n_clicks, relayoutData, value):

        df = load_data(plotdir)

        ctx = dash.callback_context
        if not ctx.triggered:
            chr, coords = initial_chr, f"{initial_start}-{initial_end}"
        elif ctx.triggered[0]['prop_id'].split('.')[0] == 'submit-button':
            chr, coords = value.split(':')
        else:
            chr, coords = value.split(':')
            if 'xaxis.range[0]' in relayoutData:
                start, end = int(relayoutData['xaxis.range[0]']), int(relayoutData['xaxis.range[1]'])
                start = max(0, start)  # Ensure start is not less than 0
                coords = f"{start}-{end}"
                value = f"{chr}:{coords}"

        start, end = map(int, coords.split('-'))

        df_view = df[(df['chrom'] == chr) & (df['start'] >= start) & (df['end'] <= end)]

        scatter_fig = go.Figure()
        scatter_fig.add_trace(go.Scatter(
            x=df_view['start'],
            y=df_view['log10_length'],
            mode='markers',
            marker=dict(
                color=df_view['targeted'].apply(lambda x: 'red' if x == 'Y' else 'blue')
            )
        ))

        scatter_fig.update_layout(
            xaxis_title="Genome Position (green = targets)",
            yaxis_title="Log Read Length",
            dragmode="pan",
            autosize=True,
        )

        scatter_fig.update_xaxes(range=[df_view['start'].min(), df_view['end'].max()])

        # draw targets
        
        track_ymax = min(df_view['log10_length']) - 0.25
        track_ymin = track_ymax- 0.5

        target_boxes= []
        for rec in targets[chr].find(start, end):
            target_boxes.append({
                'type': 'rect',
                'xref': 'x',
                'yref': 'y',
                'x0': rec.start,
                'y0': track_ymax,  # adjust y0 and y1 according to the desired track position
                'x1': rec.end,
                'y1': track_ymin,
                'fillcolor': 'green',  # adjust color and opacity as desired
                'opacity': 0.5,
                'line': {
                    'width': 0,
                },
            })
        scatter_fig.update_layout(shapes=target_boxes)

        violin_fig = px.violin(df_view,
                               x="targeted",
                               y="log10_length", 
                               width=500,
                               height=500,
                               color="targeted",
                               color_discrete_map={"Y": "red", "N": "blue"}
                               )
        
        strip_fig = px.strip(df_view,
                             x="targeted",
                             y="log10_length",
                             hover_name="target",
                             width=500,
                             height=500,
                             color="targeted",
                             color_discrete_map={"Y": "red", "N": "blue"}
                             )

        strip_fig.update_traces(marker=dict(opacity=0.5))

        return scatter_fig, violin_fig, strip_fig


    @app.callback(
        Output('genome-coordinates-input', 'value'),
        Input('scatter-plot', 'relayoutData'),
        State('genome-coordinates-input', 'value')
    )
    def update_input(relayoutData, value):
        if relayoutData and 'xaxis.range[0]' in relayoutData:
            chr, _ = value.split(':')
            start, end = int(relayoutData['xaxis.range[0]']), int(relayoutData['xaxis.range[1]'])
            start = max(0, start) # Ensure start is not less than 0
            return f"{chr}:{start}-{end}"
        return value


    @app.callback(
        Output('curve-plot', 'figure'),
        Input('update-curve-button', 'n_clicks')
    )
    def update_curve_figure(n_clicks):

        df_enrich = load_enrich(plotdir)

        fig = px.line(df_enrich, x='min_length', y='enrichment', markers=True, width=1000, height=400)
        fig.update_layout(plot_bgcolor='#ffffff')

        return fig

    app.run_server(debug=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='live coverage monitor aimed at fastqs generated by ONT devices')
    parser.add_argument('-d', '--plotdir', required=True, help='output directory')
    parser.add_argument('-t', '--targets', required=True, help='target .bed')
    args = parser.parse_args()
    dash_lrt(args.plotdir, args.targets)
    

