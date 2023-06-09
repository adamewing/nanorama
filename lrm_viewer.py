#!/usr/bin/env python

import os
import sys
import argparse
from collections import defaultdict as dd

import dash
from dash import html, dcc
import pandas as pd
from plotly.subplots import make_subplots
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
    logger.info(f'loading targets from {bed_fn}')
    targets = dd(Intersecter)

    if bed_fn is None:
        return targets

    with open(bed_fn) as bed:
        for line in bed:
            c = line.strip().split()
            chrom, start, end = c[:3]
            start = int(start)
            end = int(end)

            value = None
            if len(c) > 3:
                value = c[3]

            targets[chrom].add_interval(Interval(start, end, value=value))
    
    return targets


def dash_lrt(plotdir, target_bed, tx_bed):
    if not os.path.exists(plotdir):
        sys.exit(f'directory not found: {plotdir}')
    
    targets = load_targets(target_bed)
    tx = load_targets(tx_bed)

    # Define the initial browser view
    initial_chr = 'chr19'
    initial_start = 44047754
    initial_end = 46011212

    genome_len, target_len = load_metrics(plotdir)
    tgt_pct = '%.2f' % (target_len/genome_len*100)

    global df
    df = load_data(plotdir)
    sz = os.path.getsize(plotdir+'/read_data.csv')

    # Filter the data for the initial view
    df_view = df[(df['chrom'] == initial_chr) & (df['start'] >= initial_start) & (df['end'] <= initial_end)]

    # Create initial scatter plot
    scatter_fig = px.scatter(df_view, x="start", y="log10_length", color="targeted", color_discrete_map={"Y": "red", "N": "blue"})

    app = dash.Dash(__name__)

    app.layout = html.Div([
        html.H3(id='title', children=f'Live Read Monitor for ONT output: {os.path.basename(plotdir)}'),
        html.Div(id='metrics', children=f'genome size: {genome_len} | target size: {target_len} ({tgt_pct}%)'),
        html.H4(id='genome-desc', children=f'Enrichment Browser:'),
        html.Div([
            dcc.Input(
                id='genome-coordinates-input',
                type='text',
                value=f'{initial_chr}:{initial_start}-{initial_end}'
            ),
            html.Button('Submit', id='submit-button', n_clicks=0),
        ]),
        dcc.Graph(id='scatter-plot', figure=scatter_fig, config={'scrollZoom': True}),
        dcc.Graph(id='violin-plot', style={'display': 'inline-block'}),
        dcc.Graph(id='strip-plot', style={'display': 'inline-block'}),
        dcc.Store(id="df-size", data=sz),
        html.H4(id='curve-desc', children=f'Overall enrichment:'),
        
        dcc.Graph(id='curve-plot'),
        html.Div([
            html.Button('Update Enrichment', id='update-curve-button')
        ])
    ])


    @app.callback(
        [Output('scatter-plot', 'figure'),
        Output('violin-plot', 'figure'),
        Output('strip-plot', 'figure')],

        [Input('submit-button', 'n_clicks'),
        Input('scatter-plot', 'relayoutData'),
        Input('df-size', 'data')],

        [State('genome-coordinates-input', 'value')]
    )
    def update_figure(n_clicks, relayoutData, sz, value):
        current_sz = os.path.getsize(plotdir+'/read_data.csv')
        updated_sz = current_sz != sz

        global df

        if updated_sz:
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

        #scatter_fig = go.Figure()

        samples = df['sample'].unique()

        scatter_fig = make_subplots(rows=len(samples), cols=1, shared_xaxes=True, vertical_spacing=0)

        for i, sample in enumerate(samples):
            df_subset = df_view[df_view['sample'] == sample]

            scatter_fig.add_trace(go.Scatter(
                x=df_subset['start'],
                y=df_subset['log10_length'],
                mode='markers',
                marker=dict(
                    color=df_subset['targeted'].apply(lambda x: 'red' if x == 'Y' else 'blue'),
                    opacity=0.5
                ),
            ),
            row = i+1,
            col = 1,
            )

        scatter_fig.update_layout(
            #xaxis_title="Genome Position (green = targets)",
            #yaxis_title="Log Read Length",
            height=600,
            paper_bgcolor='rgba(0,0,0,0)',
            plot_bgcolor='rgba(0,0,0,0)',
            dragmode="pan",
            autosize=True,
            showlegend=False,
        )

        scatter_fig.update_xaxes(range=[df_view['start'].min(), df_view['end'].max()])
        scatter_fig.update_yaxes(range=[df_view['log10_length'].min()-1.2, df_view['log10_length'].max()+.2], fixedrange=True)

        # draw targets
        
        tgt_track_ymax = max(df_view['log10_length']) + 0.1
        tgt_track_ymin = min(df_view['log10_length']) - 0.1
        
        target_boxes= []
        for i in range(len(samples)):
            xref = f'x{i+1}'
            yref = f'y{i+1}'
    
            for rec in targets[chr].find(start, end):
                target_boxes.append({
                    'type': 'rect',
                    'xref': xref,
                    'yref': yref,
                    'x0': rec.start,
                    'y0': tgt_track_ymax,
                    'x1': rec.end,
                    'y1': tgt_track_ymin,
                    'fillcolor': 'green',
                    'opacity': 0.2,
                    'line': {
                        'width': 0,
                    },
                })
        

        # transcript track
    
        tx_track_ymax = tgt_track_ymin - 0.1
        tx_track_ymin = tx_track_ymax - 0.8
        tx_track_height = tx_track_ymax-tx_track_ymin

        stacklevel = {} # avoid gene overlaps

        for rec in tx[chr].find(start, end):
            overlaps = []
            for o_rec in tx[chr].find(rec.start, rec.end):
                if o_rec.value not in overlaps:
                    overlaps.append(o_rec.value)

            for i, gene_name in enumerate(overlaps):
                stacklevel[gene_name] = i

        max_stacklevel = max(stacklevel.values())

        tx_item_height = tx_track_height/(max_stacklevel+1)

        font_vis = 'rgba(0,0,0,0)'

        if end-start < 1e6:
            font_vis = 'black'

        for i in range(len(samples)):
            xref = f'x{i+1}'
            yref = f'y{i+1}'

            for rec in tx[chr].find(start, end):
                sl = stacklevel[rec.value]
                y0 = tx_track_ymax - tx_item_height*sl
                y1 = y0 - tx_item_height
                target_boxes.append({
                    'type': 'rect',
                    'xref': xref,
                    'yref': yref,
                    'x0': rec.start,
                    'y0': y0,
                    'x1': rec.end,
                    'y1': y1,
                    'fillcolor': 'orange',
                    'opacity': 0.5,
                    'line': {
                        'width': 1,
                    },
                    'label': {
                        'text': rec.value,
                        'font': {'color': font_vis}
                    },
                })

        scatter_fig.update_layout(shapes=target_boxes)

        violin_fig = px.violin(df_view,
                               x="targeted",
                               y="log10_length", 
                               width=500,
                               height=500,
                               color="targeted",
                               facet_col="sample",
                               color_discrete_map={"Y": "red", "N": "blue"}
                               )
        
        strip_fig = px.strip(df_view,
                             x="targeted",
                             y="log10_length",
                             hover_name="target",
                             width=500,
                             height=500,
                             color="targeted",
                             facet_col="sample",
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

        fig = px.line(df_enrich, x='min_length', y='enrichment', line_group='sample', color='sample', markers=True, width=1000, height=400)
        fig.update_layout(plot_bgcolor='#ffffff')

        max_e = int(df_enrich['enrichment'].max()+1)
        fig.update_layout(yaxis_range=[0,max_e])

        return fig

    if args.external:
        app.run_server(host='0.0.0.0', debug=True)

    else:
        app.run_server(debug=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='live coverage monitor aimed at fastqs generated by ONT devices')
    parser.add_argument('-d', '--plotdir', required=True, help='output directory')
    parser.add_argument('-t', '--targets', required=True, help='target .bed')
    parser.add_argument('--tx', default=None, help='transcripts (or any other annotation) .bed')
    parser.add_argument('--external', action='store_true', default=False)
    args = parser.parse_args()
    dash_lrt(args.plotdir, args.targets, args.tx)
    

