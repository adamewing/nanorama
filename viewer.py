#!/usr/bin/env python

import os
import sys
import dash
from dash import html, dcc
import pandas as pd
import plotly.graph_objs as go
import plotly.express as px
from dash.dependencies import Input, Output

import logging
FORMAT = '%(asctime)s %(message)s'
logging.basicConfig(format=FORMAT)
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def load_data(plotdir):
    logger.info(f'load {plotdir}/read_data.csv')
    return pd.read_csv(plotdir+'/read_data.csv')
    

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


plotdir = sys.argv[1]

if not os.path.exists(plotdir):
    sys.exit(f'directory not found: {plotdir}')

df = load_data(plotdir)
df_enrich = load_enrich(plotdir)


# Define the initial browser view
initial_chr = 'chr11'
initial_start = 1
initial_end = 20000000

genome_len, target_len = load_metrics(plotdir)
tgt_pct = '%.2f' % (target_len/genome_len*100)

# Filter the data for the initial view
df_view = df[(df['chrom'] == initial_chr) & (df['start'] >= initial_start) & (df['end'] <= initial_end)]

app = dash.Dash(__name__)

app.layout = html.Div([
    html.H3(id='title', children='Live Read Monitor for ONT output'),
    html.Div(id='metrics', children=f'genome size: {genome_len} | target size: {target_len} ({tgt_pct}%)'),

    html.H4(id='curve-desc', children=f'Overall enrichment:'),
    
    dcc.Graph(id='curve-plot'),
    html.Div([
        html.Button('Update Enrichment', id='update-button')
    ]),
    
    html.H4(id='genome-desc', children=f'Enrichment Browser:'),
    html.Div([
        dcc.Input(
            id='genome-coordinates-input',
            type='text',
            value=f'{initial_chr}:{initial_start}-{initial_end}'
        ),
        html.Button('Submit', id='submit-button')
    ]),
    dcc.Graph(id='scatter-plot'),
    dcc.Graph(id='violin-plot', style={'display': 'inline-block'}),
    dcc.Graph(id='strip-plot', style={'display': 'inline-block'})

])


@app.callback(
    Output('scatter-plot', 'figure'),
    Input('submit-button', 'n_clicks'),
    Input('genome-coordinates-input', 'value')
)
def update_figure(n_clicks, value):
    chr, coords = value.split(':')
    start, end = map(int, coords.split('-'))

    df_view = df[(df['chrom'] == chr) & (df['start'] >= start) & (df['end'] <= end)]

    fig = go.Figure()
    fig.add_trace(go.Scatter(
        x=df_view['start'],
        y=df_view['log10_length'],
        mode='markers',
        marker=dict(
            color=df_view['targeted'].apply(lambda x: 'red' if x == 'Y' else 'blue')
        )
    ))

    fig.update_layout(
        xaxis_title="Genome Position",
        yaxis_title="Log Read Length",
        dragmode="pan",
        autosize=True,
    )

    return fig


@app.callback(
    Output('violin-plot', 'figure'),
    Input('submit-button', 'n_clicks'),
    Input('genome-coordinates-input', 'value')
)
def update_violin_figure(n_clicks, value):
    chr, coords = value.split(':')
    start, end = map(int, coords.split('-'))

    df_view = df[(df['chrom'] == chr) & (df['start'] >= start) & (df['end'] <= end)]

    fig = px.violin(df_view, x="targeted", y="log10_length", width=500, height=500)
    fig.update_layout(plot_bgcolor='#ffffff')

    return fig


@app.callback(
    Output('strip-plot', 'figure'),
    Input('submit-button', 'n_clicks'),
    Input('genome-coordinates-input', 'value')
)
def update_strip_figure(n_clicks, value):
    chr, coords = value.split(':')
    start, end = map(int, coords.split('-'))

    df_view = df[(df['chrom'] == chr) & (df['start'] >= start) & (df['end'] <= end)]

    fig = px.strip(df_view, x="targeted", y="log10_length", hover_name="target", width=500, height=500)
    fig.update_layout(plot_bgcolor='#ffffff')

    return fig


@app.callback(
    Output('curve-plot', 'figure'),
    Input('update-button', 'n_clicks')
)
def update_curve_figure(n_clicks):

    df_enrich = load_enrich(plotdir)

    fig = px.line(df_enrich, x='min_length', y='enrichment', markers=True, width=1000, height=400)
    fig.update_layout(plot_bgcolor='#ffffff')

    return fig


if __name__ == '__main__':
    app.run_server(debug=True)

