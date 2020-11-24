#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import dash
import dash_core_components as dcc
import dash_html_components as html
import pandas as pd
import plotly.graph_objs as go
from app import app
from dash.dependencies import Input, Output

df = pd.read_csv('https://raw.githubusercontent.com/plotly/datasets/master/nhl_draft_2013_%40thejustinfisher.csv')

if 'DYNO' in os.environ:
    app_name = os.environ['DASH_APP_NAME']
else:
    app_name = 'dash-tableplot'

layout = html.Div([
    html.Div([html.H1("NHL Entry Draft 2013")], className="row", style={'textAlign': "center"}),
    html.Div([html.Span("Criteria to compare", style={'textAlign': "center", "display": "block"}),
              dcc.Dropdown(id='column-selected',multi=True, value=['Height'],
                           options=[{"label": i, 'value': i} for i in df.columns[1:5].append(df.columns[8:14])],
                           style={"display": "block", "margin-left": "auto", "margin-right": "auto", "width": "60%"})
              ], className="row", style={"padding": 10}),
    html.Div([dcc.Graph(id="my-graph")], className="row", style={"padding": 20}),
], className="container")


@app.callback(
    Output("my-graph", "figure"),
    [Input("column-selected", "value")])
def update_graph(column):
    value_header = ['Name', 'Team', 'Position']
    value_cell = [df['Name'], df['Team'], df['Position']]
    for col in column:
        value_header.append(col)
        value_cell.append(df[col])
    trace = go.Table(
        header={"values": value_header, "fill": {"color": "#FFD957"}, "align": ['center'], "height": 35,
                "line": {"width": 2, "color": "#685000"}, "font": {"size": 15}},
        cells={"values": value_cell, "fill": {"color": "#FFE89A"}, "align": ['left', 'center'],
               "line": {"color": "#685000"}})
    layout = go.Layout(title=f"Entry Draft", height=600)
    return {"data": [trace], "layout": layout}