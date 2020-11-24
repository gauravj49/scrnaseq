#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import dash_core_components as dcc
import dash_html_components as html
import pandas as pd
import plotly.graph_objs as go
from dash.dependencies import Input, Output
from app import app

df = pd.read_csv("https://raw.githubusercontent.com/plotly/datasets/master/nz_weather.csv")
available_indicators = df.columns.values[1:]

if 'DYNO' in os.environ:
    app_name = os.environ['DASH_APP_NAME']
else:
    app_name = 'dash-multipleaxesplot'

layout = html.Div([html.Div([html.H1("Rainfall for Cities in New Zealand ")], style={'textAlign': "center"}),
                   html.Div([html.Div(dcc.Dropdown(id='yaxis-column1',
                                                   options=[{'label': i, 'value': i} for i in available_indicators],
                                                   value='Auckland'), className="six columns"),
                             html.Div(dcc.Dropdown(id='yaxis-column2',
                                                   options=[{'label': i, 'value': i} for i in available_indicators],
                                                   value='Wellington'), className="six columns")], className="row"),
                   dcc.Graph(id="my-graph")
                   ], className="container", style={"padding-right": 10, 'padding-left': 10})


@app.callback(
    Output('my-graph', 'figure'),
    [Input('yaxis-column1', 'value'), Input('yaxis-column2', 'value')])
def update_figure(selected_1, selected_2):
    trace1 = go.Scatter(x=df["DATE"], y=df[selected_1], name=selected_1, mode='lines', opacity=0.6)
    trace2 = go.Scatter(x=df["DATE"], y=df[selected_2], name=selected_2, yaxis="y2", xaxis='x', mode='lines',
                        marker={"size": 8}, opacity=0.7)
    traces = [trace1, trace2]
    layout = go.Layout(title=f"Rainfall for {selected_1} and {selected_2}", margin={"l": 100, "r": 100},
                       colorway=["#287D95", "#EF533B"], legend={"x": 0.7, "y": 1, 'orientation': "h"},
                       yaxis={'title': f'Rainfall (mm)  for {selected_1}', "range": [0, 300]},
                       yaxis2={'title': f'Rainfall (mm)  for {selected_2}', 'overlaying': 'y', 'side': 'right',
                               "range": [0, 300], "showgrid": False},
                       xaxis={"title": "Date"})
    fig = go.Figure(data=traces, layout=layout)
    return fig
