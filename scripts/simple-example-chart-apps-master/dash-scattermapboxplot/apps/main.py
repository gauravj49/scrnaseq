#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import dash
import us
import dash_core_components as dcc
import dash_html_components as html
import pandas as pd
import plotly.graph_objs as go
from app import app

mapbox_access_token = "pk.eyJ1IjoicHJpeWF0aGFyc2FuIiwiYSI6ImNqbGRyMGQ5YTBhcmkzcXF6YWZldnVvZXoifQ.sN7gyyHTIq1BSfHQRBZdHA"

df1 = pd.read_csv("https://raw.githubusercontent.com/plotly/datasets/master/2011_february_us_airport_traffic.csv")
df = df1.dropna(axis=0)

if 'DYNO' in os.environ:
    app_name = os.environ['DASH_APP_NAME']
else:
    app_name = 'dash-scattermapboxplot'

layout = html.Div([
    html.Div([html.H1("Airport Locations in the United States")],
             style={'textAlign': "center", "padding-bottom": "10", "padding-top": "10"}),
    html.Div([dcc.Dropdown(id="state-selected", value=['CA'], multi=True,
                           options=[{'label': f'{us.states.lookup(i)}', 'value': i} for i in df.state.unique()],
                           style={"display": "block", "margin-left": "auto", "margin-right": "auto", "width": "50%"})]),
    html.Div(dcc.Graph(id="my-graph"))
], className="container")


@app.callback(
    dash.dependencies.Output("my-graph", "figure"),
    [dash.dependencies.Input("state-selected", "value")])
def update_figure(selected):
    trace = []
    for state in selected:
        dff = df[df["state"] == state]
        trace.append(
            go.Scattermapbox(lat=dff["lat"], lon=dff["long"], mode='markers', marker={'symbol': "airport", 'size': 10},
                             text=dff['airport'], hoverinfo='text', name=state))
    return {"data": trace,
            "layout": go.Layout(autosize=True, hovermode='closest', showlegend=False, height=700,
                                mapbox={'accesstoken': mapbox_access_token, 'bearing': 0,
                                        'center': {'lat': 38, 'lon': -94}, 'pitch': 30, 'zoom': 3,
                                        "style": 'mapbox://styles/mapbox/light-v9'})}
