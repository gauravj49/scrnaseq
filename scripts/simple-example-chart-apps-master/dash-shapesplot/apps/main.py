#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import dash
import dash_core_components as dcc
import dash_html_components as html
import pandas as pd
import plotly.graph_objs as go
from app import app

df = pd.read_csv('https://raw.githubusercontent.com/plotly/datasets/master/gapminder2007.csv')

if 'DYNO' in os.environ:
    app_name = os.environ['DASH_APP_NAME']
else:
    app_name = 'dash-shapesplot'

layout = html.Div([
    html.Div([html.H1("Country Demographic Data by Continent")], style={'textAlign': "center", "padding-bottom": "30"}),
    html.Div(
        [dcc.Dropdown(id="continent-selected", options=[{'label': i, 'value': i} for i in df['continent'].unique()[:4]],
                      value="Asia",
                      style={"display": "block", "margin-left": "auto", "margin-right": "auto", "width": "50%"})]),
    dcc.Graph(id="my-graph")
], className="container")


@app.callback(
    dash.dependencies.Output("my-graph", "figure"),
    [dash.dependencies.Input("continent-selected", "value")])
def update_figure(selected):
    dff = df[df['continent'] == selected]
    trace = go.Scatter(x=dff["gdpPercap"], y=dff["lifeExp"], text=dff["country"], hoverinfo="text", mode="markers",
                       name=selected,
                       marker={'size': 10, 'line': {'width': 0.5, 'color': '#43a2ca'}, "color": "#7bccc4"}, )
    return {"data": [trace],
            "layout": go.Layout(title="Life Expectancy vs. GDP Per Capita", hovermode='closest',
                                xaxis={"title": 'GDP Per Capita (USD)', "range": [0, 40000], "tick0": 0, "dtick": 5000,
                                       "showline": True},
                                yaxis={"title": "Life Expectancy (Years)", "range": [20, 100], "tick0": 20, "dtick": 10,
                                       "showline": False},
                                shapes=[{'type': 'circle', 'xref': "x", 'yref': "y",
                                         'x0': dff["gdpPercap"].quantile(q=0.25), 'y0': dff["lifeExp"].quantile(q=0.25),
                                         'x1': dff["gdpPercap"].quantile(q=0.75), 'y1': dff["lifeExp"].quantile(q=0.75),
                                         'opacity': 0.6, 'fillcolor': '#cbebe8', 'line': {'color': '#69ded3'}}])}