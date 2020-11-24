#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
from app import app
import dash_core_components as dcc
import dash_html_components as html
import pandas as pd
import plotly.graph_objs as go
from dash.dependencies import Input, Output

df = pd.read_csv('https://raw.githubusercontent.com/plotly/datasets/master/Emissions%20Data.csv')

if 'DYNO' in os.environ:
    app_name = os.environ['DASH_APP_NAME']
else:
    app_name = 'dash-boxplot'

layout = html.Div([html.Div([
        html.H1("Greenhouse Gas Emissions by Continent")], style={"textAlign": "center"}),
        dcc.Graph(id="my-graph"),
        html.Div([dcc.RangeSlider(id="selected-continent",min=2008,max=2011,
                 marks={2008: "2008", 2009: "2009", 2010: "2010", 2011: "2011"},value=[2008, 2010],)],
                 style={"display": "block","margin-left": "auto","margin-right": "auto","width": "60%"}),
], className="container")


@app.callback(
    Output('my-graph', 'figure'),
    [Input('selected-continent', 'value')])
def update_figure(selected):
    dff = df[(df["Year"] >= selected[0]) & (df["Year"] <= selected[1])]
    traces = []
    for continent in dff.Continent.unique():
        traces.append(go.Box(y=dff[dff["Continent"] == continent]["Emission"],name=continent,marker={"size": 4}))
    return {"data": traces,
            "layout": go.Layout(title=f"Emission Levels for {'-'.join(str(i) for i in selected)}",autosize=True,
                                margin={"l": 200, "b": 100, "r": 200},xaxis={"showticklabels": False,},
                                yaxis={"title": f"Emissions (gigatonnes of CO2)","type": "log",},)}


