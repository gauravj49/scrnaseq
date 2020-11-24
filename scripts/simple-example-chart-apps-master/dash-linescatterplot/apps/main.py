#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  4 14:07:46 2019

@author: divyachandran
"""

import dash
import dash_core_components as dcc
import dash_html_components as html
import pandas as pd
import plotly.graph_objs as go
import os
from app import app

if 'DYNO' in os.environ:
    app_name = os.environ['DASH_APP_NAME']
else:
    app_name = 'dash-linescatterplot'

df = pd.read_csv('https://raw.githubusercontent.com/plotly/datasets/master/Wage%20Rigidity%20Dataset.csv')
df.dropna(inplace=True)
df['year'] = pd.DatetimeIndex(df['Date']).year

layout = html.Div([
    html.Div([html.H1("Employment Wage Rigidity")], style={"text-align": "center"}),
    html.Div(dcc.Graph(id="my-graph")),
    html.Div([dcc.RangeSlider(id='year-slider', min=1983, max=df['year'].max(),
                              marks={1983: '1983', 1990: '1990', 2000: '2000', 2003: '2003', 2005: '2005', 2008: '2008',
                                     2010: '2010', 2013: '2013', 2016: '2016'}, value=[2000, 2005])
              ], style={"margin": 20, "padding": 30})
], className="container")


@app.callback(
    dash.dependencies.Output('my-graph', 'figure'),
    [dash.dependencies.Input('year-slider', 'value')])
def update_figure(selected_year):
    pd.options.mode.chained_assignment = None  # default='SettingWithCopyWarning'
    dff = df[(df.year >= selected_year[0]) & (df.year <= selected_year[1])]
    dff['Date'] = pd.to_datetime(dff['Date']).dt.strftime('%y/%d')
    trace1 = go.Scatter(y=dff["Hourly workers"], x=dff["Date"], mode='lines+markers', marker={"size": 3.5},
                        name="Hourly")
    trace2 = go.Scatter(y=dff['Non-hourly workers'], x=dff["Date"], mode='markers', marker={"size": 3},
                        name="Non-Hourly")
    trace3 = go.Scatter(y=dff["High school"], x=dff["Date"], mode='lines', marker={"size": 2}, name="High school")
    trace4 = go.Scatter(y=dff["Construction"], x=dff["Date"], mode='lines+markers', marker={"size": 3.5},
                        name="Construction")
    trace5 = go.Scatter(y=dff["Finance"], x=dff["Date"], mode='lines', marker={"size": 2}, name="Finance")
    trace6 = go.Scatter(y=dff["Manufacturing"], x=dff["Date"], mode='markers', marker={"size": 3}, name="Manufacturing")
    data = [trace1, trace2, trace3, trace4, trace5, trace6]
    return {"data": data,
            "layout": go.Layout(title=f"Wage Rigidity for {'-'.join(str(i) for i in selected_year)}",
                                yaxis={"title": "% of Jobstayers With a Wage Change of Zero", "range": [0, 25],
                                       "tick0": 0, "dtick": 5}, xaxis={"title": "Year", "tickangle": 45}, )}


