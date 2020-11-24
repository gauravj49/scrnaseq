#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import dash
import dash_core_components as dcc
import dash_html_components as html
import pandas as pd
import plotly.graph_objs as go
from app import app

df = pd.read_csv('https://raw.githubusercontent.com/plotly/datasets/master/Mining-BTC-180.csv')
df["month"] = pd.DatetimeIndex(df["Date"]).month

if 'DYNO' in os.environ:
    app_name = os.environ['DASH_APP_NAME']
else:
    app_name = 'dash-plotfromcsvplot'

layout = html.Div([
    html.Div([html.H1("Bitcoin Statistics Over Time")], style={'textAlign': "center"}),
    html.Div([dcc.Dropdown(id='value-selected',
                           options=[{'label': str(i).replace('-', ' '), 'value': i} for i in df.columns.values[2:9]],
                           value=["Number-transactions"],
                           multi=True,
                           style={"display": "block", "margin-left": "auto", "margin-right": "auto", "width": "80%"},
                           className="eight columns"),
              html.A("View the CSV dataset",
                     href="https://raw.githubusercontent.com/plotly/datasets/master/Mining-BTC-180.csv",
                     target="_blank", style={"width": "30%", "float": "right"}, className="four columns")],
             className="row"),
    dcc.Graph(id="my-graph"),
    html.Div([dcc.RangeSlider(id='month-selected', min=4, max=10, step=1,
                              marks={4: "April", 5: "May", 6: "June", 7: "July", 8: "August", 9: "September",
                                     10: "October"}, value=[5, 7])])
], className="container")


@app.callback(
    dash.dependencies.Output("my-graph", "figure"),
    [dash.dependencies.Input("month-selected", "value"), dash.dependencies.Input("value-selected", "value")])
def update_graph(selected1, selected2):
    dff = df[(df["month"] >= selected1[0]) & (df["month"] <= selected1[1])]
    new_list = []
    for i in selected2:
        new_list.append((str(i).replace('-', ' ')))
    trace = []
    for indicator in selected2:
        trace.append(go.Scatter(x=dff.Date, y=dff[indicator], name=indicator, mode="lines",
                                marker={'size': 15, 'line': {'width': 0.5, 'color': 'white'}}, ))

    return {"data": trace, "layout": go.Layout(title=f"{','.join(new_list[0:-1])} and {''.join(new_list[-1])} vs Date",
                                               xaxis={"title": "Date"}, yaxis={"title": f"Value"},
                                               colorway=["#C7037A", "#E20048", "#FFCB00", "#FF7C00", "#2F9609",
                                                         "#0E4770", "#A8AE0B"])}
