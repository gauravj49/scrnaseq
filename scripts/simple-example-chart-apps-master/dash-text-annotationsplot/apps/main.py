#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import dash
import dash_core_components as dcc
import dash_html_components as html
import pandas as pd
import plotly.graph_objs as go
from app import app
from dash.dependencies import Input, Output, State

df = pd.read_csv('https://raw.githubusercontent.com/plotly/datasets/master/styled-line.csv')
month = ['January', 'February', 'March', 'April', 'May', 'June', 'July', 'August', 'September', 'October', 'November',
         'December']

if 'DYNO' in os.environ:
    app_name = os.environ['DASH_APP_NAME']
else:
    app_name = 'dash-text-annotationsplot'

layout = html.Div([
    html.Div([html.H1("Monthly Temperature Highs and Lows")], style={'textAlign': "center"}),
    html.Div([html.Div([dcc.Dropdown(id='value-selected', options=[{"label": i, 'value': i} for i in df.columns[1:]],
                                     value=['High 2007', 'Low 2000'], multi=True)],
                       style={"display": "block", "margin-left": "auto", "margin-right": "auto", "width": "60%"})],
             className="row"),
    dcc.Graph(id="my-graph"),
    html.Div([html.H6("Text-Annotations", className="row",
                      style={"display": "block", "text-align": "center", "text-decoration": "underline"}),
              html.Div([dcc.Dropdown(id='x-input', options=[{"label": i, "value": i} for i in month],
                                     placeholder="Select month", value='', className="three columns"),
                        dcc.Input(id='y-input', type='number', placeholder="Input temperature", value='',
                                  className="three columns"),
                        dcc.Input(id='text-input', type='text', placeholder="Input text", value='',
                                  className="two columns"),
                        html.Button(id='submit-button', children="Submit", className="two columns"),
                        html.Button(id='remove-button', children="Remove", className="two columns"), ], className="row",
                       style={"display": "block", "margin-left": "auto", "margin-right": "auto", "width": "100%"})])
], className="container")


@app.callback(
    Output("my-graph", "figure"),
    [Input("value-selected", "value"), Input("remove-button", 'n_clicks'), Input('submit-button', 'n_clicks')],
    [State('x-input', 'value'), State('y-input', 'value'), State('text-input', 'value')])
def update_graph(selected, remove, n_clicks, x_value, y_value, text):
    dropdown = {"High 2014": "High Temperature in 2014", "Low 2014": "Low Temperature in 2014",
                "High 2007": "High Temperature in 2007", "Low 2007": "Low Temperature in 2007",
                "High 2000": "High Temperature in 2000", "Low 2000": "Low Temperature in 2000"}
    ctx = dash.callback_context
    trace = []
    for value in selected:
        trace.append(go.Scatter(x=df["Months"], y=df[value], mode="lines+markers",
                                marker={"opacity": 0.7, 'size': 5, 'line': {'width': 0.5, 'color': 'white'}},
                                name=dropdown[value]))
    layout = go.Layout(colorway=["#EF533B", "#EF963B", "#287D95", "#2CB154", "#8C299D", "#8DD936"],
                       title="Temperature Over the Months", yaxis={"title": f"Temperature (degrees F)"},
                       xaxis={"title": "Months"})
    figure = {"data": trace, "layout": layout}
    if ctx.triggered[0]['prop_id'].split('.')[0] == 'submit-button':
        layout.update({"annotations": [
            {'x': x_value.title(), 'y': y_value, 'xref': 'x', 'yref': 'y', 'text': text, 'showarrow': True,
             'align': 'center'}]})
        return figure
    else:
        return figure