#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import dash
import dash_core_components as dcc
import dash_html_components as html
import pandas as pd
import plotly.graph_objs as go
from app import app
import plotly.figure_factory as ff


df = pd.read_csv("https://raw.githubusercontent.com/divyachandran-ds/dataset/master/sales.csv")

charts = ["Box Plot", "Error Bar", "Histogram", "2D Histogram", "Distplot", "Violin Plot"]

if 'DYNO' in os.environ:
    app_name = os.environ['DASH_APP_NAME']
else:
    app_name = 'dash-statisticalplot'

layout = html.Div([
    html.H1("Statistical Charts", style={"textAlign": "center"}),
    html.P("( Sales Data of a Retail Store on Black Friday )", style={"textAlign": "center", "padding-bottom": 30}),
    html.Div([html.Div([html.Span("Type Of Chart : ")], className="six columns",
                       style={"textAlign": "right", "padding-right": 30, "padding-top": 7}),
              html.Div([dcc.Dropdown(id='chart-type', options=[{'label': i, 'value': i} for i in charts],
                                     value="Distplot")], className="six columns",
                       style={"width": "40%", "margin-left": "auto", "margin-right": "auto", "display": "block"})],
             className="row", style={"width": "80%"}),
    html.Div([dcc.Graph(id='my-graph')], className="row")
], className="container")


@app.callback(
    dash.dependencies.Output('my-graph', 'figure'),
    [dash.dependencies.Input('chart-type', 'value')])
def update_graph(chart):
    trace1 = []
    for age in ['0-17', '18-25', '26-35', '36-45', '46-50', '51-55', '55+']:
        trace1.append(go.Box(
            y=df[df["Age Group"] == age]['Purchase'], name=age))
    layout1 = go.Layout(title="Purchase vs Age group", xaxis={"title": "Age Group"}, yaxis={"title": "Sales ($)"},
                        colorway=['#b2182b', '#ef8a62', '#fddbc7', '#E7E7E7', '#d1e5f0', '#67a9cf', '#2166ac'])
    df_bar = df.groupby(["Occupation"])['Purchase'].sum().reset_index()
    trace2 = [go.Bar(x=df_bar["Occupation"], y=df_bar["Purchase"],
                     marker={'color': '#ef8a62', 'line': {'color': "#b2182b", 'width': 0.5}}, opacity=0.6,
                     error_y={'type': 'percent', 'value': 10})]
    layout2 = go.Layout(title="Purchase vs Occupation Types", xaxis={"title": "Occupation Category (20 Types)"},
                        yaxis={"title": "Sales ($)"})
    trace3 = []
    for product in ['Product 1', 'Product 2', 'Product 3']:
        trace3.append(go.Histogram(x=df[product], name=product, xbins={"size": 3}, opacity=0.8))
    layout3 = go.Layout(title="Product Category Distribution", xaxis={"title": "Products"},
                        yaxis={"title": "Frequency"}, barmode='overlay', colorway=['#e9a3c9', '#ffffbf', '#a1d76a'])
    trace4 = [go.Histogram2d(x=df["Age Group"].sort_values(), y=df['Purchase'], histnorm='probability', autobinx=False,
                             xbins={"size": 1}, autobiny=False, ybins={"size": 1000},
                             colorscale=[[0, 'rgb(12,51,131)'], [0.25, 'rgb(10,136,186)'], [0.5, 'rgb(242,211,56)'],
                                         [0.75, 'rgb(242,143,56)'], [1, 'rgb(217,30,30)']])]
    layout4 = go.Layout(title="Sales vs Age Group distribution", xaxis={"title": "Age Group"},
                        yaxis={"title": "Sales ($)"})
    x1 = df[df["City"] == "A"]["Purchase"]
    x2 = df[df["City"] == "B"]["Purchase"]
    x3 = df[df["City"] == "C"]["Purchase"]
    hist_data = [x1, x2, x3]
    group_labels = ["City A", "City B", "City C"]
    figure5 = ff.create_distplot(hist_data, group_labels, bin_size=1000, colors=['#a6cee3', '#1f78b4', '#b2df8a'])
    figure5['layout'].update(title='Sales Distribution Over Cities', yaxis={"title": " Probability Density of Sales"})
    trace6 = [{"type": 'violin', "x": df[df['Gender'] == 'M']["Resident of the current city"],
               "y": df[df['Gender'] == 'M']["Purchase"], "legendgroup": 'M', "scalegroup": 'M', "name": 'Male',
               "box": {"visible": True}, "meanline": {"visible": True},
               "line": {"color": '#C1EC00'}},
              {"type": 'violin', "x": df[df['Gender'] == 'F']["Resident of the current city"],
               "y": df[df['Gender'] == 'F']["Purchase"], "legendgroup": 'F', "scalegroup": 'F', "name": 'Female',
               "box": {"visible": True}, "meanline": {"visible": True}, "line": {"color": '#EC7899'}}]
    layout6 = go.Layout(title="Sales Distribution Over Duration of Stay", xaxis={"title": " Duration of Stay (years)"},
                        yaxis={"title": "Sales ($)"}, violinmode="group")
    if chart == "Box Plot":
        return {"data": trace1, "layout": layout1}
    elif chart == "Error Bar":
        return {"data": trace2, "layout": layout2}
    elif chart == "Histogram":
        return {"data": trace3, "layout": layout3}
    elif chart == "2D Histogram":
        return {"data": trace4, "layout": layout4}
    elif chart == "Distplot":
        return figure5
    else:
        return {"data": trace6, "layout": layout6}