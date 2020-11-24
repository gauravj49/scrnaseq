#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import dash
import dash_core_components as dcc
import dash_html_components as html
import pandas as pd
import plotly.graph_objs as go
from app import app

df = pd.read_csv('https://raw.githubusercontent.com/plotly/datasets/master/streamtube-basic.csv')
df_1 = pd.read_csv("https://raw.githubusercontent.com/divyachandran-ds/dataset/master/Energy2.csv")
df_ribbon = pd.read_csv("https://raw.githubusercontent.com/plotly/datasets/master/spectral.csv")
charts = ["3D Scatter", "3D Surface", "3D Filled Line ", "Ribbon", "3D Streamtube"]

if 'DYNO' in os.environ:
    app_name = os.environ['DASH_APP_NAME']
else:
    app_name = 'dash-3dplot'

layout = html.Div([
    html.H1("3D Charts", style={"textAlign": "center"}),
    html.Div([html.Div([html.Span("Type Of Chart : ")], className="six columns",
                       style={"textAlign": "right", "padding-right": 30, "padding-top": 7}),
              html.Div([dcc.Dropdown(id='chart-type', options=[{'label': i, 'value': i} for i in charts],
                                     value="3D Scatter")], className="six columns",
                       style={"width": "40%", "margin-left": "auto", "margin-right": "auto", "display": "block"}),
              ], className="row", style={"width": "80%"}),
    html.Div([dcc.Graph(id='my-graph')], className="row")
], className="container")

@app.callback(
    dash.dependencies.Output('my-graph', 'figure'),
    [dash.dependencies.Input('chart-type', 'value')])
def update_graph(chart):
    df_scatter = df_1[df_1["Indicator Name"] == "CO2 emissions from gaseous fuel consumption (% of total)"]
    trace1 = [go.Scatter3d(
        x=df_scatter[df_scatter["Country Name"] == 'East Asia & Pacific'].values[0][3:].tolist(),
        y=df_scatter[df_scatter["Country Name"] == 'Europe & Central Asia'].values[0][3:].tolist(),
        z=df_scatter[df_scatter["Country Name"] == 'North America'].values[0][3:].tolist(),
        mode="markers", marker={"opacity": 0.8, "color": "#4daf4a", "line": {"width": 0.5, "color": "#056B00"}})]
    layout1 = go.Layout(
        title="CO2 Emissions From Gaseous Fuel Consumption (% of total)", height=800, width=800,
        scene={"camera": {"eye": {"x": 1.3, "y": 1.3, "z": 1.3}}, "aspectmode": "cube",
               "xaxis": {"title": "East Asia & Pacific"}, "yaxis": {"title": "Europe & Central Asia"},
               "zaxis": {"title": "North America"}}, )
    trace2 = [go.Surface(
        x=df_1["Country Name"].unique(), y=df_1["Indicator Name"].unique(), z=df_1.iloc[0:, 3:].values,
        colorscale='Blackbody', colorbar={"thickness": 10, "len": 0.5, 'title': {"text": "Indicator Value"}})]
    layout2 = go.Layout(
        title="Energy Indicators", height=1000, width=1000,
        scene={"camera": {"eye": {"x": 1.45, "y": 1.45, "z": 1.45}}, "aspectmode": "cube",
               "yaxis": {"tickangle": -30, "tickmode": "array",
                         "tickvals": df_1["Indicator Name"].unique(),
                         "ticktext": ['CO2 emissions<br>from gaseous fuel<br>consumption<br>(% of total)',
                                      'Access to<br>electricity<br>(% of population)',
                                      'Total<br>natural resources<br>rents (% of GDP)',
                                      'Combustible<br>renewables<br>and waste<br>(% of total energy)',
                                      'Electric power transmission<br>and distribution losses<br>(% of output)'],
                         "tickfont": {"size": 10}, },
               "xaxis": {"title": "Country", "tickmode": "array", "tickfont": {"size": 10},
                         "tickvals": df_1['Country Name'].unique(),
                         "ticktext": ['Afghanistan', 'Arab World', 'Australia', 'Belgium', 'Bangladesh', 'Brazil',
                                      'Canada',
                                      'Colombia', 'Germany', 'East Asia<br>& Pacific', 'Europe &<br>Central Asia',
                                      'India', 'Japan', 'Latin America &<br>Caribbean', 'Middle East &<br>North Africa',
                                      'Mexico', 'North America', 'Saudi Arabia', 'Singapore',
                                      'Virgin Islands (US)', 'South Africa', 'Zimbabwe'], },
               "zaxis": {"title": "Indicator Value<br>From 2000-2014"}},
    )

    list_country = ['East Asia & Pacific', 'Europe & Central Asia', 'Middle East & North Africa', 'North America',
                    'Arab World', 'Australia', 'Latin America & Caribbean']
    fill_color = ['#d73027', '#fc8d59', '#fee090', '#ffffbf', '#e0f3f8', '#91bfdb', '#4575b4']
    trace3 = []
    for country, fill in zip(list_country, fill_color):
        df_fill = df_1[df_1["Indicator Name"] == 'Electric power transmission and distribution losses (% of output)']
        dff = df_fill.groupby("Country Name").mean().reset_index()
        df_country = dff[dff["Country Name"] == country]
        years = df_country.columns[2:].tolist()
        zeros = [0] * len(years)
        z_value = df_country[df_country["Country Name"] == country].values[0][1:].tolist()
        country_list = [country] * len(years)
        trace3.append(go.Scatter3d(
            x=years + years[::-1] + [years[0]], y=country_list * 2 + [country_list[0]],
            z=z_value + zeros + [z_value[0]],
            mode="lines", surfaceaxis=1, surfacecolor=fill, line={'color': fill, 'width': 4}, ))
    layout3 = go.Layout(
        title="Electric power transmission and distribution losses per country", height=900, width=900,
        showlegend=False,
        scene={"camera": {"eye": {"x": 1.25, "y": -1.25, "z": .5}, "projection": {"type": "orthographic"}},
               "aspectmode": "cube",
               "yaxis": {"title": "", "tickmode": "array", "tickfont": {"size": 11}, "tickvals": list_country,
                         "tickangle": 35, "ticktext": ['East Asia &<br>Pacific', 'Europe &<br>Central Asia',
                                                       'Middle East &<br>North Africa', 'North America',
                                                       'Arab World', 'Australia', 'Latin America &<br>Caribbean'], },
               "xaxis": {"title": "", "tickfont": {"size": 11}, },
               "zaxis": {"title": "Electric power transmission and<br>     distribution losses<br>       (% of output)",
                         "tickfont": {"size": 10}, }, })
    trace5 = [go.Streamtube(
        x=df['x'].values.tolist(), y=df['y'].values.tolist(), z=df['z'].values.tolist(),
        u=df['u'].values.tolist(), v=df['v'].values.tolist(), w=df['w'].values.tolist(),
        sizeref=0.5, colorscale='Viridis', cmin=0, cmax=3)]
    layout5 = go.Layout(
        scene={'camera': {'eye': {'x': -0.5, 'y': 2, 'z': 0.5}}})
    trace4 = []
    y_raw = df_ribbon.iloc[:, 0]  # wavelength
    sample_size = df_ribbon.shape[1] - 1
    for i in range(1, sample_size):
        z_raw = df_ribbon.iloc[:, i]
        x = []
        y = []
        z = []
        for j in range(0, len(z_raw)):
            z.append([z_raw[j], z_raw[j]])
            y.append([y_raw[j], y_raw[j]])
            x.append([i * 2, i * 2 + 1])
        trace4.append({'z': z, 'x': x, 'y': y, 'colorscale': "Jet", 'showscale': False, 'type': 'surface'})

    layout4 = go.Layout(
        title="Spectral Data", scene={'camera': {'eye': {'x': -0.5, 'y': 2, 'z': 0.5}}})

    if chart == "3D Scatter":
        return {"data": trace1, "layout": layout1}
    elif chart == "3D Surface":
        return {"data": trace2, "layout": layout2}
    elif chart == "3D Filled Line ":
        return {"data": trace3, "layout": layout3}
    elif chart == "Ribbon":
        return {"data": trace4, "layout": layout4}
    else:
        return {"data": trace5, "layout": layout5}
