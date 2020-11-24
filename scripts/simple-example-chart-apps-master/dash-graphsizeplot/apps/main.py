import dash_core_components as dcc
import dash_html_components as html
import pandas as pd
import plotly.graph_objs as go
from dash.dependencies import Input, Output
from app import app
import os

df1 = pd.read_csv("https://raw.githubusercontent.com/plotly/datasets/master/Aids%20Data.csv")
df = df1.groupby(['Indicator', 'Time Period']).mean().reset_index()

if 'DYNO' in os.environ:
    app_name = os.environ['DASH_APP_NAME']
else:
    app_name = 'dash-graphsizeplot'

layout = html.Div([html.Div([html.H1("AIDS statistics over time")], style={'textAlign': "center"}),
                       html.Div([dcc.Dropdown(id="selected-type", multi=True, value=['AIDS-related deaths'],
                                              options=[
                                                  {"label": "Death", "value": 'AIDS-related deaths'},
                                                  {"label": "Receiving ART",
                                                   "value": 'Coverage of people receiving ART'},
                                                  {"label": "Death averted by ART",
                                                   "value": 'Deaths averted due to ART'},
                                                  {"label": "New infections", "value": 'New HIV Infections'},
                                                  {"label": "Pregnant women needing ART",
                                                   "value": 'Pregnant women needing ARV for preventing MTCT'},
                                                  {"label": "Pregnant women received ART",
                                                   "value": 'Pregnant women who received ARV for preventing MTCT'}
                                              ], )],
                                style={"display": "block", "margin-left": "auto", "margin-right": "auto",
                                       "width": "60%", }),
                       html.Div([dcc.Graph(id="my-graph",
                                           style={"display": "block", "margin-left": "auto", "margin-right": "auto",
                                                  "width": "100%"})]),
                       html.Div([html.Span("Slide to change graphsize", style={"text-align": "center", 'padding': 10},
                                           className="row"),
                                 dcc.Slider(id="graph-size", min=200, max=900, value=450, step=5, updatemode="drag",
                                            marks={200: "200px", 300: "300px", 400: "400px", 500: "500px", 600: "600px",
                                                   700: "700px", 800: "800px", 900: "900px", },
                                            className="row")],
                                style={"display": "block", "margin-left": "auto", "margin-right": "auto",
                                       "width": "40%", "padding": 20})
                       ], className="container")


@app.callback(
    Output('my-graph', 'figure'),
    [Input('selected-type', 'value'), Input('graph-size', 'value')])
def update_figure(selected, size_selected):
    traces = []
    for select in selected:
        traces.append(
            go.Scatter(x=df[df['Indicator'] == select]['Time Period'], y=df[df['Indicator'] == select]['Data Value'],
                       mode='lines+markers', name=select, hoverinfo="name", hoverlabel={"namelength": 50},
                       marker={"size": 6, }, showlegend=False))
    layout = go.Layout(title=f"Cases vs time", xaxis={"title": "Date"}, height=size_selected, width=2 * size_selected,
                       yaxis={"title": f"Number of cases", "tickangle": -45, "tickfont": {"size": 10}}, autosize=True,
                       colorway=["#003f5c", "#955196", "#dd5182", "#ff6e54", "#ffa600", "#061460"])

    figure = {"data": traces, "layout": layout}
    return figure


