import os

import dash
import dash_core_components as dcc
import dash_html_components as html
import pandas as pd
import plotly.graph_objs as go
from dash.dependencies import Input, Output
from app import app

if 'DYNO' in os.environ:
    app_name = os.environ['DASH_APP_NAME']
else:
    app_name = 'dash-candlestickplot'

df = pd.read_csv("https://raw.githubusercontent.com/plotly/datasets/master/tesla-stock-price.csv")
df['year'] = pd.DatetimeIndex(df['date']).year

layout = html.Div([
    html.Div([html.H1("Tesla Stock Price")], style={'textAlign': "center"}),
    dcc.Graph(id="my-graph"),
    html.Div([dcc.RangeSlider(id="select-range",marks={i: '{}'.format(i) for i in df.year.unique().tolist()},
                              min=df.year.min(),max=df.year.max(),value=[2016, 2017])], style={"padding-top": 100,})
], className="container")

@app.callback(
    Output("my-graph", 'figure'),
    [Input("select-range", 'value')])
def update_figure(selected):
    dff = df[(df["year"] >= selected[0]) & (df["year"] <= selected[1])]
    trace = go.Candlestick(x=dff['date'],open=dff['open'],high=dff['high'],low=dff['low'],close=dff['close'],
                           increasing={'line': {'color': '#00CC94'}},decreasing={'line': {'color': '#F50030'}})
    return {'data': [trace],
            'layout': go.Layout(title=f"Stock Values for the period:{'-'.join(str(i) for i in selected)}",
                                xaxis={'rangeslider': {'visible': False},'autorange': "reversed",},
                                yaxis={"title": f'Stock Price (USD)'})}
