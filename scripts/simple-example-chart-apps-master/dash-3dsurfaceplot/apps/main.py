import os

import dash
import dash_core_components as dcc
import dash_html_components as html
import pandas as pd
import plotly.graph_objs as go
import us

from app import app

if 'DYNO' in os.environ:
    app_name = os.environ['DASH_APP_NAME']
else:
    app_name = 'dash-3dsurfaceplot'

df = pd.read_csv('https://raw.githubusercontent.com/plotly/datasets/master/laucnty16.csv', thousands=',')
layout = html.Div([
    html.Div([html.H1("Employment Statistics in the United States")], style={'textAlign': "center", }),
    html.Div([dcc.Dropdown(id="state-selected", options=[
        {'label': f'{us.states.lookup(str(i))}', 'value': i} for i in df['State FIPS Code'].unique()[9:]], value=45,
                           style={"display": "block", "margin-left": "auto", "margin-right": "auto",
                                  "width": "60%", })]),
    html.Div(dcc.Graph(id="my-graph")),
], className="container")


@app.callback(
    dash.dependencies.Output("my-graph", "figure"),
    [dash.dependencies.Input("state-selected", "value")]
)
def update_figure(selected):
    pd.options.mode.chained_assignment = None
    dff = df[df['State FIPS Code'] == selected]
    dff["county"] = (dff["County Name/State Abbreviation"].str.split(",", expand=True))[0]
    df1 = dff.loc[:, ["county", 'Labor Force', 'Employed', 'Unemployed']]
    df1.loc['Labor Force'] = pd.to_numeric(df1['Labor Force'], errors='ignore')
    df1.loc['Employed'] = pd.to_numeric(df1['Employed'], errors='ignore')
    df1.loc['Unemployed'] = pd.to_numeric(df1['Unemployed'], errors='ignore')
    trace = [go.Surface(y=df1.county.values, x=df1.Employed.values, z=df1.values, colorscale="YlGnBu", opacity=0.8,
                        colorbar={"title": "Number", "len": 0.5, "thickness": 15}, )]
    fig = go.Figure(data=trace,
                    layout=go.Layout(title=f'Annual Average of Labor Force Data for {us.states.lookup(str(selected))}',
                                     autosize=True, height=800,
                                     scene={"xaxis": {'title': "Annual Average of Employed  (number)",
                                                      "tickfont": {"size": 10}, 'type': "linear"},
                                            "yaxis": {"title": f"County in {us.states.lookup(str(selected))} ",
                                                      "tickfont": {"size": 10}, "tickangle": 1},
                                            "zaxis": {
                                                'title': "         Annual Average of : <br>Labour Force,Employed,Unemployed  ",
                                                "tickfont": {"size": 10}},
                                            "camera": {"eye": {"x": 2, "y": 1, "z": 1.25}}, "aspectmode": "cube", }))
    return fig
