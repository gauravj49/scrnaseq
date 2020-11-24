import re
from app import app
import os
import dash
import dash_core_components as dcc
import dash_daq as daq
import dash_html_components as html
import pandas as pd
import plotly.graph_objs as go

df = pd.read_csv('https://raw.githubusercontent.com/plotly/datasets/master/iris.csv')

if 'DYNO' in os.environ:
    app_name = os.environ['DASH_APP_NAME']
else:
    app_name = 'dash-legendplot'

layout = html.Div([html.Div([html.H1("Iris Species Comparison"), ], style={'textAlign': "center"}),
                       html.Div([html.Div([html.Span("x-Axis : ", className="six columns",
                                                     style={'textAlign': "right", "text-decoration": "underline"}),
                                           dcc.RadioItems(id="xaxis",
                                                          options=[{'label': 'Sepal Length', 'value': 'SepalLength'},
                                                                   {'label': 'Sepal Width', 'value': 'SepalWidth'}, ],
                                                          value='SepalLength',
                                                          labelStyle={"display": "inline", "padding": 10})],
                                          className="six columns"),

                                 html.Div([html.Span("y-Axis : ", className="six columns",
                                                     style={'textAlign': "right", "text-decoration": "underline"}),
                                           dcc.RadioItems(id="yaxis",
                                                          options=[{'label': 'Petal Length', 'value': 'PetalLength'},
                                                                   {'label': 'Petal Width', 'value': 'PetalWidth'}, ],
                                                          value='PetalLength',
                                                          labelStyle={"display": "inline", "padding": 10})],
                                          className="six columns")], className="row"),
                       html.Div([dcc.Graph(id="my-graph"), ], className="row"),
                       html.Div([html.Span("Legend Visible :  ", className="three columns",
                                           style={'textAlign': "left", "padding-top": 2}),
                                 daq.BooleanSwitch(id="legend", on=True, label=["Hide", "Show"], labelPosition="top",
                                                   color="#137d1c", className="two columns"),
                                 html.Span("Legend Orientation :", className="four columns",
                                           style={'textAlign': "left", "padding-top": 2}),
                                 daq.ToggleSwitch(id="position", value=True, label=["Horizontal", "Vertical"],
                                                  labelPosition="top", className="three columns")
                                 ], className="row",
                                style={"padding-top": 10, "margin-left": "auto", "margin-right": "auto",
                                       "width": "80%"}),

                       html.Div([html.P("Legend Position:Input x and y values", className="row",
                                        style={'textAlign': "center", "padding": 5,
                                               "margin-left": "auto", "margin-right": "auto", "width": "40%",
                                               "border": ".5px solid #ccc"}),
                                 html.Div([dcc.Input(id="x-value", type="number", value=1, style={"width": "30%"},
                                                     className="four columns"),
                                           dcc.Input(id="y-value", type="number", value=1, style={"width": "30%"},
                                                     className="four columns"),
                                           html.Button('Submit', id="button", className="four columns")],
                                          className="row",
                                          style={"margin-left": "auto", "margin-right": "auto", "width": "60%"})],
                                style={"padding-top": 30}),
                       ], className="container")


@app.callback(
    dash.dependencies.Output("my-graph", "figure"),
    [dash.dependencies.Input("xaxis", "value"), dash.dependencies.Input("yaxis", "value"),
     dash.dependencies.Input("legend", "on"), dash.dependencies.Input("position", "value"),
     dash.dependencies.Input("button", "n_clicks")],
    [dash.dependencies.State("x-value", "value"), dash.dependencies.State("y-value", "value")])
def update_graph(x_axis, y_axis, legend, position, n_clicks, xvalue, yvalue):
    trace = []
    for name in df.Name.unique():
        dff = df[df["Name"] == name]
        trace.append(go.Scatter(x=dff[x_axis], y=dff[y_axis], name=name, mode="markers", marker={'size': 10, }, ))
    return {"data": trace,
            "layout": go.Layout(title=f"Iris data",
                                xaxis={"title": f"{re.sub(r'(?<=[a-z])(?=[A-Z])', ' ', x_axis)} (cm)"},
                                yaxis={"title": f"{re.sub(r'(?<=[a-z])(?=[A-Z])', ' ', y_axis)} (cm)"},
                                colorway=["#E20048", "#CEF600", "#FFCB00"],
                                legend={"x": xvalue, "y": yvalue, "orientation": f'{"v" if position == True else "h"}'},
                                showlegend=legend)}


