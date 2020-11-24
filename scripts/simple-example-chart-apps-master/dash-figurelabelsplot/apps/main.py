#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import dash_core_components as dcc
import dash_html_components as html
import pandas as pd
import plotly.graph_objs as go
from app import app
import dash_daq as daq
from dash.dependencies import Input, Output


df1 = pd.read_csv('https://raw.githubusercontent.com/plotly/datasets/master/data.csv')
df = df1.iloc[0:50]

if 'DYNO' in os.environ:
    app_name = os.environ['DASH_APP_NAME']
else:
    app_name = 'dash-figurelabelsplot'

font = ["Arial", "Open Sans", "Balto", "Courier New", "PT Sans Narrow", "Times New Roman", "Comic Sans MS",
        "cursive"]


layout = html.Div([
    html.Div([html.H1("Financial Statistics by Age")], style={'textAlign': "center"}),
    html.Div([html.Div([html.Div([daq.ColorPicker(id='my-color-picker', label='Select Color : Marker & Title',
                                                  value={"hex": "#E6738D"}, )], className="row",
                                 style={"margin-left": "auto", "margin-right": "auto", "width": "100%", "padding": 20}),

                        html.Div([html.Span("Select Font Style : ", className="six columns",
                                            style={"width": "50%", "text-align": "right", "padding": 5}),
                                  dcc.Dropdown(id="select-font", options=[{'label': i, 'value': i} for i in font],
                                               value="Open Sans", placeholder="Select a font",
                                               style={"display": "block", "margin-left": "auto", "margin-right": "auto",
                                                      "width": "87%"},
                                               className="six columns")], className="row",
                                 style={"display": "block", "margin-left": "auto", "margin-right": "auto",
                                        "width": "80%", "padding": 20}),
                        html.Div([html.Span("Input Font Size : ", className="six columns",
                                            style={"width": "50%", "text-align": "right", "padding": 5}),
                                  dcc.Input(id='size-input', type='number', value=15, className="six columns",
                                            style={"margin-left": "auto", "margin-right": "auto", "width": "45%"}),
                                  ], className="row",
                                 style={"display": "block", "margin-left": "auto", "margin-right": "auto",
                                        "width": "80%", "padding": 20}),

                        ], className="six columns"),
              html.Div([html.Div([html.Div([html.Span("Category :", className="five columns",
                                                      style={"width": 100, "text-align": "right", "padding-top": 15}),
                                            html.Div(dcc.Dropdown(id="selected-type", value='MonthlyIncome',
                                                                  options=
                                                                  [{'label': "Past-due ( 30-59days )",
                                                                    'value': 'NumberOfTime30-59DaysPastDueNotWorse'},
                                                                   {'label': "Debt-ratio", 'value': 'DebtRatio'},
                                                                   {'label': "Income", 'value': 'MonthlyIncome'},
                                                                   {'label': "Open Credits/loans",
                                                                    'value': 'NumberOfOpenCreditLinesAndLoans'},
                                                                   {'label': "Past-due ( 90 days )",
                                                                    'value': 'NumberOfTimes90DaysLate'},
                                                                   {'label': "Real estate loans",
                                                                    'value': 'NumberRealEstateLoansOrLines'},
                                                                   {'label': "Past-due ( 60-89 days )",
                                                                    'value': 'NumberOfTime60-89DaysPastDueNotWorse'},
                                                                   {'label': "Dependents",
                                                                    'value': 'NumberOfDependents'}],

                                                                  ), className="seven columns",
                                                     style={"width": "60%", "margin": 0, "padding": 10})],
                                           className="row", style={"margin": 10, "padding": 10}),
                                  dcc.Graph(id="my-graph", className="row")], className="six columns"), ])])
], className="container")


@app.callback(
    Output('my-graph', 'figure'),
    [Input('selected-type', 'value'),
     Input("my-color-picker", 'value'),
     Input('select-font', 'value'),
     Input('size-input', 'value')])
def update_figure(selected_type, selected_color, selected_font, size):
    color = selected_color["hex"]
    dropdown = {
        'NumberOfTime30-59DaysPastDueNotWorse': "Past-due 30-59days (no of days)",
        'DebtRatio': "Debt-ratio",
        'MonthlyIncome': "Income (USD)",
        'NumberOfOpenCreditLinesAndLoans': "Open Credits/loans (number)",
        'NumberOfTimes90DaysLate': "Past-due 90 days (no of days)",
        'NumberRealEstateLoansOrLines': "Real estate loans (number)",
        'NumberOfTime60-89DaysPastDueNotWorse': "Past-due 60-89 days (no of days)",
        'NumberOfDependents': "Dependents (number)"}

    trace = go.Scatter(
        x=df["age"],
        y=df[selected_type],
        name=dropdown[selected_type],
        mode='markers',
        marker={'size': 10,
                "color": color,
                "opacity": 0.8,
                "line": {"color": "#4C5050", "width": 0.5}
                }, )

    return {
        "data": [trace],

        "layout": go.Layout(
            title={"text": f'{dropdown[selected_type]} Vs Age',
                   "font": {"family": selected_font,
                            "size": size + 4,
                            "color": color}
                   },
            xaxis={
                'title': 'Age (years)',
                'titlefont': {'family': selected_font,
                              "size": size,
                              "color": color}},
            yaxis={'title': f'{dropdown[selected_type]}',
                   'titlefont': {'family': selected_font,
                                 "size": size,
                                 "color": color}}
        )

    }