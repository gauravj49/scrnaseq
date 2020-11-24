#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import dash
from dash.dependencies import Input, Output
import dash_core_components as dcc
import dash_html_components as html
import pandas as pd

from app import app

df = pd.read_excel('http://www.ams.usda.gov/sites/default/files/media/GTRFigure8.xlsx')


layout = html.Div([
    html.H3('Barge Freight Rates on the Illinois River', style={'textAlign': 'center'}),
    html.Div(
        id='graph-controls',
        className='row',
        children=[
            html.Div('Port City Groups: ', className=' offset-by-two two columns', style={'marginTop': 5}),
            dcc.Dropdown(
                id='graph-dropdown',
                className='six columns',
                options=[{'value': col, 'label': col}
                         for col in ['TWC', 'MM', 'ILL', 'ST LOUIS', 'CINC', 'LOH', 'CAR-MEM', 'MEM-SO']],
                value=['TWC', 'MM', 'ILL'],
                multi=True
            )
        ]
    ),
    dcc.Graph(id='graph')
], className="container")


@app.callback(
    Output('graph', 'figure'),
    [Input('graph-dropdown', 'value')])
def update_graph(values):
    data = [dict(
        type='scatter',
        mode='lines',
        x=df.DATE,
        y=df[col],
        name=col
    ) for col in values]

    layout = dict(
        xaxis=dict(title='Date'),
        yaxis=dict(title='Freight Rate'),
        margin=dict(t=40),
        annotations=[dict(
            text='<a href="https://catalog.data.gov/dataset/grain-transportation-report-figure-8">Source Data</a>',
            x=0,
            y=1.1,
            xref='paper',
            yref='paper',
            showarrow=False
        )]
    )
    return dict(data=data, layout=layout)
