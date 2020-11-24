import os
import dash_core_components as dcc
import dash_html_components as html
import pandas as pd
import plotly.graph_objs as go
import dash_daq as daq
from dash.dependencies import Input, Output
from app import app

if 'DYNO' in os.environ:
    app_name = os.environ['DASH_APP_NAME']
else:
    app_name = 'dash-axesplot'

df = pd.read_csv("https://raw.githubusercontent.com/divyachandran-ds/Datascience/master/bchealth.csv")
font = ["Arial", "Open Sans", "Balto", "Courier New", "PT Sans Narrow", "Times New Roman", "Comic Sans MS",
        "cursive"]
layout = html.Div([
    html.Div([html.H1("British Columbia Healthcare Data")], style={"text-align": "center"}),
    html.Div([html.Div([html.P("Toggling Axes Lines, Ticks, Labels", className="row",
                        style={"text-align": "center", "border": "1px solid #ccc"}),
              html.Div([html.Span("Grid-Lines", className="four columns", style={"text-align": "right"}),
                       html.Div(daq.BooleanSwitch(id='grid-type',label=['Hide', 'Show'],labelPosition="bottom",on=True,
                                             ) , className="eight columns")], className="row", style={"padding": 5}),
              html.Div([html.Span("Tick-Labels", className="four columns", style={"text-align": "right"}),
                       html.Div(daq.BooleanSwitch(id='label-type',label=['Hide', 'Show'],labelPosition="bottom",on=True,
                                            ), className="eight columns")], className="row", style={"padding": 5}),
              html.Div([html.Span("Line", className="four columns", style={"text-align": "right"}),
                       html.Div(daq.BooleanSwitch(id='line-type',label=['Hide', 'Show'],labelPosition="bottom",on=True,
                       ),className="eight columns")], className="row", style={"padding": 5}),], className="row"),
              html.Div([html.Span("Expenses: y-axis range ", className="column",
                                  style={"text-align": "center", "border": "1px solid #ccc"}),
                        html.Div([dcc.RangeSlider(id='range1',min=0,max=100000000,step=500,updatemode="drag",value=[0, 20000000],
                                                  marks={10000000: "10M",20000000: "20M",30000000: "30M",40000000: "40M",
                                                         50000000: "50M",60000000: "60M",70000000: "70M",80000000: "80M",
                                                         90000000: "90M",100000000: "100M"})], className="column",
                       style={"margin": 0, "padding": 10})], className="row", style={"padding": 15}),
              html.Div([html.Span("Number of Practitioners: y-axis range ", className="column",
                                  style={"text-align": "center", "border": "1px solid #ccc"}),
                        html.Div([dcc.RangeSlider(id='range2',min=0,max=500,step=10,updatemode="drag",value=[2, 50],
                                                  marks={i * 50: str(i * 50) for i in range(0, 11)}),],
                                 className="column", style={"margin": 0, "padding": 10})], className="row",style={"padding": 15}),
              html.Div([html.Div([html.P("Tick : Color and Style", className="row",
                                         style={"text-align": "center", "border": "1px solid #ccc"}),
                        html.Div([html.Div(daq.ColorPicker(id='tick-color-picker',label='Tick Color',size=150,
                                                 value={"hex": "#CFD0E0"},))], className="row"),
                        html.Div([html.P("Tick Size: Length & Width", className="row",
                                         style={"text-align": "center", "border": "1px solid #ccc"}),
                                  html.Div([dcc.Input(id='length', type='number', value=10, className="six columns"),
                                            dcc.Input(id='width', type='number', value=10, className="six columns"),
                                  ], className="row", style={"padding": 3})],
                                 className="row", style={"padding-top": 15, "padding-bottom": 15}),], className="six columns",
                                       style={"width": "48%", "margin": 0, "float": "left", "padding": 3}),
                        html.Div([html.P("Axes : Color and Style", className="row",
                                         style={"text-align": "center", "border": "1px solid #ccc"}),
                                  html.Div([html.Div(daq.ColorPicker(id='axes-color-picker',label='Axes Color',
                                                                     size=150,value={"hex": "#080808"},))], className="row"),
                                  html.Div([html.P("Axes Font", className="row", style={"text-align": "center", "border": "1px solid #ccc"}),
                                            html.Div(dcc.Dropdown(id="axes-font",options=[{'label': i, 'value': i} for i in font],
                                                                  value="Open Sans",placeholder="Select a font",
                                                                  style={"padding-top": 3}), className="row")
                ], className="row", style={"padding-top": 15, "padding-bottom": 15, "padding-left": 10}),
            ], className="six columns", style={"width": "48%", "margin": 0, "float": "left", "padding": 3}),
        ], className="row", style={"padding": 5})],className="five columns", style={"width": "48%", "margin": 0, "float": "left"}),

    html.Div([html.Span("Health Service Delivery Area", className="row", style={"text-align": "center"}),
              dcc.Dropdown(id="select-hsda",options=[{'label': 'East Kootenay', 'value': '11 - East Kootenay'},
                              {'label': 'Kootenay/Boundary', 'value': '12 - Kootenay/Boundary'},
                              {'label': "Okanagan", 'value': '13 - Okanagan'},
                              {'label': 'Thompson/Cariboo', 'value': '14 - Thompson/Cariboo'},
                              {'label': 'Fraser East', 'value': '21 - Fraser East'},
                              {'label': 'Fraser North', 'value': '22 - Fraser North'},
                              {'label': 'Fraser South', 'value': '23 - Fraser South'},
                              {'label': 'Richmond', 'value': '31 - Richmond'},
                              {'label': "Vancouver", 'value': '32 - Vancouver'},
                              {'label': "North Shore/Coast Garibaldi", 'value': '33 - North Shore/Coast Garibaldi'},
                              {'label': "South Vancouver Island", 'value': '41 - South Vancouver Island'},
                              {'label': "Central Vancouver Island", 'value': '42 - Central Vancouver Island'},
                              {'label': "North Vancouver Island", 'value': '43 - North Vancouver Island'},
                              {'label': "Northwest", 'value': '51 - Northwest'},
                              {'label': "Northern Interior", 'value': '52 - Northern Interior'},
                              {'label': "Northeast", 'value': '53 - Northeast'}, ],value='31 - Richmond', ),
              dcc.Graph(id="my-graph"),], className="seven columns",
             style={"width": "48%", "margin": 0, "float": "right"})
], className="container")


@app.callback(
    Output('my-graph', 'figure'),
    [Input('grid-type', 'on'),Input('label-type', 'on'),Input('line-type', 'on'),Input("range1", 'value'),
     Input("range2", 'value'),Input('tick-color-picker', 'value'),Input('axes-color-picker', 'value'),
     Input('axes-font', 'value'),Input('select-hsda', 'value'),Input('length', 'value'),Input('width', 'value'),])
def update_figure(grid, label, line, range1, range2, tick_color, axes_color, axes_font, selected, len, width):
    color1 = tick_color["hex"]
    color2 = axes_color["hex"]
    dff = df[df["HSDA"] == selected]
    list = []
    for n in df['SPECIALTY'].unique():
        word = n.split('- ')[1]
        list.append(word.title())
    trace1 = go.Scatter(x=df['SPECIALTY'],y=dff['PAYMENTS'],mode='markers',name="Expenses",marker={"color": "#008D00",
                                                                                                   "size": 8})
    trace2 = go.Scatter(x=df['SPECIALTY'],y=dff['PRACTITIONERS'],mode='markers',marker={"color": "#FF6262","size": 8},
                        name="Practitioners",yaxis="y2")
    layout = go.Layout(title="Expenditure & Number of Practitioners vs Speciality",height=600,
                       legend={"y": 1.1, "orientation": "h"},
                       xaxis={"showgrid": grid,"showline": line,"showticklabels": label,"tickcolor": color1,
                              "ticklen": len - 4,"tickwidth": width - 4,
                              "title": {"text": f"Speciality","font": {"family": axes_font, "color": color2}},
                              "tickmode": "array","tickangle": 40,"tickvals": df['SPECIALTY'].unique(),"ticktext": list,
                              "tickfont": {"size": 7},},
                       yaxis={"showgrid": grid,"showline": line,"showticklabels": label,"tickcolor": color1,
                              "range": [range1[0], range1[1]],"ticklen": len,"tickwidth": width,
                              "title": {"text": f"Expenses (CAD)","font": {"family": axes_font, "color": color2}}},
                       yaxis2={"showgrid": grid,"showline": line,"showticklabels": label,"tickcolor": color1,
                               "range": [range2[0], range2[1]],"ticklen": len,"tickwidth": width,
                               "title": {"text": f"Number of Practitioners","font": {"family": axes_font, "color": color2}
                                },"overlaying": 'y',"side": 'right'})

    return {"data": [trace1, trace2],"layout": layout}