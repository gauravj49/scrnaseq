import dash
import dash_core_components as dcc
import dash_html_components as html
import pandas as pd
import plotly.graph_objs as go
from dash.dependencies import Input, Output
import os
from app import app

if 'DYNO' in os.environ:
    app_name = os.environ['DASH_APP_NAME']
else:
    app_name = 'dash-mapplot'

mapbox_access_token = "pk.eyJ1IjoicHJpeWF0aGFyc2FuIiwiYSI6ImNqbGRyMGQ5YTBhcmkzcXF6YWZldnVvZXoifQ.sN7gyyHTIq1BSfHQRBZdHA"


df_airports = pd.read_csv(
    'https://raw.githubusercontent.com/plotly/datasets/master/2011_february_us_airport_traffic.csv')
df_flight_paths = pd.read_csv(
    'https://raw.githubusercontent.com/plotly/datasets/master/2011_february_aa_flight_paths.csv')

map_list = ["Scatter Plots on Mapbox", "Choropleth Maps", "Scatter Plots on Maps", "Bubble Maps", "Lines on Maps"]
layout = html.Div([html.Div([html.H1("Maps"), ], style={'textAlign': "center"}),
                       html.Div([html.Span("Type of map", className="six columns",
                                           style={'textAlign': "center", "display": "block",
                                                  "text-decoration": "underline", "padding-top": 5}),
                                 dcc.Dropdown(id="map-type", options=[{'label': i, 'value': i} for i in map_list],
                                              value='Scatter Plots on Mapbox', className="six columns")],
                                className="row",
                                style={"display": "block", "margin-left": "auto",
                                       "margin-right": "auto", "width": "80%", "padding-top": 10}),
                       html.Div([dcc.Graph(id="my-graph"), ], className="row"),
                       ], className="container")


@app.callback(
    Output("my-graph", "figure"),
    [Input("map-type", "value")])
def update_graph(type):
    trace1 = [go.Scattermapbox(lat=df_airports["lat"], lon=df_airports["long"], mode='markers', hoverinfo='text',
                               marker={'symbol': "airport", 'size': 8}, text=df_airports['airport'])]
    layout1 = go.Layout(title=f'Airport locations', autosize=True, hovermode='closest', showlegend=False, height=550,
                        mapbox={'accesstoken': mapbox_access_token, 'bearing': 0, 'center': {'lat': 38, 'lon': -94},
                                'pitch': 30, 'zoom': 3, "style": 'mapbox://styles/mapbox/light-v9'}, )
    df = df_airports.groupby(["state"])["cnt"].sum().reset_index()
    trace2 = [go.Choropleth(locationmode='USA-states', locations=df["state"].tolist(), z=df["cnt"].tolist(),
                            zmin=df["cnt"].min(), zmax=60000, text=df["state"], autocolorscale=False, colorscale='Hot',
                            colorbar={"title": "Flight Count", "thickness": 15},
                            marker={'line': {'color': 'rgb(180,180,180)',
                                             'width': 0.5}}, )]
    layout2 = go.Layout(title="Statewide count of flights for February 2011",
                        geo={"scope": 'usa', "projection": {"type": "albers usa"}})
    trace3 = go.Scattergeo(locationmode='USA-states', lon=df_airports['long'], lat=df_airports['lat'],
                           text=df_airports['airport'], showlegend=False, mode='markers',
                           marker={'size': 6, 'opacity': 0.8, 'symbol': "diamond", "color": "#c51b8a",
                                   'line': {'width': 0.1, 'color': 'rgba(102, 102, 102)'}})

    layout3 = {'title': 'Airport locations',
               'geo': {'scope': 'usa', 'projection': {'type': 'albers usa'}, 'showland': True,
                       'landcolor': "rgb(250, 250, 250)", 'subunitcolor': "rgb(217, 217, 217)",
                       'countrycolor': "rgb(217, 217, 217)", 'countrywidth': 0.5, 'subunitwidth': 0.5}}

    limits = [(0, 200), (201, 500), (501, 1000), (1001, 4000), (4001, 10000), (10001, 13000), (13001, 30000)]
    color = ['#2166ac', '#67a9cf', '#d1e5f0', '#f7f7f7', '#fddbc7', '#ef8a62', '#b2182b']
    trace4 = []
    for i in range(len(limits)):
        lim = limits[i]
        dff = df_airports[(df_airports["cnt"] >= lim[0]) & (df_airports["cnt"] <= lim[1])]
        trace4.append(go.Scattergeo(locationmode='USA-states', lon=dff['long'], lat=dff['lat'], text=dff["airport"],
                                    marker=go.scattergeo.Marker(size=dff['cnt'] / 10, color=color[i],
                                                                line={"width": 0.5, "color": 'rgb(40,40,40)'},
                                                                sizemode='area'),
                                    name='Flight Count :{0} - {1}'.format(lim[0], lim[1])))

    layout4 = go.Layout(title=go.layout.Title(text='No of Flights For February 2011'), showlegend=True,
                            geo=go.layout.Geo(scope='usa', projection=go.layout.geo.Projection(type='albers usa'), ))
    trace = []
    for i in range(len(df_flight_paths)):
        trace.append(
            go.Scattergeo(locationmode='USA-states',
                          lon=[df_flight_paths['start_lon'][i], df_flight_paths['end_lon'][i]],
                          lat=[df_flight_paths['start_lat'][i], df_flight_paths['end_lat'][i]], mode='lines',
                          showlegend=False,
                          line=go.scattergeo.Line(width=1, color='#7a0177', ),
                          opacity=float(df_flight_paths['cnt'][i]) / float(df_flight_paths['cnt'].max()), ))

    layout5 = go.Layout(title='Flight paths February 2011', geo=go.layout.Geo(scope="north america",
                                                                              projection=go.layout.geo.Projection(
                                                                                  type="orthographic"), showland=True,
                                                                              landcolor='rgb(243, 243, 243)',
                                                                              countrycolor='rgb(204, 204, 204)', ), )
    if type == "Scatter Plots on Mapbox":
        return {"data": trace1, "layout": layout1}
    elif type == "Choropleth Maps":
        return {"data": trace2, "layout": layout2}
    elif type == "Scatter Plots on Maps":
        return {"data": [trace3], "layout": layout3}
    elif type == "Bubble Maps":
        return {"data": trace4, "layout": layout4}
    else:
        return {"data": [trace3] + trace, "layout": layout5}


