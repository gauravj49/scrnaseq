import dash_core_components as dcc
import dash_html_components as html
import pandas as pd
import plotly.graph_objs as go
from dash.dependencies import Input, Output
from app import app
import os

if 'DYNO' in os.environ:
    app_name = os.environ['DASH_APP_NAME']
else:
    app_name = 'dash-basicplot'

df = pd.read_csv('https://raw.githubusercontent.com/divyachandran-ds/dataset/master/diabetes.csv')

layout = html.Div([
    html.Div([html.H1("Basic Charts")], className="row", style={'textAlign':"center", "padding": 5,"margin-left":"auto",
                                                                "margin-right": "auto", "width": "60%"}),
    html.Div([html.Span("Type Of Plot", style={'textAlign': "center", "display": "block"}),
             dcc.Dropdown(id='plot-selected',options=[{"label": i, 'value': i}
                          for i in ["Scatter Plot", "Line plot", "Bar plot", "Horizontal Bar", "Pie plot", "Table","Filled Area Plot"]],
                          value="Scatter Plot",)],className="row", style={"display": "block","margin-left": "auto",
                                                                          "margin-right": "auto","width": "40%","padding": 20}),
    html.Div([html.Div([html.Span("X-Axis : Scatter,Horizontal Barplot,Filledarea,Table",
                                  style={'textAlign': "center","display": "block"}),
                       dcc.Dropdown(id='xaxis-selected',options=[{"label": i, 'value': i} for i in df.columns[0:8]],
                                    value='Age',)], className="six columns"),
             html.Div([html.Span("Y-Axis : Scatter,Line,Barplot,Table", style={'textAlign': "center", "display": "block"}),
                       dcc.Dropdown(id='yaxis-selected',options=[{"label": i, 'value': i} for i in df.columns[0:8]],
                                    value='Glucose',)], className="six columns"),],className="row",style={"padding": 20}),
    html.Div([dcc.Graph(id="my-graph"),], className="row", style={"padding": 30}),
], className="container")


@app.callback(
    Output("my-graph", "figure"),
    [Input("xaxis-selected", "value"),Input("yaxis-selected", "value"),Input("plot-selected", "value")]

)
def update_graph(x_axis, y_axis, plot):
    text = {'Pregnancies':"Pregnancies",'Glucose':"Glucose (mmol/L)",'BloodPressure':"Blood Pressure (mm Hg) ",
            'SkinThickness':"Skin Thickness (mm)",'Insulin':"Insulin (mu U/ml)",'BMI':"BMI (weight in kg/(height in m)^2)",
            'DiabetesPedigreeFunction': "Diabetes Pedigree Function",'Age': "Age (years) ",'Outcome': "Outcome"}

    trace1 = go.Scatter(x=df[x_axis],y=df[y_axis],mode="markers",marker={"color": "#FF5757","opacity": 0.7,'size': 7,},)
    layout1 = go.Layout(title=f"{plot}",xaxis={"title": f"{text[x_axis]}"},yaxis={"title": f"{text[y_axis]}"})
    trace2 = go.Scatter(x=df["Age"].sort_values(ascending=True),y=df[y_axis],mode="lines",
                        marker={"color": "#0B660B","opacity": 0.7,'size': 5,'line': {'width': 0.5, 'color': 'white'}},)
    layout2 = go.Layout(title=f"{plot}",xaxis={"title": "Age (years)"},yaxis={"title": f"{text[y_axis]}"})
    trace3 = go.Bar(x=df["Outcome"].unique(),y=df[y_axis],width=0.3,marker={"color": "#660B3E"})
    layout3 = go.Layout(title=f"{plot}",xaxis={"title": "Diabetes Outcome"},yaxis={"title": f"{text[y_axis]}"})
    trace4 = go.Bar(x=df[x_axis],y=df["Outcome"].unique(),orientation='h',marker={"color": "#54770D"},width=0.3,)
    layout4 = go.Layout(title=f"{plot}",xaxis={"title": f"{text[x_axis]}"},yaxis={"title": f"Diabetes Outcome"})
    list = []
    for i in df.columns[0:8]:
        list.append(df[i].mean())
    trace5 = go.Pie(labels=df.columns[0:8],values=list)
    layout5 = go.Layout(title=f"{plot}")
    trace6 = go.Table(header={"values": [x_axis, y_axis],"fill": {"color": "#008366"},"align": ['left', 'center'],},
                      cells={"values":[df[x_axis], df[y_axis]],"fill":{"color": "#84FEE3"},"align":['left', 'center'],})
    layout6 = go.Layout(title=f"{plot}")
    trace7 = go.Scatter(x=df["Age"].sort_values(ascending=True),y=df[x_axis],marker={"color": "#EF2D9B"},fill='tozeroy',
                        fillcolor="#FF6BBF")
    layout7 = go.Layout(title=f"{plot}",xaxis={"title": "Age (years)"},yaxis={"title": f"{text[x_axis]}"})

    if plot == "Scatter Plot":
        return {"data": [trace1],"layout": layout1}
    elif plot == "Line plot":
        return {"data": [trace2],"layout": layout2}
    elif plot == "Bar plot":
        return {"data": [trace3],"layout": layout3}
    elif plot == "Horizontal Bar":
        return {"data": [trace4],"layout": layout4}
    elif plot == "Pie plot":
        return {"data": [trace5],"layout": layout5}
    elif plot == "Table":
        return {"data": [trace6],"layout": layout6}
    else:
        return {"data": [trace7],"layout": layout7}
