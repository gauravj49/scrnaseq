import os
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output
from app import app
from apps import main, code

if 'DYNO' in os.environ:
    app_name = os.environ['DASH_APP_NAME']
else:
    app_name = 'dash-animationplot/'

app.layout = html.Div([
    dcc.Location(id='url', refresh=False),
    html.Div(id='page-content', children=main.layout)
])


@app.callback(Output('page-content', 'children'),
              [Input('url', 'pathname')])
def display_page(pathname):
    if pathname is None or pathname.replace(app_name, '').strip('/') == '':
        return main.layout
    elif pathname.replace(app_name, '').strip('/') == 'code':
        return code.layout
    else:
        return main.layout


if __name__ == '__main__':
    app.run_server(debug=True)
