#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  1 11:02:59 2019

@author: divyachandran
"""

import dash

app = dash.Dash(__name__)
server = app.server
app.config.suppress_callback_exceptions = True