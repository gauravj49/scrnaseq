#!/bin/bash

set -e

cd $1
rm -rf .git
git init
git add *
git commit -am "commit"
git remote add plotly dokku@dash-simple-apps.plotly.host:$1
git push -f plotly master
rm -rf .git
cd ..
