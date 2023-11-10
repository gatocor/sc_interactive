#!/bin/bash

var=0.1

cp -r ../src/* ./app
cp -r ../requirements.txt ./
docker build . -t scinteractive:$var
#docker push dsblab/single_cell_analysis:$var
rm -rf ./app/*
rm requirements.txt
touch ./app/__init__.py
