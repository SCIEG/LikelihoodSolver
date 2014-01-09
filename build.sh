#!/bin/bash

rm -rf target
cp -R additions target
cd src
make
cd ..