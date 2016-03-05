#!/bin/bash

rm -rf target
cp -R additions target
cd src

if which make >/dev/null; then
	make
else
	mingw32-make
fi

cd ..