#!/bin/bash 

# Fetch gaas

mkdir -p tmp
if [ ! -d tmp/GAAS ]; then
    git clone git@github.com:NBISweden/gaas.git tmp/GAAS
else 
    pushd tmp/GAAS
    git pull
    popd
fi

