#!/bin/bash 

# Retrieve common SSH keys from the repo

mkdir -p tmp
if [ ! -d tmp/annotation-cluster ]; then
    git clone git@github.com:NBISweden/annotation-cluster.git tmp/annotation-cluster
else 
    pushd tmp/annotation-cluster
    git pull
    popd
fi

