#!/bin/bash -e

# Exit if installed
if test -d "$HOME/conda" ; then
    exit 0
fi

cd /tmp
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash ./Miniconda3* -b -f -p "$HOME/conda"
. "$HOME/conda/bin/activate"
conda init
conda env create -f "$HOME/nbis/GAAS/annotation/tools/webapollo/apollo/gmod.yaml" || echo "env create returned non-zero"
rm -f ./Miniconda3*
