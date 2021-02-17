#!/bin/bash

cd nbis/GAAS || exit 1

git status --porcelain=1 | grep '^ M ' | while read -r mod path; do
    git checkout "$path"
done
