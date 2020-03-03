#!/bin/bash

if [[ $# != 1 ]];then
    echo "usage example: ./script Bio::Seq"
    exit
fi

perl -M${1} -le "\$mname=\"${1}.pm\";\$mname=~s#::#/#g;print \"$1 INSTALLED AT \$INC{\$mname}\";" 2>/dev/null || echo "${1} NOT INSTALLED"
