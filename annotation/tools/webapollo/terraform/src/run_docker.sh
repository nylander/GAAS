#!/bin/bash
usage="
USAGE: $0 --admin admin_username --password admin_password

"
admin_username=
admin_password=

if [ $# -lt 1 ]; then
    echo "$usage"
    exit 1
fi

isNonOptionArg=false
while [ "$1" != "" ]; do
    if [ "$isNonOptionArg" == "true" ]; then
        isNonOptionArg=false
    elif [ "$1" == "--" ]; then
        isNonOptionArg=true
    elif [ "${1:0:1}" == "-" ]; then
        case $1 in
            -h| --help) echo "$usage"; exit;;
            --admin) admin_username=$2;shift;;
            --password) admin_password=$2;shift;;
            -*) echo "Error! Wrong argument: $1">&2; exit;;
        esac
    fi
    shift
done

if [ "$admin_username" == "" -o "$admin_password" == "" ];then
    echo "$usage"
    exit 1
fi

pushd /media/storage
mkdir -p software
pushd software
git clone https://github.com/NBISweden/GAAS


sudo docker run --name webapollo -d -e CATALINA_OPTS=-Xmx14G -e APOLLO_ADMIN_EMAIL=$admin_username -e APOLLO_ADMIN_PASSWORD=$admin_password  -v /mnt/sequences/data/:/data -v /mnt/sequences/postgres-data/:/var/lib/postgresql -p 8888:8080 gmod/apollo:latest
