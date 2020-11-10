#!/bin/bash

here=$(realpath "$(dirname "$0")")

mkdir -p /mnt/sequences/data
mkdir -p /mnt/sequences/postgres-data

docker run --rm --name webapollo -d -e CATALINA_OPTS="-Xmx14G" -e APOLLO_ADMIN_EMAIL=admin@nbis.se -e APOLLO_ADMIN_PASSWORD="$(cat "$here/apolloadminpassword")" -v /mnt/sequences/data/:/data -v /mnt/sequences/postgres-data/:/var/lib/postgresql -p 8888:8080 gmod/apollo:latest

