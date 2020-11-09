#!/bin/bash
pushd /media/storage
mkdir -p software
pushd software
git clone https://github.com/NBISweden/GAAS

sudo docker run --name webapollo -d -e CATALINA_OPTS=-Xmx14G -e APOLLO_ADMIN_EMAIL=adminuser -e APOLLO_ADMIN_PASSWORD=gW/aGeajRShQ3g  -v /mnt/sequences/data/:/data -v /mnt/sequences/postgres-data/:/var/lib/postgresql -p 8888:8080 gmod/apollo:latest
