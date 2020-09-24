#!/bin/bash

docker run --name webapollo -d -e CATALINA_OPTS="-Xmx14G" -v /mnt/sequences/data/:/data -v /mnt/sequences/postgres-data/:/var/lib/postgresql -p 8888:8080 gmod/apollo:latest

