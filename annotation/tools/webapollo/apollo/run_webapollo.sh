#!/bin/bash

# If given a single paramter, try to use that.

if [ "$#" -eq 1 ]; then
  if [ -d "$1" ] ; then
    instancedir=$(realpath "$1")
  else
    instancedir=/mnt/$1
  fi
else
  instancedir=/mnt/*
fi

external_port=8888

sysmem=$(sed -ne 's/^MemTotal:\s*\([0-9][0-9]*\)\s.*/\1/p' /proc/meminfo)
heapsize=$(echo "$sysmem" \* 8 / 10 | bc)

if [ -d "$instancedir" ]; then
  :
else
  echo "Not a single directory under /mnt/ or incorrect parameter, bailing out"
  exit 1
fi

if [ -r "$instancedir/external_port" ]; then
  external_port=$(cat "$instancedir/external_port")
fi

if [ -r "$instancedir/adminpassword" ]; then
  adminpass=$(cat "$instancedir/adminpassword")
else
  adminpass=$(openssl rand -base64 48 | tr  -c -d '[:alnum:]')
  echo "$adminpass" > "$instancedir/adminpassword"
fi

docker run --name webapollo -d -e CATALINA_OPTS="-Xmx${heapsize}k  -e APOLLO_ADMIN_EMAIL=admin@nbis.se -e APOLLO_ADMIN_PASSWORD=${adminpass}" -v "${instancedir}/data/:/data" -v "${instancedir}"/postgres-data/:/var/lib/postgresql -p "${external_port}:8080 gmod/apollo:latest

