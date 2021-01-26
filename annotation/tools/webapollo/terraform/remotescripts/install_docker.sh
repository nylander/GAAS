#!/bin/bash 

# Create and join the docker group.
sudo addgroup --system docker
sudo adduser "$(id -u -n)" docker

# install docker from snap
sudo snap install docker

sudo snap start docker
exit 0
