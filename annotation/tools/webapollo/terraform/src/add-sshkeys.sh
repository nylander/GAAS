#!/bin/bash 

# Add SSH keys

if [ -d ~/authorized-keys -a ~/.ssh ];then
    cat ~/authorized-keys/*.pub >> ~/.ssh/authorized_keys
    chmod 600 ~/.ssh/authorized_keys
    rm -rf ~/authorized-keys
fi

