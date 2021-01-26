#!/bin/bash

LABEL="$1"
LABEL=${LABEL:-sequence}
LABEL=$(echo $LABEL | tr -c -d '[a-z0-9]')

mountpath=/mnt/$LABEL

if mount | grep -q "$mountpath" ; then
    exit 0
fi

sudo mkdir -p "$mountpath"
sudo mkfs.ext4 -L "$LABEL" /dev/vdb
echo "LABEL=$LABEL $mountpath ext4 defaults,noatime,_netdev,nofail 0 2" | sudo tee --append /etc/fstab
sudo mount -a
sudo chown "$USER" "$mountpath"
echo "$LABEL" > "$mountpath/species"
