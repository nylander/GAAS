#!/bin/bash

mountpath=/media/storage
sudo mkdir $mountpath
sudo mkfs.ext4 /dev/vdb
echo "/dev/vdb  $mountpath ext4 defaults,noatime,_netdev,nofail 0 2" | sudo tee --append /etc/fstab
sudo mount -a
sudo chown $USER $mountpath
