variable floating_ip_pool {
    description = "Name of the floating IP pool"
    default = "Public External IPv4 Network"
}
variable internal_network_name {
    description = "Name of the internal network"
    default = "SNIC 2020/20-30 Internal IPv4 Network"
}
variable ssh_user {
    description = "SSH user name (use the default user for the OS image)"
    default     = "ubuntu"
}
variable private_ssh_key {
    description = "Local path to the private SSH key"
    default =  "./private/webap.key"
}
variable public_ssh_key {
    description = "Local path to the public SSH key"
    default =  "./private/webap.key.pub"
}
variable image_name {
    description = "Name of the image for the instance"
    default =  "Ubuntu 18.04 LTS (Bionic Beaver) - latest"
}
variable flavor_name {
    description = "Name of the image for the instance"
    default =  "ssc.large"
}
variable instance_name {
    description = "Name for the instance"
    default = "webapollo"
}
variable keypair_name {
    description = "Name of keypair"
    default = "key-webap"
}
variable secgroup_name {
    description = "Name of seurity group"
    default = "secgroup-webap"
}
