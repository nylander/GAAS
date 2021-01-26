variable external_network_name {
    description = "Name of the external network"
    default = "Public External IPv4 Network"
}
variable internal_network_name {
    description = "Name of the internal network"
    default = "net_webapollo"
}
variable internal_subnet_name {
    description = "Name of the internal subnet"
    default = "subnet_webapollo"
}
variable port_name {
    description = "Name of the port"
    default = "port_webapollo"
}
variable router_name {
    description = "Name of the router_name"
    default = "router_webapollo"
}
variable ssh_user {
    description = "SSH user name (use the default user for the OS image)"
    default     = "ubuntu"
}
variable private_ssh_key {
    description = "Local path to the private SSH key"
    default =  "./private/webapollo.key"
}
variable public_ssh_key {
    description = "Local path to the public SSH key"
    default =  "./private/webapollo.key.pub"
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
    default = "key-webapollo"
}
variable secgroup_name {
    description = "Name of seurity group"
    default = "secgroup-webapollo"
}
variable size_volume {
    description = "size of the volume to be attached to the instance"
    default = 100
}
variable git_common_keys {
    description = "git repo containing common ssh keys"
    default = "git@github.com:NBISweden/annotation-cluster.git"
}
variable admin_username {
    description = "Name of the admin user to webapollo"
    default = "admin@nbis.se"
}
variable admin_password {
    description = "Password of the admin user to webapollo"
}
variable project_suffix {
    description = "project suffix added to all resource names"
    default = ""
}
