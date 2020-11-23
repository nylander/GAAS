resource "openstack_compute_keypair_v2" "my-cloud-key" {
  name       = var.keypair_name
  public_key =  file(var.public_ssh_key)

}

resource "openstack_networking_network_v2" "net_webap" {
  name           = var.internal_network_name
  admin_state_up = "true"
}

resource "openstack_networking_subnet_v2" "subnet_webap" {
  name       = var.internal_subnet_name
  network_id = openstack_networking_network_v2.net_webap.id
  cidr       = "10.0.0.0/24"
  ip_version = 4
  dns_nameservers = ["8.8.8.8"]
  enable_dhcp = "true"
}

resource "openstack_networking_router_v2" "router_webap" {
  name                = var.router_name
  admin_state_up      = "true"
  external_network_id = var.external_network_id
}

resource "openstack_networking_router_interface_v2" "router_interface" {
  router_id = openstack_networking_router_v2.router_webap.id
  subnet_id = openstack_networking_subnet_v2.subnet_webap.id
}



resource "openstack_compute_secgroup_v2" "secgroup_webap" {
  name        = var.secgroup_name
  description = "security group"

  rule {
    from_port   = 22
    to_port     = 22
    ip_protocol = "tcp"
    cidr        = "0.0.0.0/0"
  }
  rule {
    from_port   = 443
    to_port     = 443
    ip_protocol = "tcp"
    cidr        = "0.0.0.0/0"
  }

  rule {
    from_port   = 80
    to_port     = 80
    ip_protocol = "tcp"
    cidr        = "0.0.0.0/0"
  }
  rule {
    from_port   = 8080
    to_port     = 8080
    ip_protocol = "tcp"
    cidr        = "0.0.0.0/0"
  }
  rule {
    from_port   = 8888
    to_port     = 8888
    ip_protocol = "tcp"
    cidr        = "0.0.0.0/0"
  }
}

resource "openstack_networking_port_v2" "port_webap" {
  name               = var.port_name
  network_id         = openstack_networking_network_v2.net_webap.id
  admin_state_up     = "true"
  security_group_ids = [openstack_compute_secgroup_v2.secgroup_webap.id]

  fixed_ip {
    subnet_id  = openstack_networking_subnet_v2.subnet_webap.id
    ip_address = "10.0.0.10"
  }
}

resource "openstack_networking_floatingip_v2" "ip_webap" {
  pool = var.floating_ip_pool
}

resource "openstack_blockstorage_volume_v2" "vol_webap" {
  name = "vol_webap"
  size = var.size_volume
}




resource "openstack_compute_instance_v2" "instance_webap" {
  depends_on = [openstack_networking_subnet_v2.subnet_webap]
  name            = var.instance_name
  image_name      = var.image_name
  flavor_name     = var.flavor_name
  key_pair        = openstack_compute_keypair_v2.my-cloud-key.id
  security_groups = [openstack_compute_secgroup_v2.secgroup_webap.id ]

  network {
    name = var.internal_network_name
  }

}

resource "openstack_compute_volume_attach_v2" "volume_attachement_1" {
  instance_id = openstack_compute_instance_v2.instance_webap.id
  volume_id   = openstack_blockstorage_volume_v2.vol_webap.id
}


resource "openstack_compute_floatingip_associate_v2" "ip_webap" {
  floating_ip = openstack_networking_floatingip_v2.ip_webap.address
  instance_id = openstack_compute_instance_v2.instance_webap.id
  fixed_ip    = openstack_compute_instance_v2.instance_webap.network[0].fixed_ip_v4
}

resource "null_resource" "provision" {
  depends_on = [openstack_compute_floatingip_associate_v2.ip_webap]
  connection {
    type = "ssh"
    user = var.ssh_user
    private_key = file(var.private_ssh_key)
    agent  = "true"
    timeout  = "5m"
    host = openstack_networking_floatingip_v2.ip_webap.address
  }

  provisioner "local-exec" {
    command  =  "bash src/get-sshkeys.sh"
  }
  provisioner "file" {
    source      = "tmp/annotation-cluster/ansible-ubuntu-18.04/roles/common/files/authorized-keys"
    destination = "/home/ubuntu"
  }
  provisioner "file" {
    source      = "src"
    destination = "/home/ubuntu"
  }
  provisioner "remote-exec" {
    inline = [
      "bash /home/ubuntu/src/add-sshkeys.sh",
      "bash /home/ubuntu/src/mount_volume.sh",
      "bash /home/ubuntu/src/install_docker.sh",
      "bash /home/ubuntu/src/run_docker.sh --admin ${var.admin_username} --password ${var.admin_password} "
    ]
  }
}


