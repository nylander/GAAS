resource "openstack_compute_keypair_v2" "my-cloud-key" {
  name       = var.keypair_name
  public_key =  "${file(var.public_ssh_key)}"

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

resource "openstack_networking_floatingip_v2" "myip" {
  pool = var.floating_ip_pool
}

resource "openstack_compute_instance_v2" "nj_webap" {
  name            = var.instance_name
  image_name      = var.image_name
  flavor_name     = var.flavor_name
  key_pair        = openstack_compute_keypair_v2.my-cloud-key.name
  security_groups = ["default", var.secgroup_name]

  network {
    name = var.internal_network_name
  }

}


resource "openstack_compute_floatingip_associate_v2" "myip" {
  floating_ip = openstack_networking_floatingip_v2.myip.address
  instance_id = openstack_compute_instance_v2.nj_webap.id
  fixed_ip    = openstack_compute_instance_v2.nj_webap.network[0].fixed_ip_v4
}

resource "null_resource" "provision" {
  depends_on = ["openstack_compute_floatingip_associate_v2.myip"]
  connection {
    type = "ssh"
    user = var.ssh_user
    private_key = "${file(var.private_ssh_key)}"
    agent  = true
    timeout  = "5m"
    host = "${openstack_networking_floatingip_v2.myip.address}"
  }

  provisioner "remote-exec" {
    scripts = [
      "src/install_docker.sh",
      "src/run_docker.sh"
    ]
  }
}


