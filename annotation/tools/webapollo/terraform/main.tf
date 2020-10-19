resource "openstack_compute_keypair_v2" "my-cloud-key" {
  name       = "nj-webap"
  public_key =  "${file("private/webap.key.pub")}"

}

resource "openstack_compute_secgroup_v2" "secgroup_webap" {
  name        = "secgroup_webap"
  description = "Security group for webapollo"

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
  pool = "Public External IPv4 Network"
}

resource "openstack_compute_instance_v2" "nj_webap" {
  name            = "webap-vm-1"
  image_name      = "Ubuntu 18.04 LTS (Bionic Beaver) - latest"
  flavor_name     = "ssc.large"
  key_pair        = openstack_compute_keypair_v2.my-cloud-key.name
  security_groups = ["default", "secgroup_webap"]

  network {
    name = "SNIC 2020/20-30 Internal IPv4 Network"
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
    user = "ubuntu"
    private_key = "${file("private/webap.key")}"
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


