terraform {
  required_providers {
    openstack = {
      source = "terraform-provider-openstack/openstack"
    }
    docker = {
      source = "terraform-providers/docker"
    }
  }
  required_version = ">= 0.13"
}
