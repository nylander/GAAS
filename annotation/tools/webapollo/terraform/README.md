# Deploy an instance on SSC with webapollo using Terraform 

## Usage

1. create an ssh-key by 

    ```
    mkdir -p private
    ssh-keygen -q -t rsa  -N "" -f private/webap.key 
    ```

2. Initialize terraform by 

    `terraform init`

3. Create a plain text file `terraform.tfvars` and define parameters, e.g. 

    ```
    external_network_id="8a2049af-7ff7-4303-b794-a6387dfee03d" # id of the external network
    floating_ip_pool = "Public External IPv4 Network"          # name of the external floatig ip pool


    size_volume = 100 # define size of block storage to be mounted in GB 

    admin_username = "adminuser"
    admin_password = "verystrongpassword"
    ```

4. Set up Openstack credentials by sourcing an RC file from an SSC project

    `source <rc-file.sh>`

5. Apply the deployment by 

    `terraform apply`


After successful deployment, the webapollo server can be accessed at http://instance-ip:8888/

