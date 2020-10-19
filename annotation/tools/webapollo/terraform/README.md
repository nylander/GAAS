# Deploy an instance on SSC with webapollo using Terraform 

## Usage

1. create an ssh-key by 

    ```
    mkdir -p private
    ssh-keygen -q -t rsa  -N "" -f private/webap.key 
    ```

2. Initialize terraform by 

    `terraform init`

3. Set up Openstack credentials by sourcing an RC file from an SSC project

    `source <rc-file.sh>`

4. Apply the deployment by 

    `terraform apply`


After successful deployment, the webapollo server can be accessed at http://instance-ip:8888/

