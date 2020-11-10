#!/bin/bash

here=$(realpath "$(dirname "$0")")
mkdir -p ~/config/serve
rsync -L -au "$here/web/" ~/config/serve/
docker run --rm -p 80:8080 -p 443:8443 -v ~/config:/config:ro -v "$here/nginx.conf":/etc/nginx/nginx.conf:ro --name nginx -d nginx
