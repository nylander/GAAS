#!/bin/bash

certbot renew --webroot  --config-dir ~/certbot/config --logs-dir ~/certbot/logs --work-dir ~/certbot/work

rsync -L ~/certbot/config/live/webapollo.nbis.se/privkey.pem ~/config/key.pem
rsync -L ~/certbot/config/live/webapollo.nbis.se/cert.pem ~/config/cert.pem
rsync -L ~/certbot/config/live/webapollo.nbis.se/fullchain.pem ~/config/ca.pem

