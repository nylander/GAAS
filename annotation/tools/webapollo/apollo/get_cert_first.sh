#!/bin/bash

certbot certonly --webroot -d webapollo.nbis.se  --config-dir ~/certbot/config --logs-dir ~/certbot/logs --work-dir ~/certbot/work

