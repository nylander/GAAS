#!/bin/bash

if type ruby >/dev/null 2>/dev/null; then
  :
else
  echo 	"Missing ruby"
fi

if type groovy >/dev/null 2>/dev/null; then
  :
else
  echo  "Missing groovy"
fi

if [ -d "$APOLLO_DATA_DIR" ]; then
  : 
else
  echo "APOLLO_DATA_DIR is not set or does not point to a directory"
fi

