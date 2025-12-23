#!/bin/bash

# Install this wrapper as run_clair3.sh in your $PATH
# e.g. install -m 755 docs/run_clair3_docker.sh ~/.local/bin/run_clair3.sh


# Run Clair3 inside a Docker container
docker run -it \
  -v $(pwd):$(pwd) \
  hkubal/clair3:latest \
  /opt/bin/run_clair3.sh $@