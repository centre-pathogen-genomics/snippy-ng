#!/bin/bash

# Install this wrapper as run_clair3.sh in your $PATH
# e.g. install -m 755 docs/run_clair3_docker.sh ~/.local/bin/run_clair3.sh

# Example usage with snippy-ng:
# You must specify the model path inside the container
# e.g. 
# snippy-ng long --ref tests/data/JKD6159.fasta \
#  --reads tests/data/JKD6159.fastq.gz \
#  --clair3-model /opt/models/r1041_e82_400bps_sup_v430_bacteria_finetuned

# Run Clair3 inside a Docker container
docker run -it \
  -v $(pwd):$(pwd) \
  hkubal/clair3:latest \
  /opt/bin/run_clair3.sh $@