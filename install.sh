#!/usr/bin/env bash
set -euo pipefail

# Download and run the main installer
curl -sSL https://github.com/centre-pathogen-genomics/snippy-ng/releases/latest/download/install.sh | bash -s -- "$@"