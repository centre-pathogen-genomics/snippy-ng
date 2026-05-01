#!/bin/bash

set -eo pipefail

# Install this wrapper as run_clair3.sh in your $PATH
# e.g. install -m 755 docs/installation/run_clair3_docker.sh ~/.local/bin/run_clair3.sh
#
# The wrapper mounts the parent directories of common Clair3 path arguments so
# absolute host paths continue to work inside the container.

IMAGE="${CLAIR3_DOCKER_IMAGE:-hkubal/clair3:latest}"
PLATFORM_ARG=()
if [[ -n "${CLAIR3_DOCKER_PLATFORM:-}" ]]; then
  PLATFORM_ARG=(--platform "${CLAIR3_DOCKER_PLATFORM}")
fi

declare -a docker_mounts=()
declare -a seen_mounts=()

canonicalize_path() {
  python3 -c 'import os, sys; print(os.path.abspath(os.path.expanduser(sys.argv[1])))' "$1"
}

add_mount() {
  local candidate="$1"
  [[ -z "$candidate" ]] && return 0

  local target="$candidate"
  if [[ "$target" != /* ]]; then
    target="$(canonicalize_path "$PWD/$target")"
  else
    target="$(canonicalize_path "$target")"
  fi

  local mount_path="$target"
  if [[ ! -d "$mount_path" ]]; then
    mount_path="$(dirname "$mount_path")"
  fi

  local seen
  for seen in "${seen_mounts[@]}"; do
    [[ "$seen" == "$mount_path" ]] && return 0
  done

  seen_mounts+=("$mount_path")
  docker_mounts+=(-v "$mount_path:$mount_path")
}

extract_arg_value() {
  local key="$1"
  shift
  local prev=""
  local arg
  for arg in "$@"; do
    if [[ "$prev" == "$key" ]]; then
      printf '%s\n' "$arg"
      return 0
    fi
    if [[ "$arg" == "$key="* ]]; then
      printf '%s\n' "${arg#*=}"
      return 0
    fi
    prev="$arg"
  done
  return 1
}

for key in --model_path --bam_fn --ref_fn --bed_fn --vcf_fn --output; do
  if value="$(extract_arg_value "$key" "$@")"; then
    add_mount "$value"
  fi
done

# Keep current working directory mounted for relative paths and log files.
add_mount "$PWD"

docker run -it --rm \
  "${PLATFORM_ARG[@]}" \
  "${docker_mounts[@]}" \
  "$IMAGE" \
  /opt/bin/run_clair3.sh "$@"
