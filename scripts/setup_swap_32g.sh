#!/usr/bin/env bash
set -euo pipefail
if swapon --show | grep -q /swap_scfactory_32g.img; then
  echo "swap already active"
  swapon --show
  exit 0
fi
sudo fallocate -l 32G /swap_scfactory_32g.img || sudo dd if=/dev/zero of=/swap_scfactory_32g.img bs=1M count=32768 status=progress
sudo chmod 600 /swap_scfactory_32g.img
sudo mkswap /swap_scfactory_32g.img
sudo swapon /swap_scfactory_32g.img
if ! grep -q "/swap_scfactory_32g.img" /etc/fstab; then
  echo "/swap_scfactory_32g.img none swap sw 0 0" | sudo tee -a /etc/fstab >/dev/null
fi
swapon --show
free -h
