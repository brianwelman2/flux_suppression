#!/bin/bash

cd /net/tina/vault2-tina/welman/flux_suppression/

# Delete all gains solutions created from kalcal
find gains/ -name "sigma_*" -exec rm -r {} \;


# Delete all restored fluxes
find fluxes/ -name "sigma_*" -exec rm -r {} \;
