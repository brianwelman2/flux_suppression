#!/bin/bash

cd /home/brian/Degrees/masters/code/python/projects/flux_suppression

MODE=$1
TINT=$2

TINA_PATH=/net/tina/vault2-tina/welman/flux_suppression/images

tigger \
$TINA_PATH/true-img/true-$MODE.fits \
$TINA_PATH/noisy-img/noisy-$MODE.fits \
$TINA_PATH/cubical-img/tint_$TINT/cubical-*-$MODE.fits \
$TINA_PATH/kalcal-img/kalcal-*-$MODE.fits