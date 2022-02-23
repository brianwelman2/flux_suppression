#!/bin/bash

cd /home/brian/Degrees/masters/code/python/projects/flux_suppression

MODE=$1
ANT=$2
TINA_PATH=/net/tina/vault2-tina/welman/flux_suppression/images
 
tigger \
$TINA_PATH/true-img/$ANT/src-1/$ANT-true-src-1-$MODE.fits \
$TINA_PATH/true-img/$ANT/src-2/$ANT-true-src-2-$MODE.fits \
$TINA_PATH/true-img/$ANT/src-100/$ANT-true-src-100-$MODE.fits