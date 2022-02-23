#!/bin/bash

cd /home/brian/Degrees/masters/code/python/projects/flux_suppression

MODE=$1
SIGMA_F=$2
ANT=$3

if [ $ANT = kat7 ]
then
    tigger \
    images/kalcal-img/$ANT/src-1/fp-100/sigma_$SIGMA_F/$ANT-fp-100-src-1-$SIGMA_F-$MODE.fits \
    images/kalcal-img/$ANT/src-2/fp-100/sigma_$SIGMA_F/$ANT-fp-100-src-2-$SIGMA_F-$MODE.fits \
    images/kalcal-img/$ANT/src-2/fp-70/sigma_$SIGMA_F/$ANT-fp-70-src-2-$SIGMA_F-$MODE.fits
else
    tigger \
    images/kalcal-img/$ANT/src-1/fp-100/sigma_$SIGMA_F/$ANT-fp-100-src-1-$SIGMA_F-$MODE.fits \
    images/kalcal-img/$ANT/src-2/fp-100/sigma_$SIGMA_F/$ANT-fp-100-src-2-$SIGMA_F-$MODE.fits \
    images/kalcal-img/$ANT/src-2/fp-70/sigma_$SIGMA_F/$ANT-fp-70-src-2-$SIGMA_F-$MODE.fits \
    images/kalcal-img/$ANT/src-100/fp-100/sigma_$SIGMA_F/$ANT-fp-100-src-100-$SIGMA_F-$MODE.fits \
    images/kalcal-img/$ANT/src-100/fp-70/sigma_$SIGMA_F/$ANT-fp-70-src-100-$SIGMA_F-$MODE.fits \
    images/kalcal-img/$ANT/src-100/fp-50/sigma_$SIGMA_F/$ANT-fp-50-src-100-$SIGMA_F-$MODE.fits \
    images/kalcal-img/$ANT/src-100/fp-30/sigma_$SIGMA_F/$ANT-fp-30-src-100-$SIGMA_F-$MODE.fits
fi