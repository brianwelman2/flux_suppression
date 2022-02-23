#!/usr/bin/bash

cd /home/brian/Degrees/masters/code/python/projects/flux_suppression/

# Path to folder in tina
TINA_PATH=welman@tina.ru.ac.za:/home/welman/masters/projects/flux_suppression

# Transfer folders and files
scp -r $TINA_PATH/plots/ ~/Degrees/masters/code/python/projects/flux_suppression/plots/