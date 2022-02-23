#!/usr/bin/bash

cd /home/brian/Degrees/masters/code/python/projects/flux_suppression/

# Path to folder in tina
TINA_PATH=welman@tina.ru.ac.za:/home/welman/masters/projects/flux_suppression

# Transfer folders and files
scp -r configs/ $TINA_PATH/.
scp -r logs/ $TINA_PATH/.
scp -r scripts/ $TINA_PATH/.
scp -r skymodels/ $TINA_PATH/.

# Make directories and copy to tina
find cubical gains images plots -type d > dirs.txt
scp dirs.txt $TINA_PATH/.
rm dirs.txt