#!/bin/bash

for i in $@; do
  #inkscape --export-area-drawing --without-gui --export-pdf="$(basename $i .svg).pdf" $i
  rsvg-convert -f pdf -o $(basename $i .svg).pdf $(basename $i .svg).svg
done
