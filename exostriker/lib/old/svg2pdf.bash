#!/bin/bash

for i in $@; do
  inkscape --export-area-drawing --without-gui --export-pdf="$(basename $i .svg).pdf" $i
done
