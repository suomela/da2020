#!/bin/bash

set -e

for a in token leader all wave-1 wave-2 wave-3 apsp slow; do
    echo "$a"
    rm -f "video/$a.mp4"
    ffmpeg -framerate 30 -pattern_type glob -i "figs/$a-*-*.png" -vf format=yuv420p "video/$a.mp4"
done
