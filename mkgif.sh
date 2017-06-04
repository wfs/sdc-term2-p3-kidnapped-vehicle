#!/bin/bash

anGIF=`basename "$0"`

if [[ $# -ne 6 ]]; then
    echo ""
    echo "$anGIF requires six inputs"
    echo ""
    echo "Example: $anGIF 30 640 0:03 5 input.mp4 output.gif"
    echo ""
    echo "         30         => Frame rate"
    echo "         640        => Target width in pixels"
    echo "         0:03       => Start time in min:sec"
    echo "         5          => Duration from start in sec"
    echo "         input.mp4  => Source movie file"
    echo "         output.gif => Destination animated gif"
    echo ""
    exit 0
fi

palette="/tmp/palette.png"
filters="fps=$1,scale=$2:-1:flags=lanczos"

ffmpeg -v warning -ss $3 -t $4 -i $5 -vf "$filters,palettegen" -y $palette
ffmpeg -v warning -ss $3 -t $4 -i $5 -i $palette -lavfi "$filters [x]; [x][1:v] paletteuse" -y $6
