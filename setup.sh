#!/usr/bin/env bash

sudo apt-get install ffmpeg gnuplot openmpi-bin -y

if [[ "$?" == "0" ]]; then
    echo "SETUP SUCCESSFUL"
else
    echo "SETUP FAILED"
fi

exit 0
