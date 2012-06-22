#!/bin/sh

run=1
while [ run != 0 ]; do
    ./game_of_strife.py
    run=$?
done
