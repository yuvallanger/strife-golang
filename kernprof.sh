#!/bin/bash

rm -vf kernprof.h5
kernprof.py -l run.py -d kernprof.h5
