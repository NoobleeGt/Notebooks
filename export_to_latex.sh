#!/bin/bash

cd Simulation

jupyter nbconvert --to latex *_temp.ipynb
rm *_temp.ipynb
