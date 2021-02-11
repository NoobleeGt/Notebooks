#!/bin/bash

cd Electrotechnique

jupyter nbconvert --to latex *_temp.ipynb
rm *_temp.ipynb
