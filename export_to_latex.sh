#!/bin/bash

cd StateSpaceControl

jupyter nbconvert --to latex *_temp.ipynb
rm *_temp.ipynb
