#!/bin/bash

cd AutomatiqueAvancee

jupyter nbconvert --to latex *_temp.ipynb
rm *_temp.ipynb
