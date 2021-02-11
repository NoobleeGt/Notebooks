#!/bin/bash

function convert {
  local dires=()

  for elem in $(ls)
  do
    if [ -d $elem ]
    then
      dires+=($elem)
    else
      if [[ $elem == *.svg* ]];
      then
        name=${elem%".svg"}
        inkscape -z -e "$name.png" -w 1024 $elem
      fi
    fi
  done

  for dir in ${dires[@]}
  do
    cd $dir
    convert $dir
  done
  cd ../
}

convert 
