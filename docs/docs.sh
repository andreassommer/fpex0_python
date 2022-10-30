#!/bin/bash

# get own directory (i.e. the directory, this script resides in)
OWNDIR=$(dirname -- "${BASH_SOURCE[0]}")

# display help if called without arguments
if [[ $1 = --help ]] || [[ $1 = help ]]; then
   echo ""
   echo "Usage:  docs"
   echo "Builds documentation at docs/ from the project root."
   echo "Depends on pdoc."
   echo ""
   exit 0
fi

# run pdoc
cd $OWNDIR
pdoc ../fpex0 -o ./

# command not found?
out=$?
if [[ $out == "127" ]]; then
    echo ""
    echo "Could not find pdoc."
    echo "Please install, e.g. via 'pip install pdoc'."
elif [[ $out == "0" ]]; then
    echo "Project documentation built at $OWNDIR."
fi