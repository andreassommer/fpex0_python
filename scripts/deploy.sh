#!/bin/bash


# display help if called without arguments
if [[ $1 = --help ]] || [[ $1 = help ]]; then
   echo "Usage:  deploy [gitremote]"
   echo "   [gitremote] --> git remote to beused (default is used if not specified)"
   echo "Tag is extracted from ../pyproject.toml"
   exit 0
fi

# extract version from pyproject.toml
version=v$( grep "version" ../pyproject.toml | grep -o '".*"' | sed 's/"//g' )


# first argument (if available): remote
if [[ "$1" = "" ]]; then
   gitremote=$( git remote | head -n 1 )
else
   gitremote=$1
fi


# set working directory
parent_path=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )
cd "$parent_path"

# get HEAD's git commit id
commit_hash=$(git rev-parse HEAD)

# commands to be executed
gitUPDATEREF=(git update-ref refs/tags/$version $commit_hash)
gitPUSH=(git push $gitremote $version)

# display commands
echo "Invoking following commands:"
echo ${gitUPDATEREF[@]}
echo ${gitPUSH[@]}
echo "Press Ctrl-C to abort, enter to continue."
read -n 1 -s

# push invokes Package ci-cd workflow on github, jobs build and deploy
# run the commands
"${gitUPDATEREF[@]}"
"${gitPUSH[@]}"


