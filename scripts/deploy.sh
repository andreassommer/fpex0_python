#!/bin/bash


# display help if called without arguments
if [[ $1 = --help ]] || [[ $1 = help ]] || [[ $# -eq 0 ]]; then
   echo "Usage:  deploy version"
   echo "   version --> version number"
   exit 0
fi


# first argument: version numberg, preponing v
version=v$1

# second argument (if available): remote
if [[ "$2" = "" ]]
   gitremote=$( git remote | head -n 1 )
elif
   gitremote=$2
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
echo $gitUPDATEREF
echo $gitPUSH
echo "Press Ctrl-C to abort, enter to continue.'
read -n 1 -s

# push invokes Package ci-cd workflow on github, jobs build and deploy
# run the commands
#"${gitUPDATEREF[@]}"
#"${gitPUSH[@]}"


