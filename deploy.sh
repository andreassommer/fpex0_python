#!/bin/bash


# get own directory (i.e. the directory, this script resides in)
OWNDIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )


# parse options, create a test versiontag as default
prefix=tv
while [[ $# ]]; do
   case $1 in
      -t | --test)
      prefix=tv
      shift
      ;;
      
      -r | --release)
      prefix=rv
      shift
      ;;

      *)
      # first non-option argument
      break
   esac
done

# test if pyproject.toml exists
PYPROJECTFILEcandidate=$OWNDIR/pyproject.toml
if [[ ! -f "$PYPROJECTFILEcandidate" ]]; then
  PYPROJECTFILE=NOT_FOUND
  versiontag=NOT_FOUND
else
  # extract version from pyproject.toml
  PYPROJECTFILE=$PYPROJECTFILEcandidate
  versiontag=${prefix}$( grep "version" $PYPROJECTFILE | grep -o '".*"' | sed 's/"//g' )
fi

# display info
echo ""
echo "Project file: $PYPROJECTFILE"
echo "Version tag : $versiontag"
echo ""

# display help if called without arguments
if [[ $1 = --help ]] || [[ $1 = help ]]; then
   echo ""
   echo "Usage:  deploy [option] [gitremote]"
   echo "   [option]    --> One of the following two options:"
   echo "                    * -r or --release for release tag (rv*)"
   echo "                    * -t or --test for test-release tag (tv*)"
   echo "   [gitremote] --> git remote to be used (default is used if not specified)"
   echo ""
   echo "NOTES:"
   echo "  * If no option is given, a test-release tag is created."
   echo "  * You can delete the tag from github using the following command:"
   echo "       git push --delete gitremotename $versiontag "
   echo "  * Local repository tag can be deleted by invoking:"
   echo "       git tag -d $versiontag"
   echo ""
   exit 0
fi

# if no version tag is found, error
if [[ "$versiontag" = "NOT_FOUND" ]]; then
   echo "Could not find project file."
   echo "Expected location: $PYPROJECTFILEcandidate"
   exit 1
fi

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
gitUPDATEREF=(git update-ref refs/tags/$versiontag $commit_hash)
gitPUSH=(git push $gitremote $versiontag)

# display commands
echo ""
echo "Invoking following commands:"
echo "  ${gitUPDATEREF[@]}"
echo "  ${gitPUSH[@]}"
echo ""
echo "Press Ctrl-C to abort, enter to continue."
read -n 1 -s

# push invokes Package ci-cd workflow on github, jobs build and deploy
# run the commands
"${gitUPDATEREF[@]}"
"${gitPUSH[@]}"
