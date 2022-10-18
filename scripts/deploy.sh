#!/bin/bash
# - doc -
# takes: 1 appropriate python-package version string
# 
# creates a tag named $1 and pushs it to the repo
# at the repo a workflow is invoked if triggering on pushed tag
# -  -  -

# python3 -m updateVersion
# version=$(cat version)
# --> changed because pyproject.toml is hardly readable

# set working directory
parent_path=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )
cd "$parent_path"

version=$1
# create tag
commit_hash=$(cat ../.git/refs/heads/master)
git update-ref refs/tags/v$version $commit_hash

# push invokes Package ci-cd workflow on github, jobs build and deploy
git push origin v$version