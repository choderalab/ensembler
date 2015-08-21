#!/bin/bash
echo $TRAVIS_PULL_REQUEST $TRAVIS_BRANCH

if [[ "$TRAVIS_PULL_REQUEST" != "false" ]]; then
    echo "This is a pull request. No deployment will be done."; exit 0
fi


if [[ "$TRAVIS_BRANCH" != "master" ]]; then
    echo "No deployment on BRANCH='$TRAVIS_BRANCH'"; exit 0
fi


binstar -t $BINSTAR_TOKEN upload --force -u omnia -p ensembler-dev $HOME/miniconda/conda-bld/linux-64/ensembler-dev-*
