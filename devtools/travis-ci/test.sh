#!/bin/bash
# This installs the program and runs unit tests
set -e
conda config --add channels omnia
conda build devtools/conda-recipe
conda install --yes --use-local ensembler-dev
conda install --yes nose nose-timer
pushd .; cd /
nosetests ensembler -vv --nocapture --exe -a unit --with-timer

if [[ "$TRAVIS_PULL_REQUEST" == "false" ]]; then
    nosetests ensembler -vv --nocapture --exe -a modeller --with-timer
else
    echo "This is a pull request. Secure environment variables are not available, so will not attempt to run Modeller tests."
fi

popd
