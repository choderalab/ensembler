# This installs the program and runs unit tests
conda config --add channels https://conda.anaconda.org/omnia
conda build devtools/conda-recipe
conda install --yes --use-local ensembler-dev
conda install --yes nose
pushd .; cd /
nosetests ensembler -v --exe -a unit

if [[ "$TRAVIS_PULL_REQUEST" == "false" ]]; then
    nosetests ensembler -v --exe -a modeller
else
    echo "This is a pull request. Secure environment variables are not available, so will not attempt to run Modeller tests."
fi

popd