# This runs unit tests

if [[ "$TRAVIS_PULL_REQUEST" == "false" ]]
then
    conda install --yes --use-local ensembler
    conda install --yes nose
    pushd .; cd /
    nosetests ensembler -v --exe -a modeller
    popd
else
    echo "This is a pull request. Secure environment variables are not available, so will not attempt to run Modeller tests."
fi
