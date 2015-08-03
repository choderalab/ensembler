# This runs unit tests
conda build devtools/conda-recipe

if [[ "$TRAVIS_PULL_REQUEST" == "false" ]]
then
    conda install --yes --use-local ensembler
    conda install --yes nose
    cd /
    nosetests ensembler -v --exe -a modeller
    cd
else
    echo "This is a pull request. Secure environment variables are not available, so will not attempt to run Modeller tests."
fi