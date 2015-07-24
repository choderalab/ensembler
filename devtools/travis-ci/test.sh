# This runs unit tests
conda build devtools/conda-recipe

if [[ "$TRAVIS_PULL_REQUEST" == "true" ]]
then
    echo "This is a pull request. Secure environment variables are not available, so will not attempt to run Modeller tests."
else
    conda build devtools/conda-recipe
    conda install --yes --use-local ensembler
    conda install --yes nose
    nosetests ensembler -v --exe -a modeller
fi