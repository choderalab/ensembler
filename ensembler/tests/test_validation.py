import os
import yaml
from ensembler.core import default_project_dirnames
from ensembler.validation import molprobity_validation, parse_molprobity_oneline_analysis_output
from ensembler.tests.integrationtest_utils import integrationtest_context
from nose.plugins.attrib import attr

sample_molprobity_oneline_analysis_output_text = """#INPUT  : /dirpath
#pdbFileName:x-H_type:chains:residues:nucacids:resolution:rvalue:rfree:clashscore:clashscoreB<40:minresol:maxresol:n_samples:pct_rank:pct_rank40:cbeta>0.25:numCbeta:rota<1%:numRota:ramaOutlier:ramaAllowed:ramaFavored:numRama:numbadbonds:numbonds:pct_badbonds:pct_resbadbonds:numbadangles:numangles:pct_badangles:pct_resbadangles:MolProbityScore:Mol_pct_rank
KC1D_HUMAN_D0_4KB8_D.pdb:nuclear:1:268:0::::115:0:0:9999:1784:0:100:38:251:13:235:16:38:212:266:133:2197:6.05:40.30:250:2973:8.41:57.46:3.828:4
"""

sample_molprobity_oneline_analysis_output_parsed = {'KC1D_HUMAN_D0_4KB8_D': {'numRota': 235, 'pct_badangles': 8.41, 'chains': 1, 'numangles': 2973, 'pct_rank40': 100, 'numRama': 266, 'nucacids': 0, 'pct_resbadbonds': 40.3, 'numbonds': 2197, 'pct_badbonds': 6.05, 'n_samples': 1784, 'ramaOutlier': 16, 'maxresol': 9999.0, '#pdbFileName': 'KC1D_HUMAN_D0_4KB8_D.pdb', 'x-H_type': 'nuclear', 'Mol_pct_rank': 4, 'MolProbityScore': 3.828, 'rvalue': '', 'residues': 268, 'clashscoreB<40': 0.0, 'clashscore': 115.0, 'pct_resbadangles': 57.46, 'rota<1%': 13, 'numbadangles': 250, 'ramaAllowed': 38, 'rfree': '', 'pct_rank': 0, 'minresol': 0.0, 'numCbeta': 251, 'ramaFavored': 212, 'numbadbonds': 133, 'resolution': '', 'cbeta>0.25': 38}}

ref_pct_badbonds = {
    'KC1D_HUMAN_D0_4KB8_D': 6.05,
    'KC1D_HUMAN_D0_4HNF_A': 6.55,
}


@attr('unit')
def test_parse_molprobity_oneline_analysis_output():
    assert parse_molprobity_oneline_analysis_output(sample_molprobity_oneline_analysis_output_text) == sample_molprobity_oneline_analysis_output_parsed


@attr('non_conda_dependencies')
def test_molprobity_validation():
    target_id = 'EGFR_HUMAN_D0'
    template_ids = ['KC1D_HUMAN_D0_4KB8_D', 'KC1D_HUMAN_D0_4HNF_A']
    with integrationtest_context('refined_explicit'):
        molprobity_validation(target_id)
        for template_id in template_ids:
            results_filepath = os.path.join(
                default_project_dirnames.models, target_id, template_id, 'molprobity-refine_explicit_md.yaml'
            )
            assert os.path.exists(results_filepath)
            with open(results_filepath) as results_file:
                results_dict = yaml.load(results_file)
                assert results_dict.get('pct_badbonds') == ref_pct_badbonds[template_id]

        target_results_filepath = os.path.join(
            default_project_dirnames.models, target_id, 'validation_scores_sorted-molprobity-refine_explicit_md'
        )
        with open(target_results_filepath) as target_results_file:
            target_results = target_results_file.read().splitlines()
            assert target_results[0] == 'KC1D_HUMAN_D0_4KB8_D 3.828'
            assert target_results[1] == 'KC1D_HUMAN_D0_4HNF_A 4.035'
