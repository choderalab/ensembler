import os
import shutil
import msmseeder.initproject
import msmseeder.tests
import msmseeder.core
from msmseeder.utils import enter_temp_dir
from mock import Mock


def test_initproject():
    with enter_temp_dir():
        msmseeder.initproject.initproject('.')
        assert os.path.exists(msmseeder.core.default_project_dirnames.targets)
        assert os.path.exists(msmseeder.core.default_project_dirnames.templates)
        assert os.path.exists(msmseeder.core.default_project_dirnames.structures)
        assert os.path.exists(msmseeder.core.default_project_dirnames.models)
        assert os.path.exists(msmseeder.core.default_project_dirnames.packaged_models)
        assert os.path.exists(msmseeder.core.default_project_dirnames.structures_pdb)
        assert os.path.exists(msmseeder.core.default_project_dirnames.structures_sifts)
        assert os.path.exists(msmseeder.core.default_project_dirnames.templates_structures_resolved)
        assert os.path.exists(msmseeder.core.default_project_dirnames.templates_structures_modeled_loops)
        assert os.path.exists('meta0.yaml')


def test_gen_init_metadata():
    metadata = msmseeder.initproject.gen_init_metadata('.')
    metadata_keys = [
        'datestamp',
        'init_path',
        'python_version',
        'python_full_version',
        'msmseeder_version',
        'msmseeder_commit'
    ]
    for key in metadata_keys:
        assert metadata.get(key) is not None


def test_extract_targets_from_targetexplorer_json():
    targets_json = {
        "results": [
            {
                "ac": "P51566",
                "domains": [
                    {
                        "sequence": "YQILSKMGEGTFGQVLECFDNKNKEVVAIKVIRSINKYREAAMIEIDVLQRLTRHDVGGSRCVQIRNWFDYRNHICIVFEKLGPSLYDFLRKNSYRSFPIDLVRELGRQLLESVAYMHDLRLIHTDLKPENILLVSSEYIKIPDYKFLSRPTKDGSYFKNLPKSSAIKLIDFGSTTFEHQDHNYIVSTRHYRAPEVILGVGWNYPCDLWSIGCILVELCSGEALFQTHENLEHLAMMERVLGPLPPHMVLRADRRSEKYFRRGAKLDWPEGATSRDSLKAVWKLPRLPNLIMQHVDHSAGDLIDLLQGLLRYDPTERFKAREALNHPFF",
                        "targetid": "AFC1_ARATH_D0"
                    }
                ],
                "entry_name": "AFC1_ARATH",
                "family": "CMGC",
                "nbioassays": 0,
                "npdbs": 0,
                "npubs": 3,
                "sequence": "MQSSVYRDKASSIAMILETQRNVEFPHRIVDKRPRKRPRLTWDAAPPLLPPPPPPTVFQPPLYYGPEFASGLVPNFVYPNMFYNGLPRQGSPPWRPDDKDGHYVFVVGDTLTPRYQILSKMGEGTFGQVLECFDNKNKEVVAIKVIRSINKYREAAMIEIDVLQRLTRHDVGGSRCVQIRNWFDYRNHICIVFEKLGPSLYDFLRKNSYRSFPIDLVRELGRQLLESVAYMHDLRLIHTDLKPENILLVSSEYIKIPDYKFLSRPTKDGSYFKNLPKSSAIKLIDFGSTTFEHQDHNYIVSTRHYRAPEVILGVGWNYPCDLWSIGCILVELCSGEALFQTHENLEHLAMMERVLGPLPPHMVLRADRRSEKYFRRGAKLDWPEGATSRDSLKAVWKLPRLPNLIMQHVDHSAGDLIDLLQGLLRYDPTERFKAREALNHPFFTRSREQSIPPFNPNPHPFLYNQKN"
            },
            {
                "ac": "P51567",
                "domains": [
                    {
                        "sequence": "YKIYSKMGEGTFGQVLECWDRERKEMVAVKIVRGVKKYREAAMIEIEMLQQLGKHDKGGNRCVQIRNWFDYRNHICIVFEKLGSSLYDFLRKNNYRSFPIDLVREIGWQLLECVAFMHDLRMIHTDLKPENILLVSSDYVKIPEYKGSRLQRDVCYKRVPKSSAIKVIDFGSTTYERQDQTYIVSTRHYRAPEVILGLGWSYPCDVWSVGCIIVELCTGEALFQTHENLEHLAMMERVLGPFPQQMLKKVDRHSEKYVRRGRLDWPDGATSRDSLKAVLKLPRLQNLIMQHVDHSAGELINMVQGLLRFDPSERITAREALRHPFF",
                        "targetid": "AFC2_ARATH_D0"
                    }
                ],
                "entry_name": "AFC2_ARATH",
                "family": "CMGC",
                "nbioassays": 0,
                "npdbs": 0,
                "npubs": 2,
                "sequence": "MEMERVHEFPHTHMDRRPRKRARLGWDVLPQATKAQVGMFCGQEIGNISSFASSGAPSDNSSSLCVKGVARNGSPPWREDDKDGHYIFELGDDLTPRYKIYSKMGEGTFGQVLECWDRERKEMVAVKIVRGVKKYREAAMIEIEMLQQLGKHDKGGNRCVQIRNWFDYRNHICIVFEKLGSSLYDFLRKNNYRSFPIDLVREIGWQLLECVAFMHDLRMIHTDLKPENILLVSSDYVKIPEYKGSRLQRDVCYKRVPKSSAIKVIDFGSTTYERQDQTYIVSTRHYRAPEVILGLGWSYPCDVWSVGCIIVELCTGEALFQTHENLEHLAMMERVLGPFPQQMLKKVDRHSEKYVRRGRLDWPDGATSRDSLKAVLKLPRLQNLIMQHVDHSAGELINMVQGLLRFDPSERITAREALRHPFFARRR"
            },
        ]
    }
    targets = msmseeder.initproject.extract_targets_from_targetexplorer_json(targets_json)
    assert targets[0] == ('AFC1_ARATH_D0', 'YQILSKMGEGTFGQVLECFDNKNKEVVAIKVIRSINKYREAAMIEIDVLQRLTRHDVGGSRCVQIRNWFDYRNHICIVFEKLGPSLYDFLRKNSYRSFPIDLVRELGRQLLESVAYMHDLRLIHTDLKPENILLVSSEYIKIPDYKFLSRPTKDGSYFKNLPKSSAIKLIDFGSTTFEHQDHNYIVSTRHYRAPEVILGVGWNYPCDLWSIGCILVELCSGEALFQTHENLEHLAMMERVLGPLPPHMVLRADRRSEKYFRRGAKLDWPEGATSRDSLKAVWKLPRLPNLIMQHVDHSAGDLIDLLQGLLRYDPTERFKAREALNHPFF')
    assert targets[1] == ('AFC2_ARATH_D0', 'YKIYSKMGEGTFGQVLECWDRERKEMVAVKIVRGVKKYREAAMIEIEMLQQLGKHDKGGNRCVQIRNWFDYRNHICIVFEKLGSSLYDFLRKNNYRSFPIDLVREIGWQLLECVAFMHDLRMIHTDLKPENILLVSSDYVKIPEYKGSRLQRDVCYKRVPKSSAIKVIDFGSTTYERQDQTYIVSTRHYRAPEVILGLGWSYPCDVWSVGCIIVELCTGEALFQTHENLEHLAMMERVLGPFPQQMLKKVDRHSEKYVRRGRLDWPDGATSRDSLKAVLKLPRLQNLIMQHVDHSAGELINMVQGLLRFDPSERITAREALRHPFF')


def test_attempt_symlink_structure_files():
    pdbid = '4CFE'
    structure_paths = [os.path.abspath(os.path.join('tests', 'resources'))]
    with enter_temp_dir():
        os.mkdir('pdb')
        project_pdb_filepath = os.path.join('pdb', pdbid + '.pdb.gz')
        structure_type = 'pdb'
        msmseeder.initproject.attempt_symlink_structure_files(pdbid, '.', structure_paths, structure_type)
        assert os.path.exists(project_pdb_filepath)


# def test_pdbfix_template():
#     template_filepath = os.path.abspath(os.path.join('tests', 'resources', 'mock_template.pdb'))
#     with enter_temp_dir():
#         os.makedirs(msmseeder.core.default_project_dirnames.templates_structures_resolved)
#         os.makedirs(msmseeder.core.default_project_dirnames.templates_structures_modeled_loops)
#         shutil.copy(template_filepath, os.path.join(msmseeder.core.default_project_dirnames.templates_structures_resolved, 'mock_template.pdb'))
#         template = Mock()
#         template.chainid = 'A'
#         template.templateid = 'mock_template'
#         template.resolved_seq = 'YQNLSPVGSGGSVCAAFD'
#         template.full_seq = 'YQNLSPVGSGAYGSVCAAFD'
#         msmseeder.initproject.pdbfix_template(template)
#
#         pdbfixed_template_filepath = os.path.join(msmseeder.core.default_project_dirnames.templates_structures_modeled_loops, 'mock_template.pdb')
#         assert os.path.exists(pdbfixed_template_filepath)