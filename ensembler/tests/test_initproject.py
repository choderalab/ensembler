import os
import ensembler.initproject
import ensembler.tests
import ensembler.core
import ensembler.UniProt
from ensembler.utils import enter_temp_dir
from lxml import etree


def test_initproject():
    with enter_temp_dir():
        ensembler.initproject.initproject('.')
        assert os.path.exists(ensembler.core.default_project_dirnames.targets)
        assert os.path.exists(ensembler.core.default_project_dirnames.templates)
        assert os.path.exists(ensembler.core.default_project_dirnames.structures)
        assert os.path.exists(ensembler.core.default_project_dirnames.models)
        assert os.path.exists(ensembler.core.default_project_dirnames.packaged_models)
        assert os.path.exists(ensembler.core.default_project_dirnames.structures_pdb)
        assert os.path.exists(ensembler.core.default_project_dirnames.structures_sifts)
        assert os.path.exists(ensembler.core.default_project_dirnames.templates_structures_resolved)
        assert os.path.exists(ensembler.core.default_project_dirnames.templates_structures_modeled_loops)
        assert os.path.exists('meta0.yaml')


def test_gen_init_metadata():
    metadata = ensembler.initproject.gen_init_metadata('.')
    metadata_keys = [
        'datestamp',
        'init_path',
        'python_version',
        'python_full_version',
        'ensembler_version',
        'ensembler_commit'
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
    targets = ensembler.initproject.extract_targets_from_targetexplorer_json(targets_json)
    assert targets[0] == ('AFC1_ARATH_D0',
                          'YQILSKMGEGTFGQVLECFDNKNKEVVAIKVIRSINKYREAAMIEIDVLQRLTRHDVGGSRCVQIRNWFDYRNHICIVFEKLGPSLYDFLRKNSYRSFPIDLVRELGRQLLESVAYMHDLRLIHTDLKPENILLVSSEYIKIPDYKFLSRPTKDGSYFKNLPKSSAIKLIDFGSTTFEHQDHNYIVSTRHYRAPEVILGVGWNYPCDLWSIGCILVELCSGEALFQTHENLEHLAMMERVLGPLPPHMVLRADRRSEKYFRRGAKLDWPEGATSRDSLKAVWKLPRLPNLIMQHVDHSAGDLIDLLQGLLRYDPTERFKAREALNHPFF')
    assert targets[1] == ('AFC2_ARATH_D0',
                          'YKIYSKMGEGTFGQVLECWDRERKEMVAVKIVRGVKKYREAAMIEIEMLQQLGKHDKGGNRCVQIRNWFDYRNHICIVFEKLGSSLYDFLRKNNYRSFPIDLVREIGWQLLECVAFMHDLRMIHTDLKPENILLVSSDYVKIPEYKGSRLQRDVCYKRVPKSSAIKVIDFGSTTYERQDQTYIVSTRHYRAPEVILGLGWSYPCDVWSVGCIIVELCTGEALFQTHENLEHLAMMERVLGPFPQQMLKKVDRHSEKYVRRGRLDWPDGATSRDSLKAVLKLPRLQNLIMQHVDHSAGELINMVQGLLRFDPSERITAREALRHPFF')


def test_attempt_symlink_structure_files():
    pdbid = '4CFE'
    structure_paths = [os.path.abspath(os.path.join('tests', 'resources'))]
    with enter_temp_dir():
        os.mkdir('pdb')
        project_pdb_filepath = os.path.join('pdb', pdbid + '.pdb.gz')
        structure_type = 'pdb'
        ensembler.initproject.attempt_symlink_structure_files(pdbid, '.', structure_paths, structure_type)
        assert os.path.exists(project_pdb_filepath)


def test_log_unique_domain_names():
    with open(os.path.join('tests', 'resources', 'uniprot-CK1-kinases.xml')) as uniprotxml_file:
        uniprotxml_string = ensembler.UniProt.remove_uniprot_xmlns(uniprotxml_file.read())
        uniprotxml = etree.fromstring(uniprotxml_string)
    ensembler.initproject.log_unique_domain_names('domain:"Protein kinase" AND reviewed:yes AND family:"CK1" AND taxonomy:9606', uniprotxml)