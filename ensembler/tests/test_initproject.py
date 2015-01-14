import os

from lxml import etree
from nose.plugins.attrib import attr

import ensembler
import ensembler.cli
import ensembler.cli_commands
import ensembler.initproject
import ensembler.tests
import ensembler.core
import ensembler.UniProt
from ensembler.utils import enter_temp_dir

from ensembler.tests.integration_test_utils import integration_test_context


@attr('unit')
def test_initproject():
    with enter_temp_dir():
        ensembler.initproject.InitProject('.')
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


@attr('unit')
def test_gen_init_metadata():
    initproject_obj = ensembler.initproject.InitProject('.', run_main=False)
    metadata = initproject_obj._gen_init_metadata()
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


@attr('unit')
def test_extract_targets_from_targetexplorer_json():
    targetexplorer_json = {
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
    gather_targets_obj = ensembler.initproject.GatherTargetsFromTargetExplorer('', run_main=False)
    targets = gather_targets_obj._extract_targets_from_json(targetexplorer_json)
    assert (targets[0].id, str(targets[0].seq)) == ('AFC1_ARATH_D0',
                          'YQILSKMGEGTFGQVLECFDNKNKEVVAIKVIRSINKYREAAMIEIDVLQRLTRHDVGGSRCVQIRNWFDYRNHICIVFEKLGPSLYDFLRKNSYRSFPIDLVRELGRQLLESVAYMHDLRLIHTDLKPENILLVSSEYIKIPDYKFLSRPTKDGSYFKNLPKSSAIKLIDFGSTTFEHQDHNYIVSTRHYRAPEVILGVGWNYPCDLWSIGCILVELCSGEALFQTHENLEHLAMMERVLGPLPPHMVLRADRRSEKYFRRGAKLDWPEGATSRDSLKAVWKLPRLPNLIMQHVDHSAGDLIDLLQGLLRYDPTERFKAREALNHPFF')
    assert (targets[1].id, str(targets[1].seq)) == ('AFC2_ARATH_D0',
                          'YKIYSKMGEGTFGQVLECWDRERKEMVAVKIVRGVKKYREAAMIEIEMLQQLGKHDKGGNRCVQIRNWFDYRNHICIVFEKLGSSLYDFLRKNNYRSFPIDLVREIGWQLLECVAFMHDLRMIHTDLKPENILLVSSDYVKIPEYKGSRLQRDVCYKRVPKSSAIKVIDFGSTTYERQDQTYIVSTRHYRAPEVILGLGWSYPCDVWSVGCIIVELCTGEALFQTHENLEHLAMMERVLGPFPQQMLKKVDRHSEKYVRRGRLDWPDGATSRDSLKAVLKLPRLQNLIMQHVDHSAGELINMVQGLLRFDPSERITAREALRHPFF')


@attr('unit')
def test_attempt_symlink_structure_files():
    pdbid = '4CFE'
    structure_paths = [os.path.abspath(os.path.join('tests', 'resources'))]
    with enter_temp_dir():
        os.mkdir('pdb')
        project_pdb_filepath = os.path.join('pdb', pdbid + '.pdb.gz')
        structure_type = 'pdb'
        ensembler.initproject.attempt_symlink_structure_files(pdbid, '.', structure_paths, structure_type)
        assert os.path.exists(project_pdb_filepath)


@attr('unit')
def test_log_unique_domain_names():
    with open(os.path.join('tests', 'resources', 'uniprot-CK1-kinases.xml')) as uniprotxml_file:
        uniprotxml_string = ensembler.UniProt.remove_uniprot_xmlns(uniprotxml_file.read())
        uniprotxml = etree.fromstring(uniprotxml_string)
    ensembler.initproject.log_unique_domain_names('domain:"Protein kinase" AND reviewed:yes AND family:"CK1" AND taxonomy:9606', uniprotxml)


@attr('integration')
def test_command_gather_targets_from_uniprot():
    with integration_test_context(set_up_project_stage='init'):
        ref_fasta = """\
>EGFR_HUMAN_D0
FKKIKVLGSGAFGTVYKGLWIPEGEKVKIPVAIKELREATSPKANKEILDEAYVMASVDN
PHVCRLLGICLTSTVQLITQLMPFGCLLDYVREHKDNIGSQYLLNWCVQIAKGMNYLEDR
RLVHRDLAARNVLVKTPQHVKITDFGLAKLLGAEEKEYHAEGGKVPIKWMALESILHRIY
THQSDVWSYGVTVWELMTFGSKPYDGIPASEISSILEKGERLPQPPICTIDVYMIMVKCW
MIDADSRPKFRELIIEFSKMARDPQRYL
>KC1D_HUMAN_D0
YRLGRKIGSGSFGDIYLGTDIAAGEEVAIKLECVKTKHPQLHIESKIYKMMQGGVGIPTI
RWCGAEGDYNVMVMELLGPSLEDLFNFCSRKFSLKTVLLLADQMISRIEYIHSKNFIHRD
VKPDNFLMGLGKKGNLVYIIDFGLAKKYRDARTHQHIPYRENKNLTGTARYASINTHLGI
EQSRRDDLESLGYVLMYFNLGSLPWQGLKAATKRQKYERISEKKMSTPIEVLCKGYPSEF
ATYLNFCRSLRFDDKPDYSYLRQLFRNLF
"""

        args = {
            '--gather_from': 'uniprot',
            '--query': 'mnemonic:EGFR_HUMAN OR mnemonic:KC1D_HUMAN',
            '--uniprot_domain_regex': '^Protein kinase',
            '--help': False,
        }
        ensembler.cli_commands.gather_targets.dispatch(args)
        test_fasta = open(os.path.join(ensembler.core.default_project_dirnames.targets, 'targets.fa')).read()
        assert test_fasta == ref_fasta


@attr('integration')
def test_gather_templates_from_pdb():
    ref_templates_resolved_seq = """\
>KC1D_HUMAN_D0_4KB8_A
YRLGRKIGSGSFGDIYLGTDIAAGEEVAIKLECVKTKHPQLHIESKIYKMMQGGVGIPTI
RWCGAEGDYNVMVMELLGPSLEDLFNFCSRKFSLKTVLLLADQMISRIEYIHSKNFIHRD
VKPDNFLMGLGKKGNLVYIIDFGLAKKYGTARYASINTHLGIEQSRRDDLESLGYVLMYF
NLGSLPWQGLKERISEKKMSTPIEVLCKGYPSEFATYLNFCRSLRFDDKPDYSYLRQLFR
NLF
>KC1D_HUMAN_D0_4KB8_D
YRLGRKIGDIYLGTDIAAGEEVAIKLECPQLHIESKIYKMMQGGVGIPTIRWCGAEGDYN
VMVMELLGPSLEDLFNFCSRKFSLKTVLLLADQMISRIEYIHSKNFIHRDVKPDNFLMGL
GKKGNLVYIIDFGLAKKYRDAQHIPYRENKNLTGTARYASINTHLGIEQSRRDDLESLGY
VLMYFNLGSLPWQGLKAATKRQKYERISEKKMSTPIEVLCKGYPSEFATYLNFCRSLRFD
DKPDYSYLRQLFRNLF
"""
    with integration_test_context(set_up_project_stage='targets'):
        pdbids = ['4KB8']
        chainids = {'4KB8': ['A', 'D']}
        uniprot_domain_regex = '^Protein kinase'
        ensembler.initproject.gather_templates_from_pdb(pdbids, uniprot_domain_regex, chainids=chainids)
        assert open(os.path.join(ensembler.core.default_project_dirnames.templates, 'templates-resolved-seq.fa')).read() == ref_templates_resolved_seq