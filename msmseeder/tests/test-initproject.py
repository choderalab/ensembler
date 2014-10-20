import os
import msmseeder.initproject
import msmseeder.tests
import msmseeder.core


def test_initproject():
    with msmseeder.tests.utils.enter_temp_directory():
        msmseeder.initproject.initproject('.')
        for dirpath in msmseeder.core.project_dirnames:
            assert os.path.exists(dirpath)
        assert os.path.exists('meta.yaml')


def test_gen_init_metadata():
    metadata = msmseeder.initproject.gen_init_metadata('.')
    init_metadata = metadata.get('init')
    assert init_metadata is not None
    metadata_keys = [
        'datestamp',
        'init_path',
        'python_version',
        'python_full_version',
        'msmseeder_version',
        'msmseeder_commit'
    ]
    for key in metadata_keys:
        assert init_metadata.get(key) is not None


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
    assert targets[0] == ('AFC1_ARATH_D0', 'YQILSKMGEGTFGQVLECFDNKNKEVVAIKVIRSINKYREAAMIEIDVLQRLTRHDVGGS\nRCVQIRNWFDYRNHICIVFEKLGPSLYDFLRKNSYRSFPIDLVRELGRQLLESVAYMHDL\nRLIHTDLKPENILLVSSEYIKIPDYKFLSRPTKDGSYFKNLPKSSAIKLIDFGSTTFEHQ\nDHNYIVSTRHYRAPEVILGVGWNYPCDLWSIGCILVELCSGEALFQTHENLEHLAMMERV\nLGPLPPHMVLRADRRSEKYFRRGAKLDWPEGATSRDSLKAVWKLPRLPNLIMQHVDHSAG\nDLIDLLQGLLRYDPTERFKAREALNHPFF')
    assert targets[1] == ('AFC2_ARATH_D0', 'YKIYSKMGEGTFGQVLECWDRERKEMVAVKIVRGVKKYREAAMIEIEMLQQLGKHDKGGN\nRCVQIRNWFDYRNHICIVFEKLGSSLYDFLRKNNYRSFPIDLVREIGWQLLECVAFMHDL\nRMIHTDLKPENILLVSSDYVKIPEYKGSRLQRDVCYKRVPKSSAIKVIDFGSTTYERQDQ\nTYIVSTRHYRAPEVILGLGWSYPCDVWSVGCIIVELCTGEALFQTHENLEHLAMMERVLG\nPFPQQMLKKVDRHSEKYVRRGRLDWPDGATSRDSLKAVLKLPRLQNLIMQHVDHSAGELI\nNMVQGLLRFDPSERITAREALRHPFF')


def test_attempt_symlink_structure_files():
    pdbid = '4CFE'
    structure_paths = [os.path.abspath(os.path.join('tests', 'resources'))]
    with msmseeder.tests.utils.enter_temp_directory():
        project_pdb_filepath = os.path.join(pdbid + '.pdb.gz')
        structure_type = 'pdb'
        msmseeder.initproject.attempt_symlink_structure_files(pdbid, '.', structure_paths, structure_type)
        assert os.path.exists(project_pdb_filepath)