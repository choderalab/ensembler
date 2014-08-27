def test_build_model():
    import os
    import Bio.SeqIO
    import msmseeder.modelling
    from mock import Mock
    # given
    target = Mock()
    template = Mock()
    target.id = 'mock_target'
    target.seq = 'YILGDTLGVGGKVKVGKH'
    template.id = 'mock_template'
    template.seq = 'YQNLSPVGSGAYGSVCAAFD'

    # when
    output_text = msmseeder.modelling.build_model(target,
                                    template,
                                    template_structure_dir='tests/resources',
                                    aln_filepath='tests/resources/alignment.pir',
                                    seqid_filepath='tests/resources/sequence-identity.txt',
                                    model_pdbfilepath='tests/resources/model.pdb.gz',
                                    restraint_filepath='tests/resources/restraints.rsr.gz')

    # then

    # import gzip
    # with gzip.open('tests/resources/model.pdb.gz') as pdb_file:
    #     for i in 5:
    #         print pdb_file.readline(),

    # Example model.pdb.gz contents (not asserting this as it may be dependent
    # upon Modeller version as well as modelling stochasticity):
    #
    # ..EXPDTA    THEORETICAL MODEL, MODELLER 9.12 2014/08/26 13:15:44
    # REMARK   6 MODELLER OBJECTIVE FUNCTION:       326.6798
    # REMARK   6 MODELLER BEST TEMPLATE % SEQ ID:  27.778
    # ATOM      1  N   TYR     1      48.812  50.583  13.949  1.00110.28           N
    # ATOM      2  CA  TYR     1      49.070  50.334  15.387  1.00110.28           C

    assert os.path.exists('tests/resources/model.pdb.gz')
    assert os.path.getsize('tests/resources/model.pdb.gz') > 0

    os.remove('tests/resources/model.pdb')
    os.remove('tests/resources/alignment.pir')
    os.remove('tests/resources/sequence-identity.txt')
    os.remove('tests/resources/model.pdb.gz')
    os.remove('tests/resources/restraints.rsr.gz')
