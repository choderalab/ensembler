def test_extract_residues_from_4CFE():
    import os, StringIO
    import MSMSeeder.PDB

    # 4CFE contains a 'TPO' residue
    pdb_input_filepath = os.path.join('tests', 'resources', '4CFE.pdb')
    desired_chainID = 'A'
    desired_resnums = [str(x) for x in range(16, 269)]
    ofile = StringIO.StringIO()

    nlines_extracted = MSMSeeder.PDB.extract_residues_by_resnum(ofile, pdb_input_filepath, desired_resnums, desired_chainID)
    print len(desired_resnums)
    print nlines_extracted
    ofile.close()

    assert nlines_extracted == len(desired_resnums)

def test_extract_residues_from_3HLL():
    import os, StringIO
    import MSMSeeder.PDB

    # 3HLL contains resnums '56A' and '93B'
    pdb_input_filepath = os.path.join('tests', 'resources', '3HLL.pdb')
    desired_chainID = 'A'
    desired_resnums = [str(x) for x in range(24, 172) + range(183, 309)]
    desired_resnums[desired_resnums.index('56')] = '56A'
    desired_resnums[desired_resnums.index('93')] = '93B'
    ofile = StringIO.StringIO()

    nlines_extracted = MSMSeeder.PDB.extract_residues_by_resnum(ofile, pdb_input_filepath, desired_resnums, desired_chainID)
    print len(desired_resnums)
    print nlines_extracted
    print ofile.getvalue()
    ofile.close()

    assert nlines_extracted == len(desired_resnums)

