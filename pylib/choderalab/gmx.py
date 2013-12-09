# Library for helping with Gromacs set-up and simulation
#
# Daniel L. Parton <partond@mskcc.org> - 13 Mar 2013
#
import re

# 3-letter amino acid codes
aa3codes = ['ALA','CYS','ASP','GLU','PHE','GLY','HIS','HIE','HID','HIP','ILE','LYS','LEU','MET','ASN','PRO','GLN','ARG','SER','THR','VAL','TRP','TYR']

# data structure containing aa hydrogen names for gmx and openmm versions of amber99sb-ildn
# Used for mapping between the two
gmx_omm_Hname_mapping = {
'gmx':
    {
    'GLY':['HA1','HA2'],
    'SER':['HB1','HB2'],
    'LEU':['HB1','HB2'],
    'ILE':['HG11','HG12'],
    'ASN':['HB1','HB2'],
    'GLN':['HB1','HB2','HG1','HG2'],
    'ARG':['HB1','HB2','HG1','HG2','HD1','HD2'],
    'HID':['HB1','HB2'],
    'HIE':['HB1','HB2'],
    'HIP':['HB1','HB2'],
    'HIS':['HB1','HB2'], # Note that HIS doesn't appear in the amber force field, but OpenMM uses this to represent all histidines (HID,HIE and HIP are listed as variants in the amber*.xml file)
    'TRP':['HB1','HB2'],
    'PHE':['HB1','HB2'],
    'TYR':['HB1','HB2'],
    'ASP':['HB1','HB2'],
    'GLU':['HB1','HB2','HG1','HG2'],
    'LYS':['HB1','HB2','HG1','HG2','HD1','HD2','HE1','HE2'],
    'PRO':['HD1','HD2','HG1','HG2','HB1','HB2'],
    'CYS':['HB1','HB2'],
    'MET':['HB1','HB2','HG1','HG2'],
    'ALA':[],
    'THR':[],
    'VAL':[]
    },
'openmm':
    {
    'GLY':['HA2','HA3'],
    'SER':['HB2','HB3'],
    'LEU':['HB2','HB3'],
    'ILE':['HG12','HG13'],
    'ASN':['HB2','HB3'],
    'GLN':['HB2','HB3','HG2','HG3'],
    'ARG':['HB2','HB3','HG2','HG3','HD2','HD3'],
    'HID':['HB2','HB3'],
    'HIE':['HB2','HB3'],
    'HIP':['HB2','HB3'],
    'HIS':['HB2','HB3'],
    'TRP':['HB2','HB3'],
    'PHE':['HB2','HB3'],
    'TYR':['HB2','HB3'],
    'ASP':['HB2','HB3'],
    'GLU':['HB2','HB3','HG2','HG3'],
    'LYS':['HB2','HB3','HG2','HG3','HD2','HD3','HE2','HE3'],
    'PRO':['HD2','HD3','HG2','HG3','HB2','HB3'],
    'CYS':['HB2','HB3'],
    'MET':['HB2','HB3','HG2','HG3'],
    'ALA':[],
    'THR':[],
    'VAL':[]
    }
}

def Hmap_omm2gmx(resname,omm_atomname):
    '''Converts OpenMM style amber99sb-ildn Hydrogen names to Gromacs style.'''
    if omm_atomname in gmx_omm_Hname_mapping['openmm'][resname]:
        atom_index = gmx_omm_Hname_mapping['openmm'][resname].index(omm_atomname)
        gmx_atomname = gmx_omm_Hname_mapping['gmx'][resname][atom_index]
        return gmx_atomname
    else:
        return omm_atomname

def Hmap_gmx2omm(resname,gmx_atomname):
    '''Converts Gromacs style amber99sb-ildn Hydrogen names to OpenMM style.'''
    if gmx_atomname in gmx_omm_Hname_mapping['gmx'][resname]:
        atom_index = gmx_omm_Hname_mapping['gmx'][resname].index(gmx_atomname)
        omm_atomname = gmx_omm_Hname_mapping['openmm'][resname][atom_index]
        return omm_atomname
    else:
        return gmx_atomname

def Hmap_pdb_omm2gmx(pdb_line):
    '''Converts OpenMM style amber99sb-ildn Hydrogen names in OpenMM .pdb file to Gromacs style.'''
    if 'HA2 GLY' in pdb_line:
        return pdb_line.replace('HA2 GLY', 'HA1 GLY')
    elif 'HA3 GLY' in pdb_line:
        return pdb_line.replace('HA3 GLY', 'HA2 GLY')
    elif 'HB2 SER' in pdb_line:
        return pdb_line.replace('HB2 SER', 'HB1 SER')
    elif 'HB3 SER' in pdb_line:
        return pdb_line.replace('HB3 SER', 'HB2 SER')
    elif 'HB2 LEU' in pdb_line:
        return pdb_line.replace('HB2 LEU', 'HB1 LEU')
    elif 'HB3 LEU' in pdb_line:
        return pdb_line.replace('HB3 LEU', 'HB2 LEU')
    elif 'HG12 ILE' in pdb_line:
        return pdb_line.replace('HG12 ILE', 'HG11 ILE')
    elif 'HG13 ILE' in pdb_line:
        return pdb_line.replace('HG13 ILE', 'HG12 ILE')
    elif 'HB2 ASN' in pdb_line:
        return pdb_line.replace('HB2 ASN', 'HB1 ASN')
    elif 'HB3 ASN' in pdb_line:
        return pdb_line.replace('HB3 ASN', 'HB2 ASN')
    elif 'HB2 GLN' in pdb_line:
        return pdb_line.replace('HB2 GLN', 'HB1 GLN')
    elif 'HB3 GLN' in pdb_line:
        return pdb_line.replace('HB3 GLN', 'HB2 GLN')
    elif 'HG2 GLN' in pdb_line:
        return pdb_line.replace('HG2 GLN', 'HG1 GLN')
    elif 'HG3 GLN' in pdb_line:
        return pdb_line.replace('HG3 GLN', 'HG2 GLN')
    elif 'HB2 ARG' in pdb_line:
        return pdb_line.replace('HB2 ARG', 'HB1 ARG')
    elif 'HB3 ARG' in pdb_line:
        return pdb_line.replace('HB3 ARG', 'HB2 ARG')
    elif 'HG2 ARG' in pdb_line:
        return pdb_line.replace('HG2 ARG', 'HG1 ARG')
    elif 'HG3 ARG' in pdb_line:
        return pdb_line.replace('HG3 ARG', 'HG2 ARG')
    elif 'HD2 ARG' in pdb_line:
        return pdb_line.replace('HD2 ARG', 'HD1 ARG')
    elif 'HD3 ARG' in pdb_line:
        return pdb_line.replace('HD3 ARG', 'HD2 ARG')
    elif 'HB2 HID' in pdb_line:
        return pdb_line.replace('HB2 HID', 'HB1 HID')
    elif 'HB3 HID' in pdb_line:
        return pdb_line.replace('HB3 HID', 'HB2 HID')
    elif 'HB2 HIE' in pdb_line:
        return pdb_line.replace('HB2 HIE', 'HB1 HIE')
    elif 'HB3 HIE' in pdb_line:
        return pdb_line.replace('HB3 HIE', 'HB2 HIE')
    elif 'HB2 HIP' in pdb_line:
        return pdb_line.replace('HB2 HIP', 'HB1 HIP')
    elif 'HB3 HIP' in pdb_line:
        return pdb_line.replace('HB3 HIP', 'HB2 HIP')
    elif 'HB2 HIS' in pdb_line:
        return pdb_line.replace('HB2 HIS', 'HB1 HIS')
    elif 'HB3 HIS' in pdb_line:
        return pdb_line.replace('HB3 HIS', 'HB2 HIS')
    elif 'HB2 TRP' in pdb_line:
        return pdb_line.replace('HB2 TRP', 'HB1 TRP')
    elif 'HB3 TRP' in pdb_line:
        return pdb_line.replace('HB3 TRP', 'HB2 TRP')
    elif 'HB2 PHE' in pdb_line:
        return pdb_line.replace('HB2 PHE', 'HB1 PHE')
    elif 'HB3 PHE' in pdb_line:
        return pdb_line.replace('HB3 PHE', 'HB2 PHE')
    elif 'HB2 TYR' in pdb_line:
        return pdb_line.replace('HB2 TYR', 'HB1 TYR')
    elif 'HB3 TYR' in pdb_line:
        return pdb_line.replace('HB3 TYR', 'HB2 TYR')
    elif 'HB2 ASP' in pdb_line:
        return pdb_line.replace('HB2 ASP', 'HB1 ASP')
    elif 'HB3 ASP' in pdb_line:
        return pdb_line.replace('HB3 ASP', 'HB2 ASP')
    elif 'HB2 GLU' in pdb_line:
        return pdb_line.replace('HB2 GLU', 'HB1 GLU')
    elif 'HB3 GLU' in pdb_line:
        return pdb_line.replace('HB3 GLU', 'HB2 GLU')
    elif 'HG2 GLU' in pdb_line:
        return pdb_line.replace('HG2 GLU', 'HG1 GLU')
    elif 'HG3 GLU' in pdb_line:
        return pdb_line.replace('HG3 GLU', 'HG2 GLU')
    elif 'HB2 LYS' in pdb_line:
        return pdb_line.replace('HB2 LYS', 'HB1 LYS')
    elif 'HB3 LYS' in pdb_line:
        return pdb_line.replace('HB3 LYS', 'HB2 LYS')
    elif 'HG2 LYS' in pdb_line:
        return pdb_line.replace('HG2 LYS', 'HG1 LYS')
    elif 'HG3 LYS' in pdb_line:
        return pdb_line.replace('HG3 LYS', 'HG2 LYS')
    elif 'HD2 LYS' in pdb_line:
        return pdb_line.replace('HD2 LYS', 'HD1 LYS')
    elif 'HD3 LYS' in pdb_line:
        return pdb_line.replace('HD3 LYS', 'HD2 LYS')
    elif 'HE2 LYS' in pdb_line:
        return pdb_line.replace('HE2 LYS', 'HE1 LYS')
    elif 'HE3 LYS' in pdb_line:
        return pdb_line.replace('HE3 LYS', 'HE2 LYS')
    elif 'HD2 PRO' in pdb_line:
        return pdb_line.replace('HD2 PRO', 'HD1 PRO')
    elif 'HD3 PRO' in pdb_line:
        return pdb_line.replace('HD3 PRO', 'HD2 PRO')
    elif 'HG2 PRO' in pdb_line:
        return pdb_line.replace('HG2 PRO', 'HG1 PRO')
    elif 'HG3 PRO' in pdb_line:
        return pdb_line.replace('HG3 PRO', 'HG2 PRO')
    elif 'HB2 PRO' in pdb_line:
        return pdb_line.replace('HB2 PRO', 'HB1 PRO')
    elif 'HB3 PRO' in pdb_line:
        return pdb_line.replace('HB3 PRO', 'HB2 PRO')
    elif 'HB2 CYS' in pdb_line:
        return pdb_line.replace('HB2 CYS', 'HB1 CYS')
    elif 'HB3 CYS' in pdb_line:
        return pdb_line.replace('HB3 CYS', 'HB2 CYS')
    elif 'HB2 MET' in pdb_line:
        return pdb_line.replace('HB2 MET', 'HB1 MET')
    elif 'HB3 MET' in pdb_line:
        return pdb_line.replace('HB3 MET', 'HB2 MET')
    elif 'HG2 MET' in pdb_line:
        return pdb_line.replace('HG2 MET', 'HG1 MET')
    elif 'HG3 MET' in pdb_line:
        return pdb_line.replace('HG3 MET', 'HG2 MET')
    else:
        return pdb_line

def Hmap_gro_gmx2omm(gro_line):
    '''Converts Gromacs style amber99sb-ildn Hydrogen names in Gromacs .top file to OpenMM style.'''
    if 'GLY    HA1' in gro_line:
        return re.sub('GLY    HA1', 'GLY    HA2', gro_line)
    elif 'GLY    HA2' in gro_line:
        return re.sub('GLY    HA2', 'GLY    HA3', gro_line)
    elif 'SER    HB1' in gro_line:
        return re.sub('SER    HB1', 'SER    HB2', gro_line)
    elif 'SER    HB2' in gro_line:
        return re.sub('SER    HB2', 'SER    HB3', gro_line)
    elif 'LEU    HB1' in gro_line:
        return re.sub('LEU    HB1', 'LEU    HB2', gro_line)
    elif 'LEU    HB2' in gro_line:
        return re.sub('LEU    HB2', 'LEU    HB3', gro_line)
    elif 'ILE   HG11' in gro_line:
        return re.sub('ILE   HG11', 'ILE   HG12', gro_line)
    elif 'ILE   HG12' in gro_line:
        return re.sub('ILE   HG12', 'ILE   HG13', gro_line)
    elif 'ASN    HB1' in gro_line:
        return re.sub('ASN    HB1', 'ASN    HB2', gro_line)
    elif 'ASN    HB2' in gro_line:
        return re.sub('ASN    HB2', 'ASN    HB3', gro_line)
    elif 'GLN    HB1' in gro_line:
        return re.sub('GLN    HB1', 'GLN    HB2', gro_line)
    elif 'GLN    HB2' in gro_line:
        return re.sub('GLN    HB2', 'GLN    HB3', gro_line)
    elif 'GLN    HG1' in gro_line:
        return re.sub('GLN    HG1', 'GLN    HG2', gro_line)
    elif 'GLN    HG2' in gro_line:
        return re.sub('GLN    HG2', 'GLN    HG3', gro_line)
    elif 'ARG    HB1' in gro_line:
        return re.sub('ARG    HB1', 'ARG    HB2', gro_line)
    elif 'ARG    HB2' in gro_line:
        return re.sub('ARG    HB2', 'ARG    HB3', gro_line)
    elif 'ARG    HG1' in gro_line:
        return re.sub('ARG    HG1', 'ARG    HG2', gro_line)
    elif 'ARG    HG2' in gro_line:
        return re.sub('ARG    HG2', 'ARG    HG3', gro_line)
    elif 'ARG    HD1' in gro_line:
        return re.sub('ARG    HD1', 'ARG    HD2', gro_line)
    elif 'ARG    HD2' in gro_line:
        return re.sub('ARG    HD2', 'ARG    HD3', gro_line)
    elif 'HID    HB1' in gro_line:
        return re.sub('HID    HB1', 'HID    HB2', gro_line)
    elif 'HID    HB2' in gro_line:
        return re.sub('HID    HB2', 'HID    HB3', gro_line)
    elif 'HIE    HB1' in gro_line:
        return re.sub('HIE    HB1', 'HIE    HB2', gro_line)
    elif 'HIE    HB2' in gro_line:
        return re.sub('HIE    HB2', 'HIE    HB3', gro_line)
    elif 'HIP    HB1' in gro_line:
        return re.sub('HIP    HB1', 'HIP    HB2', gro_line)
    elif 'HIP    HB2' in gro_line:
        return re.sub('HIP    HB2', 'HIP    HB3', gro_line)
    elif 'HIS    HB1' in gro_line:
        return re.sub('HIS    HB1', 'HIS    HB2', gro_line)
    elif 'HIS    HB2' in gro_line:
        return re.sub('HIS    HB2', 'HIS    HB3', gro_line)
    elif 'TRP    HB1' in gro_line:
        return re.sub('TRP    HB1', 'TRP    HB2', gro_line)
    elif 'TRP    HB2' in gro_line:
        return re.sub('TRP    HB2', 'TRP    HB3', gro_line)
    elif 'PHE    HB1' in gro_line:
        return re.sub('PHE    HB1', 'PHE    HB2', gro_line)
    elif 'PHE    HB2' in gro_line:
        return re.sub('PHE    HB2', 'PHE    HB3', gro_line)
    elif 'TYR    HB1' in gro_line:
        return re.sub('TYR    HB1', 'TYR    HB2', gro_line)
    elif 'TYR    HB2' in gro_line:
        return re.sub('TYR    HB2', 'TYR    HB3', gro_line)
    elif 'ASP    HB1' in gro_line:
        return re.sub('ASP    HB1', 'ASP    HB2', gro_line)
    elif 'ASP    HB2' in gro_line:
        return re.sub('ASP    HB2', 'ASP    HB3', gro_line)
    elif 'GLU    HB1' in gro_line:
        return re.sub('GLU    HB1', 'GLU    HB2', gro_line)
    elif 'GLU    HB2' in gro_line:
        return re.sub('GLU    HB2', 'GLU    HB3', gro_line)
    elif 'GLU    HG1' in gro_line:
        return re.sub('GLU    HG1', 'GLU    HG2', gro_line)
    elif 'GLU    HG2' in gro_line:
        return re.sub('GLU    HG2', 'GLU    HG3', gro_line)
    elif 'LYS    HB1' in gro_line:
        return re.sub('LYS    HB1', 'LYS    HB2', gro_line)
    elif 'LYS    HB2' in gro_line:
        return re.sub('LYS    HB2', 'LYS    HB3', gro_line)
    elif 'LYS    HG1' in gro_line:
        return re.sub('LYS    HG1', 'LYS    HG2', gro_line)
    elif 'LYS    HG2' in gro_line:
        return re.sub('LYS    HG2', 'LYS    HG3', gro_line)
    elif 'LYS    HD1' in gro_line:
        return re.sub('LYS    HD1', 'LYS    HD2', gro_line)
    elif 'LYS    HD2' in gro_line:
        return re.sub('LYS    HD2', 'LYS    HD3', gro_line)
    elif 'LYS    HE1' in gro_line:
        return re.sub('LYS    HE1', 'LYS    HE2', gro_line)
    elif 'LYS    HE2' in gro_line:
        return re.sub('LYS    HE2', 'LYS    HE3', gro_line)
    elif 'PRO    HD1' in gro_line:
        return re.sub('PRO    HD1', 'PRO    HD2', gro_line)
    elif 'PRO    HD2' in gro_line:
        return re.sub('PRO    HD2', 'PRO    HD3', gro_line)
    elif 'PRO    HG1' in gro_line:
        return re.sub('PRO    HG1', 'PRO    HG2', gro_line)
    elif 'PRO    HG2' in gro_line:
        return re.sub('PRO    HG2', 'PRO    HG3', gro_line)
    elif 'PRO    HB1' in gro_line:
        return re.sub('PRO    HB1', 'PRO    HB2', gro_line)
    elif 'PRO    HB2' in gro_line:
        return re.sub('PRO    HB2', 'PRO    HB3', gro_line)
    elif 'CYS    HB1' in gro_line:
        return re.sub('CYS    HB1', 'CYS    HB2', gro_line)
    elif 'CYS    HB2' in gro_line:
        return re.sub('CYS    HB2', 'CYS    HB3', gro_line)
    elif 'MET    HB1' in gro_line:
        return re.sub('MET    HB1', 'MET    HB2', gro_line)
    elif 'MET    HB2' in gro_line:
        return re.sub('MET    HB2', 'MET    HB3', gro_line)
    elif 'MET    HG1' in gro_line:
        return re.sub('MET    HG1', 'MET    HG2', gro_line)
    elif 'MET    HG2' in gro_line:
        return re.sub('MET    HG2', 'MET    HG3', gro_line)
    else:
        return gro_line


def Hmap_top_gmx2omm(top_line):
    '''Converts Gromacs style amber99sb-ildn Hydrogen names in Gromacs .top file to OpenMM style.'''
    if 'GLY    HA1' in top_line:
        return re.sub('GLY    HA1', 'GLY    HA2', top_line)
    elif 'GLY    HA2' in top_line:
        return re.sub('GLY    HA2', 'GLY    HA3', top_line)
    elif 'SER    HB1' in top_line:
        return re.sub('SER    HB1', 'SER    HB2', top_line)
    elif 'SER    HB2' in top_line:
        return re.sub('SER    HB2', 'SER    HB3', top_line)
    elif 'LEU    HB1' in top_line:
        return re.sub('LEU    HB1', 'LEU    HB2', top_line)
    elif 'LEU    HB2' in top_line:
        return re.sub('LEU    HB2', 'LEU    HB3', top_line)
    elif 'ILE   HG11' in top_line:
        return re.sub('ILE   HG11', 'ILE   HG12', top_line)
    elif 'ILE   HG12' in top_line:
        return re.sub('ILE   HG12', 'ILE   HG13', top_line)
    elif 'ASN    HB1' in top_line:
        return re.sub('ASN    HB1', 'ASN    HB2', top_line)
    elif 'ASN    HB2' in top_line:
        return re.sub('ASN    HB2', 'ASN    HB3', top_line)
    elif 'GLN    HB1' in top_line:
        return re.sub('GLN    HB1', 'GLN    HB2', top_line)
    elif 'GLN    HB2' in top_line:
        return re.sub('GLN    HB2', 'GLN    HB3', top_line)
    elif 'GLN    HG1' in top_line:
        return re.sub('GLN    HG1', 'GLN    HG2', top_line)
    elif 'GLN    HG2' in top_line:
        return re.sub('GLN    HG2', 'GLN    HG3', top_line)
    elif 'ARG    HB1' in top_line:
        return re.sub('ARG    HB1', 'ARG    HB2', top_line)
    elif 'ARG    HB2' in top_line:
        return re.sub('ARG    HB2', 'ARG    HB3', top_line)
    elif 'ARG    HG1' in top_line:
        return re.sub('ARG    HG1', 'ARG    HG2', top_line)
    elif 'ARG    HG2' in top_line:
        return re.sub('ARG    HG2', 'ARG    HG3', top_line)
    elif 'ARG    HD1' in top_line:
        return re.sub('ARG    HD1', 'ARG    HD2', top_line)
    elif 'ARG    HD2' in top_line:
        return re.sub('ARG    HD2', 'ARG    HD3', top_line)
    elif 'HID    HB1' in top_line:
        return re.sub('HID    HB1', 'HID    HB2', top_line)
    elif 'HID    HB2' in top_line:
        return re.sub('HID    HB2', 'HID    HB3', top_line)
    elif 'HIE    HB1' in top_line:
        return re.sub('HIE    HB1', 'HIE    HB2', top_line)
    elif 'HIE    HB2' in top_line:
        return re.sub('HIE    HB2', 'HIE    HB3', top_line)
    elif 'HIP    HB1' in top_line:
        return re.sub('HIP    HB1', 'HIP    HB2', top_line)
    elif 'HIP    HB2' in top_line:
        return re.sub('HIP    HB2', 'HIP    HB3', top_line)
    elif 'HIS    HB1' in top_line:
        return re.sub('HIS    HB1', 'HIS    HB2', top_line)
    elif 'HIS    HB2' in top_line:
        return re.sub('HIS    HB2', 'HIS    HB3', top_line)
    elif 'TRP    HB1' in top_line:
        return re.sub('TRP    HB1', 'TRP    HB2', top_line)
    elif 'TRP    HB2' in top_line:
        return re.sub('TRP    HB2', 'TRP    HB3', top_line)
    elif 'PHE    HB1' in top_line:
        return re.sub('PHE    HB1', 'PHE    HB2', top_line)
    elif 'PHE    HB2' in top_line:
        return re.sub('PHE    HB2', 'PHE    HB3', top_line)
    elif 'TYR    HB1' in top_line:
        return re.sub('TYR    HB1', 'TYR    HB2', top_line)
    elif 'TYR    HB2' in top_line:
        return re.sub('TYR    HB2', 'TYR    HB3', top_line)
    elif 'ASP    HB1' in top_line:
        return re.sub('ASP    HB1', 'ASP    HB2', top_line)
    elif 'ASP    HB2' in top_line:
        return re.sub('ASP    HB2', 'ASP    HB3', top_line)
    elif 'GLU    HB1' in top_line:
        return re.sub('GLU    HB1', 'GLU    HB2', top_line)
    elif 'GLU    HB2' in top_line:
        return re.sub('GLU    HB2', 'GLU    HB3', top_line)
    elif 'GLU    HG1' in top_line:
        return re.sub('GLU    HG1', 'GLU    HG2', top_line)
    elif 'GLU    HG2' in top_line:
        return re.sub('GLU    HG2', 'GLU    HG3', top_line)
    elif 'LYS    HB1' in top_line:
        return re.sub('LYS    HB1', 'LYS    HB2', top_line)
    elif 'LYS    HB2' in top_line:
        return re.sub('LYS    HB2', 'LYS    HB3', top_line)
    elif 'LYS    HG1' in top_line:
        return re.sub('LYS    HG1', 'LYS    HG2', top_line)
    elif 'LYS    HG2' in top_line:
        return re.sub('LYS    HG2', 'LYS    HG3', top_line)
    elif 'LYS    HD1' in top_line:
        return re.sub('LYS    HD1', 'LYS    HD2', top_line)
    elif 'LYS    HD2' in top_line:
        return re.sub('LYS    HD2', 'LYS    HD3', top_line)
    elif 'LYS    HE1' in top_line:
        return re.sub('LYS    HE1', 'LYS    HE2', top_line)
    elif 'LYS    HE2' in top_line:
        return re.sub('LYS    HE2', 'LYS    HE3', top_line)
    elif 'PRO    HD1' in top_line:
        return re.sub('PRO    HD1', 'PRO    HD2', top_line)
    elif 'PRO    HD2' in top_line:
        return re.sub('PRO    HD2', 'PRO    HD3', top_line)
    elif 'PRO    HG1' in top_line:
        return re.sub('PRO    HG1', 'PRO    HG2', top_line)
    elif 'PRO    HG2' in top_line:
        return re.sub('PRO    HG2', 'PRO    HG3', top_line)
    elif 'PRO    HB1' in top_line:
        return re.sub('PRO    HB1', 'PRO    HB2', top_line)
    elif 'PRO    HB2' in top_line:
        return re.sub('PRO    HB2', 'PRO    HB3', top_line)
    elif 'CYS    HB1' in top_line:
        return re.sub('CYS    HB1', 'CYS    HB2', top_line)
    elif 'CYS    HB2' in top_line:
        return re.sub('CYS    HB2', 'CYS    HB3', top_line)
    elif 'MET    HB1' in top_line:
        return re.sub('MET    HB1', 'MET    HB2', top_line)
    elif 'MET    HB2' in top_line:
        return re.sub('MET    HB2', 'MET    HB3', top_line)
    elif 'MET    HG1' in top_line:
        return re.sub('MET    HG1', 'MET    HG2', top_line)
    elif 'MET    HG2' in top_line:
        return re.sub('MET    HG2', 'MET    HG3', top_line)
    else:
        return top_line

em_mdp_contents = '''; VARIOUS PREPROCESSING OPTIONS =
title                    =
cpp                      = /lib/cpp
include                  =
define                   =

; RUN CONTROL PARAMETERS =
integrator               = steep
; start time and timestep in ps =
tinit                    = 0.0
nsteps                   = 100
; number of steps for center of mass motion removal =
nstcomm                  = 1
; ENERGY MINIMIZATION OPTIONS
; Force tolerance and initial step-size
emtol                    = 10
emstep                   = 0.01
; Max number of iterations in relax_shells
niter                    = 20
; Step size (1/ps^2) for minimization of flexible constraints
fcstep                   = 0
; Frequency of steepest descents steps when doing CG
nstcgsteep               = 100
nbfgscorr                = 10'''

top_hoh_text = '''[ moleculetype ]
; molname   nrexcl
HOH     2

[ atoms ]
; id  at type     res nr  res name  at name  cg nr  charge    mass
  1   OW          1       HOH       O        1      -0.834    16.00000
  2   HW          1       HOH       H1       1       0.417     1.00800
  3   HW          1       HOH       H2       1       0.417     1.00800

#ifndef FLEXIBLE

[ settles ]
; OW    funct   doh dhh
1       1       0.09572 0.15139

[ exclusions ]
1   2   3
2   1   3
3   1   2

#else

[ bonds ]
; i     j       funct   length  force_constant
1       2       1       0.09572 502416.0   0.09572        502416.0 
1       3       1       0.09572 502416.0   0.09572        502416.0 
        

[ angles ]
; i     j       k       funct   angle   force_constant
2       1       3       1       104.52  628.02      104.52  628.02  

#endif

'''

eq_mdp_contents='''; VARIOUS PREPROCESSING OPTIONS = 
title                    = 
cpp                      = /lib/cpp
include                  = 
define                   =

; RUN CONTROL PARAMETERS = 
integrator               = sd
; start time and timestep in ps = 
tinit                    = 0.0
dt                       = 0.002
nsteps                   = 50000
; number of steps for center of mass motion removal = 
nstcomm                  = 1


; OUTPUT CONTROL OPTIONS = 
; Output frequency for coords (x), velocities (v) and forces (f) = 
nstxout                  = 100000
nstvout                  = 100000
nstfout                  = 0
; Output frequency for energies to log file and energy file = 
nstlog                   = 10000
nstenergy                = 500
; Output frequency and precision for xtc file = 
nstxtcout                = 10000
xtc_precision            = 10000
; This selects the subset of atoms for the xtc file. You can = 
; select multiple groups. By default all atoms will be written. = 
xtc-grps                 = 
; Selection of energy groups = 
energygrps               = Protein Water_and_ions

; NEIGHBORSEARCHING PARAMETERS = 
; nblist update frequency = 
nstlist                  = 5
; ns algorithm (simple or grid) = 
ns_type                  = grid
; Periodic boundary conditions: xyz or none = 
pbc                      = xyz
; nblist cut-off         = 
rlist                    = 1.0
domain-decomposition     = no

; OPTIONS FOR ELECTROSTATICS AND VDW = 
; Method for doing electrostatics = 
coulombtype              = PME
; rcoulomb_switch          = 0.0 (not used with coulombtype=PME)
rcoulomb                 = 1.0
; Dielectric constant (DC) for cut-off or DC of reaction field = 
epsilon_r                = 1 
; Method for doing Van der Waals = 
vdw_type                 = Switch
; cut-off lengths        = 
rvdw_switch              = 0.95
rvdw                     = 1.0
; Apply long range dispersion corrections for Energy and Pressure = 
DispCorr                 = No
; Spacing for the PME/PPPM FFT grid = 
fourierspacing           = 0.12
; FFT grid size, when a value is 0 fourierspacing will be used = 
fourier_nx               = 10
fourier_ny               = 10
fourier_nz               = 10
; EWALD/PME/PPPM parameters = 
pme_order                = 4
ewald_rtol               = 1e-05
epsilon_surface          = 0
optimize_fft             = no

; OPTIONS FOR WEAK COUPLING ALGORITHMS = 
; Temperature coupling   = 
tcoupl                   = no
; NOTE: using Langevin dynamics as thermostat
; Groups to couple separately = 
tc-grps                  = Protein Water_and_ions
; Time constant (ps) and reference temperature (K) = 
; This is the inverse of the damping coefficient (in ps)
tau_t                    = 0.2 0.2
ref_t                    = 300 300
; ld_seed = -1 means that the seed is calculated from the process ID
ld_seed                  = -1
; Pressure coupling      = 
Pcoupl                   = parrinello-rahman
; Pressure coupling type - isotropic,semiisotropic,anisotropic,surface-tension,triclinic
Pcoupltype               = isotropic
; Time constant (ps), compressibility (1/bar) and reference P (bar) = 
tau-p                    = 1.0 1.0 
compressibility          = 4.5E-5 4.5E-5
ref-p                    = 1.0 1.0 


; GENERATE VELOCITIES FOR STARTUP RUN = 
gen_vel                  = yes
gen_temp                 = 300
gen_seed                 = 591

; OPTIONS FOR BONDS     = 
constraints              = hbonds 
; Type of constraint algorithm = 
constraint_algorithm     = Lincs
; Do not constrain the start configuration = 
unconstrained_start      = no
; Relative tolerance of shake = 
shake_tol                = 0.0001
; Highest order in the expansion of the constraint coupling matrix = 
lincs_order              = 4
; Lincs will write a warning to the stderr if in one step a bond = 
; rotates over more degrees than = 
lincs_warnangle          = 30
; Convert harmonic bonds to morse potentials = 
morse                    = no

'''

prod_mdp_contents='''; VARIOUS PREPROCESSING OPTIONS = 
title                    = 
cpp                      = /lib/cpp
include                  = 
define                   =

; RUN CONTROL PARAMETERS = 
integrator               = sd
; start time and timestep in ps = 
tinit                    = 0.0
dt                       = 0.002
nsteps                   = 5000000
; number of steps for center of mass motion removal = 
nstcomm                  = 1


; OUTPUT CONTROL OPTIONS = 
; Output frequency for coords (x), velocities (v) and forces (f) = 
nstxout                  = 100000
nstvout                  = 100000
nstfout                  = 0
; Output frequency for energies to log file and energy file = 
nstlog                   = 10000
nstenergy                = 500
; Output frequency and precision for xtc file = 
nstxtcout                = 10000
xtc_precision            = 10000
; This selects the subset of atoms for the xtc file. You can = 
; select multiple groups. By default all atoms will be written. = 
xtc-grps                 = 
; Selection of energy groups = 
energygrps               = Protein Water_and_ions

; NEIGHBORSEARCHING PARAMETERS = 
; nblist update frequency = 
nstlist                  = 10
; ns algorithm (simple or grid) = 
ns_type                  = grid
; Periodic boundary conditions: xyz or none = 
pbc                      = xyz
; nblist cut-off         = 
rlist                    = 1.0
domain-decomposition     = no

; OPTIONS FOR ELECTROSTATICS AND VDW = 
; Method for doing electrostatics = 
coulombtype              = PME
; rcoulomb_switch          = 0.0 (not used with coulombtype=PME)
rcoulomb                 = 1.0
; Dielectric constant (DC) for cut-off or DC of reaction field = 
epsilon_r                = 1 
; Method for doing Van der Waals = 
vdw_type                 = Shift
; cut-off lengths        = 
rvdw_switch              = 0.95
rvdw                     = 1.0
; Apply long range dispersion corrections for Energy and Pressure = 
DispCorr                 = No
; Spacing for the PME/PPPM FFT grid = 
fourierspacing           = 0.12
; FFT grid size, when a value is 0 fourierspacing will be used = 
fourier_nx               = 10
fourier_ny               = 10
fourier_nz               = 10
; EWALD/PME/PPPM parameters = 
pme_order                = 4
ewald_rtol               = 1e-05
epsilon_surface          = 0
optimize_fft             = no

; OPTIONS FOR WEAK COUPLING ALGORITHMS = 
; Temperature coupling   = 
tcoupl                   = no
; NOTE: using Langevin dynamics as thermostat
; Groups to couple separately = 
tc-grps                  = Protein Water_and_ions
; Time constant (ps) and reference temperature (K) = 
; This is the inverse of the damping coefficient (in ps)
tau_t                    = 0.2 0.2
ref_t                    = 300 300
; ld_seed = -1 means that the seed is calculated from the process ID
ld_seed                  = -1
; Pressure coupling      = 
Pcoupl                   = parrinello-rahman
; Pressure coupling type - isotropic,semiisotropic,anisotropic,surface-tension,triclinic
Pcoupltype               = isotropic
; Time constant (ps), compressibility (1/bar) and reference P (bar) = 
tau-p                    = 1.0 1.0 
compressibility          = 4.5E-5 4.5E-5
ref-p                    = 1.0 1.0 


; GENERATE VELOCITIES FOR STARTUP RUN = 
gen_vel                  = yes
gen_temp                 = 300
gen_seed                 = 591

; OPTIONS FOR BONDS     = 
constraints              = hbonds 
; Type of constraint algorithm = 
constraint_algorithm     = Lincs
; Do not constrain the start configuration = 
unconstrained_start      = no
; Relative tolerance of shake = 
shake_tol                = 0.0001
; Highest order in the expansion of the constraint coupling matrix = 
lincs_order              = 4
; Lincs will write a warning to the stderr if in one step a bond = 
; rotates over more degrees than = 
lincs_warnangle          = 30
; Convert harmonic bonds to morse potentials = 
morse                    = no

'''

