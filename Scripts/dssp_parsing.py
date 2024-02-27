"""
Code for running and parsing DSSP secondary structure prediction analysis. Port 
the results into a MN object's attribute to redraw cartoon representations. 
"""

# PREAMBLE

# FUNCTIONS
def run_dssp(dssp_bin_path, structure_file_path, output_file = 'mkdssp_output.cif', working_dir = './'):
    """
    run the user-designated mkdssp command on the input structure file. output
    written to mkdssp_output.cif or user-defined file. 

    :parameter dssp_bin_path: string, global or local path to the mkdssp 
        executable file.
    :parameter structure_file_path: string, global or local path to the 
        structure file used as input to mkdssp. format: pdb or cif file types; 
        need to debug file format issues. 
    :parameter output_file: string, local or global path for results from mkdssp
        to be written. Extension of this file name is important. If '.dssp',
        then output file will have the original dssp file format, which should
        be avoided due to character length limits for fields. This will raise
        issues for structures with large number of residues.
        If anything else (including '.cif'), then the new mmcif file
        format will be used.
        default: 'mkdssp_output.cif'
    :parameter working_dir: string, local or global path to designate where the
        working directory should be. 
    :return output_file or 1: return value depends on success of the mkdssp run.
        if the mkdssp command errored out or returns a non-zero returncode, then
        return 1 (int). 
        else, return the string file path written by mkdssp.

    """
    import subprocess

    # try running the command
    try:
        stdout = subprocess.run(f'{dssp_bin_path} {structure_file_path} {output_file}',capture_output=True,shell=True,cwd=working_dir,check=True)
    # if the subprocess.run errors with a CalledProcessError, then the mkdssp 
    # call failed. return 1 as the exit status, denoting the program failure
    except subprocess.CalledProcessError as e:
        print(f'DSSP calculation failed with error: {e}')
        return 1

    # similarly, if the returncode != 0, then some other error occurred and 
    # we should not attempt to parse the output file for 2ndary structure info
    # so return 1 as the exit status, denoting the error/failure
    if stdout.returncode != 0:
        return 1
    
    return output_file


def parse_dssp_output(results):
    """ 
    NEVER FINISHED. NOT NEEDED SINCE THE DSSP FILE FORMAT SHOULD BE AVOIDED; it
    suffers from similar issues as the PDB file format, hard-set column issues
    
    read the input dssp-formatted results from mkdssp and gather the 2ndary 
    structure labels associated with the residue index
    :parameter results: string, mkdssp results to be parsed
    :return struct_dict: dictionary, keys have {resid}_{chainid} format and 
        values are the assigned 2ndary structure string
    """
    import re

    struct_dict = {}
    flag = False
    for line in results.split('\n'):
        if re.search('#', line):
            flag = True
            continue
        if flag:
            resnum = line[0:5].strip()
            chainid= line[10:12].strip()
            struct = line[14:17].strip()
            struct_dict[f'{resnum}_{chainid}'] = struct

    return struct_dict


def parse_mmcif_output(result_file_path):
    """
    read the mmcif-formatted results from mkdssp and gather the 2ndary structure
    labels associated with the residue index
    :parameter results: string, path to the temporary mmcif file written to 
        storage since biotite does not seem to be able to handle file-like 
        objects
    :return struct_dict: numpy array, size (nAtoms,) of ints where values are
        the assigned 2ndary structure integers as defined in the 
        struct_type_dict.
    """
    import numpy as np
    import biotite.structure.io.pdbx as pdbx

    # define mmcif 2ndary structure dictionary
    # list of possible conf_type_id values: 
    # https://manpages.debian.org/bookworm/dssp/mkdssp.1.en.html
    # integers used to denote these types of structures, reorganized a bit but
    # won't lose details about structures like the current implementation in MN 
    # does. 
    # NOTE: the mmcif file format output from mkdssp does not differentiate a 
    # \beta-bridge residue and a beta-sheet/strand residue. The original dssp 
    # output format did differentiate the two. No great loss in my opinion,
    # still gonna denote the loss by skipping the index. 
    # 0 --> unstructured/loop
    # 1 --> H, HELX_RH_AL_P, \alpha helix
    # 2 --> G, HELX_RH_3T_P, 3_{10} helix
    # 3 --> I, HELX_RH_PI_P, \pi helix
    # 4 --> P, HELX_LH_PP_P, Helix_PPII 
    # 5 --> B, STRN, \beta-bridge residue
    # 5 --> E, STRN, strand
    # 7 --> T, TURN_TY1_P, turn
    # 8 --> S, BEND, bend
    struct_type_dict = {'HELX_RH_AL_P': 1,
                        'HELX_RH_3T_P': 2,
                        'HELX_RH_PI_P': 3,
                        'HELX_LH_PP_P': 4,
                        'STRN'        : 5,
                        'TURN_TY1_P'  : 7,
                        'BEND'        : 8}

    # read mmcif file
    cif_file = pdbx.PDBxFile.read(result_file_path)
    # grab the dictionary associated with the _struct_conf field in the mmcif 
    # file
    sec_struct = cif_file.get_category('struct_conf')
    # make sure its not an empty dict
    if not sec_struct:
        print('No 2ndary structure information found in {result_file_path} file.')
        return {}

    # sec_struct is a dict with keys associated with each column in the section
    # we are interested in 2ndary structure id and type as well as residue start 
    # and end numbers, chain id. 
    structures = list(zip(sec_struct['id'],
                          sec_struct['conf_type_id'],
                          sec_struct['beg_auth_asym_id'],
                          sec_struct['beg_auth_seq_id'],
                          sec_struct['end_auth_seq_id']))

    resids = []
    chains = []
    structures_list = []
    for structure in structures:
        res_range = range(int(structure[3]),int(structure[4])+1)
        resids += list(res_range)
        chains.append(structure[2])
        structures_list.append([structure[0],   # dssp-assigned 2ndary structure id
                                structure[1],   # type of 2ndary structure
                                structure[2],   # chainid
                                res_range])     # resid range

    resids = set(resids)
    chains = set(chains)

    # create an array of zeros to be filled with 2ndary structure integers
    nAtoms = int(cif_file.get_category('atom_site')['id'][-1])
    sec_struct_attribute = np.zeros(nAtoms, dtype=int)
    
    # Fill the sec_struct_attribute array with integers associated with 2ndary
    # structure indices (defined in struct_type_dict) as listed in the 
    # structures. 
    # loop over all atom i's, (chain id, resid) to match to results in the 
    # structures object 
    for i, (chainid, resid) in enumerate(zip(cif_file.get_category('atom_site')['label_asym_id'], cif_file.get_category('atom_site')['label_seq_id'].astype(int))):
        
        # quick and dirty check to see if the atom's chainid and resid are in 
        # the respective sets. if not, move on to the next atom.
        if chainid not in chains or resid not in resids:
            continue

        # loop over structures to determine the 2ndary structure type
        for structure in structures_list:
            # if chainid and resid match with a structure's, update the 
            # struct_type_dict element with the appropriate integer value
            if chainid == structure[2] and resid in structure[3]:
                sec_struct_attribute[i] = struct_type_dict[structure[1]]
                continue

    return sec_struct_attribute


##########################################
# test harness for functions defined here
##########################################
if __name__ == '__main__': 
    # run the run_dssp function; this test expects the mkdssp executable to be
    # in the path so that it can be ran without a explicit path
    output_file = run_dssp('mkdssp','./dssp_tests/model_structure.pdb',output_file='./dssp_tests/test_dssp.dssp')

    # read in the new dssp file to compare with the provided old file
    with open(output_file,'r') as new_file:
        new_dssp = new_file.read()

    with open('./dssp_tests/dssp.dssp','r') as old_file:
        old_dssp = old_file.read()

    # compare the function's file versus the file written directly from mkdssp
    # HEADER line format is strangely different. I can't immediately explain why
    # all other lines are exact matches.
    assert new_dssp[384:] == old_dssp[384:]
    print('Done checking the dssp file format')

    # run the mkdssp command with output formatted using the mmcif file format
    # can't 
    output_file = run_dssp('mkdssp','./dssp_tests/model_structure.pdb',output_file='./dssp_tests/test_dssp.cif')

    # load the two cif files into biotite structure objects, grab important 
    # 2ndary structure information dict and compare the two.
    import biotite.structure.io.pdbx as pdbx
    import numpy as np

    new_dssp_cif = pdbx.PDBxFile.read(output_file)
    new_sec_struct = new_dssp_cif.get_category('struct_conf')['conf_type_id']

    old_dssp_cif = pdbx.PDBxFile.read('./dssp_tests/dssp.cif')
    old_sec_struct = old_dssp_cif.get_category('struct_conf')['conf_type_id']

    assert np.all(new_sec_struct == old_sec_struct)
    print('Done checking the mmcif file format')

    # now parse the original_dssp file to gather the 2ndary structure attribute
    # values
    model_ss_values = parse_mmcif_output('./dssp_tests/dssp.cif')

    with open('./dssp_tests/ss_mapping.dat','w') as out:
        out.write(''.join([str(i)+'\n' for i in np.array(model_ss_values)]))



