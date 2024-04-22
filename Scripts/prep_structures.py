
import MDAnalysis
import numpy as np
import json
from pathlib import Path
import shutil
import re
import pickle
import sys

af_dirs = list(Path(sys.argv[1]).glob('*/AlphaFold/'))
outputDir = sys.argv[2]

af_dirs.sort()

order = {}

for af_dir in af_dirs:
    structure_file = list(af_dir.glob('*min_00.pdb'))[0]
    ranking_file = list(af_dir.glob('ranking_combined.json'))[0]
    
    # parse directory string (e.g. Sphmag02G125400.1)
    try:
        chromosomeID, geneID, temp = re.findall('(\d+)',af_dir.parent.stem)
    except:
        # unnumbered chromosomes are handled here
        chromosomeID = 'U'
        geneID = re.findall('(\d+)',af_dir.parent.stem)[0]

    # get the subdirectory or make one
    chromosomeDict = order.get(chromosomeID,{})
    
    with open(ranking_file, 'r') as json_in:
        geneData = json.load(json_in)
        modelID  = geneData['order'][0]
        pTMs     = geneData['ptms'][modelID]
        pLDDTs   = geneData['plddts'][modelID]

    chromosomeDict[geneID] = [pTMs, pLDDTs]

    order[chromosomeID] = chromosomeDict

    shutil.move(ranking_file, outputDir + 'ranking_jsons/' + ranking_file.parent.parent.stem + '.json')

    # load structure_file into MDA
    u = MDAnalysis.Universe(str(structure_file))
    u_all = u.select_atoms('all')
    # center the structure to (0,0,0)
    u_all.positions -= u_all.center_of_mass()
    positions = u_all.positions
    # calculate the princ components and align structure to those
    U,S,Vt = np.linalg.svd(positions)
    new_positions = positions @ Vt.T
    u_all.positions = new_positions
    # write centered, rotated structure to new directory
    u_all.write(outputDir + 'AlphaFold_structures/' + structure_file.parent.parent.name + '.pdb')

    print(f'Done working on {structure_file.parent.parent.name}')

with open('order.pkl','wb') as pkl_out:
    pickle.dump(order, pkl_out)

