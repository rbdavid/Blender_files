
import molecularnodes as mn
import bpy
import pickle
from pathlib import Path

# load order pickle file
# also 
# annotation/results pickle file? to add those to properties?

with open('C:\\Users\\rbdch\\Dropbox\\for_blender\\SAFA_visualizations\\sdiv_gene_order.pkl','rb') as pkl_in:
    order = pickle.load(pkl_in)

structure_file_path = 'D:\\BlenderFiles\\SAFA_visualization_work\\StructureFiles\\'

x_boundaries = (0,2)  # units of meters
y_init = 0  # units of rows
z_init = 0  # units of rows

# track the previous grid's 1/2 max of x-dimension to get ideal placement of next object
previous_halfmax_x = 0.
# track the previous row's 1/2 max of y dimensions to get ideal spacing between rows
previous_halfmax_y = 0.
# track the previous row's 1/2 max of z dimensions to get ideal spacing between rows
previous_halfmax_z = 0.
# set the spacing between structures
delta_scaling = 1.05
# set the blender object scaling transformation
size_scaling = 0.2
# create collector for y-dimension values
rows_y_dims = []
rows_z_dims = []

# loop over all entries in the dictionary
chromosomeIDs = [key for key in order.keys()]
chromosomeIDs.sort()
for chromosomeID in chromosomeIDs:
    # loop over all proteins associated with chromosomeID, sorted by pTMs values
    proteinIDs = list(order[chromosomeID].items())
    proteinIDs.sort(key = lambda x: x[1][0], reverse=True)
    for (protein, scores_list) in proteinIDs[:50]:
        if chromosomeID == 'U': 
            name = f'Sphmag{chromosomeID}{protein}'
        else: 
            name = f'Sphmag{chromosomeID}G{protein}'
        structure_file = structure_file_path + name + '.1.pdb'
        # load molecule object into blender
        mol = mn.io.local.load(structure_file, name=name, style='customRBD')
        b_object = bpy.data.objects[name]
        # scale the object down by a constant amount
        b_object.scale[0] = size_scaling
        b_object.scale[1] = size_scaling
        b_object.scale[2] = size_scaling

        # get the in-canvas dimensions of the object
        dimensions = b_object.dimensions*size_scaling
        # translate the object to its appropriate grid position
        # we know starting location will be (0,0,0) and then we want to span 
        # structures out along the x-axis until a max value is hit, then jump
        # down a grid space to the next row that starts at (0,y,0)
        
        # check to see if we've reached the soft boundary where we jump to the 
        # next row
        if previous_halfmax_x > x_boundaries[1]:
            previous_halfmax_x  = x_boundaries[0]
            previous_halfmax_z -= delta_scaling*max(rows_z_dims)
            rows_z_dims = []
       
        # if the object isn't the first in a row
        if previous_halfmax_x != x_boundaries[0]:
            # add first half of object's x dim to the location
            b_object.location[0] += previous_halfmax_x + dimensions[0]/2.
        else:
            b_object.location[0] = x_boundaries[0]
        # add second half of object's x dim to the previous_halfmax_x
        previous_halfmax_x = b_object.location[0] + dimensions[0]/2.
        
        b_object.location[1] += previous_halfmax_y
        rows_y_dims.append(dimensions[1]/2.)

        b_object.location[2] += previous_halfmax_z
        rows_z_dims.append(dimensions[2]/2.)

        # add metadata as properties of the structure object
        b_object['pLDDT'] = scores_list[1]
        b_object['pTMs']  = scores_list[0]
        b_object['Chromosome']  = str(chromosomeID)

        #b_object.modifiers['MolecularNodes'].node_group.nodes['MN_color_attribute_random.001'].inputs[3].default_value = 21
    
    previous_halfmax_x  = x_boundaries[0]
    previous_halfmax_y -= max(rows_y_dims)*delta_scaling
    previous_halfmax_z  = 0
    rows_y_dims = []

