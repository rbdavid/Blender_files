"""
functions to align objects to a standardized view point and then zoom in on them
"""

# PREAMBLE
import bpy
import numpy as np

def svd_align(vertices):
    """
    remove CoM translation, rotate the vertices to align principal components of
    vertices along the xyz dimensions, add back CoG translation, return vertices
    and the applied rotation matrix.

    use singular value decomposition to calculate the rotation matrix for 
    aligning the vertices to the xyz dimensions

    :parameter vertices: 2d array, assumed to be Nx3 associated with the atomic
        positions of a single structure. 
    :return new_vertices: 2d array, same shape as input vertices but now rotated.
    :return new_vertices: 2d array, shape 3x3, rotation matrix used to calculate
        the new vertices.
    """
    import numpy as np
    
    # calc the center of geometry
    CoG = np.mean(vertices,axis=0)
    # remove center of geometry from the vertices
    vertices = vertices - CoG
    # do the SVD to get the unitary array needed to rotate the vertices so that
    # the principal components are aligned with the xyz dimensions. 
    U, S, Vt = np.linalg.svd(vertices)
    # multiply matrices to do the alignment
    new_vertices = vertices @ Vt.T
    # add the center of geometry back to get original origin
    new_vertices = new_vertices + CoG
    
    return new_vertices, Vt

    """
    # getting the structure object
    structure = bpy.data.objects[object_name]
    # getting the structure object's position attribute list (shape: nAtoms of 
    # Vectors)
    positions = structure.data.attributes.get('position')
    # convert to a nAtoms x 3 numpy array
    pos_array = np.array(list(map(lambda x: x.vector, positions.data.values())))
    # calculate and apply the rotation matrix to align the principal axes to the
    # xy-plane; keep the rotation matrix so it can be applied again if need be
    new_pos_array, rot_matrix = svd_align(pos_array)
    # update the position attribute values with vectors from the rotated 
    # position array
    positions.data.foreach_set('vector',new_pos_array.reshape(-1).copy(order='c'))
    # force update the position data so that geometry nodes and the 3d viewport
    # are updated.
    structure.data.update()
    
    """

def rotate_vertices(vertices, rot_array, rot_type = 'rotation_matrix'):
    """
    rotate the given vertices by the rotation_matrix
    remove CoG, rotate, add CoG back.

    :parameter vertices: 2d array, assumed to be Nx3 associated with the atomic
        positions of a single structure. 
    :parameter rot_array: numpy array, varying shapes based on the rotation type
    :parameter rot_type: string, accepted values: 'rotation_matrix', 
        'euler_angles', and 'quaternion'
    :return new_vertices: 2d array, same shape as input vertices but now rotated
    """

    import numpy as np

    CoG = np.mean(vertices,axis=0)
    vertices = vertices - CoG
    
    # check if a rotation matrix is given
    if rot_type.lower() == 'rotation_matrix' and rot_array.shape() == (3,3):
        # matrix multiply the vertices by the rot_array transpose
        vertices = vertices @ rot_array.T
        vertices = vertices + CoG
    #elif rot_type.lower() == 'quaternion' and rot_array.shape() == (4,):
    #    # do the thing here
    #elif rot_type.lower() == 'euler_angles' and rot_array.shape() == (4,):
    #    # do the thing here

    return vertices


def zoom_to_fit(camera_object,visual_objects_list,scaling = 1.05):
    """
    given a camera object and a list of visible objects, zoom the camera in/out
    so that white-space is minimized (as determined by the scaling value) 
    without altering the camera's directionality

    :parameter camera_object: bpy_types.Object, assumed to be associated with a
        scene's active camera but doesn't have to be. 
    :parameter visual_objects_list: list of bpy_types.Object, the objects to be
        zoomed in/out on. 
    :parameter scaling: float, magnitude that the camera zoom will be scaled by;
        basically a fudge factor so the zoom isn't toooo zealous.

    NOTE: the bpy.ops.view3d.camera_to_view_selected() returns {'CANCELLED'}
          when a protein structure is drawn using Spheres but works with
          Cartoon. I cannot explain why. To avoid this issue, only apply this 
          function when protein(s) are drawn using Cartoon reps. 
    """
    # grab the original camera vector
    cam_loc_norm = cam.location.normalized()

    # set the important objects' select value to True
    for obj in visual_objects_list:
        obj.select_set(True)

    # run a zealous zoom/crop to focus on the selected objects
    # NOTE: this seems to cause a slight change in camera vector
    bpy.ops.view3d.camera_to_view_selected()
    # deselect all the objects; for bookkeeping sake
    bpy.ops.object.select_all(action='DESELECT')

    # calculate the zealous zoom magnitude
    new_cam_mag = cam.location.magnitude

    # to back away from the zealous zoom and correct any potential changes to
    # the original camera vector, scale the cam_loc_norm by scaling*new_cam_mag
    cam.location = scaling * new_cam_mag * cam_loc_norm

    return


