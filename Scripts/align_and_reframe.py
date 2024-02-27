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
    # transpose Vt to do the right matrix mult
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

    # NOTE: there's two paths to follow from here. Alter the underlying 
    #       attribute values (so the actual atom positions) OR apply a rotation
    #       on the blender object (without changing the atom position values)

    # PATH ONE: update the position attribute values with vectors from the 
    # rotated position array
    positions.data.foreach_set('vector',new_pos_array.reshape(-1).copy(order='c'))
    # force update the position data so that geometry nodes and the 3d viewport
    # are updated.
    structure.data.update()

    # PATH TWO: convert the Vt 3x3 rotation matrix to a quaternion and apply on
    # the blender object's quaternion rotation values
    mu_Vt = Matrix(Vt)  # note blender python pre-imports mathutils
    
    # if you wanna use quaternions
    structure.rotation_mode = 'QUATERNION'
    structure.rotation_quaternion = mu_Vt.to_quaternion()

    # if you wanna use euler angles
    structure.object.rotation_mode = 'XYZ'
    structure.rotation_euler = mu_Vt.to_euler()

    # NOTE: I THINK THIS SECOND PATH IS MORE APPROPRIATE FOR BLENDER WORKFLOWS
    
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

    NOTE: AVOID USING THIS WHEN POSSIBLE, I THINK. SHOULD AVOID ALTERING THE 
          UNDERLYING ATOM ATTRIBUTES AND RATHER APPLY ROTATIONS USING THE 
          OBJECT TRANSFORM OPTIONS
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


def zoom_to_fit(camera_object,visual_objects_list,scaling = 1.05, maintain_vector = True):
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
    :parameter maintain_vector: boolean, set to True if you wish the camera to
        keep the original viewport vector; False if shifts are acceptable. 

    NOTE: the bpy.ops.view3d.camera_to_view_selected() returns {'CANCELLED'}
          when a protein structure is drawn using Spheres but works with
          Cartoon. I cannot explain why. To avoid this issue, only apply this 
          function when protein(s) are drawn using Cartoon reps. 
    """
    if maintain_vector:
        # grab the original camera vector
        cam_loc_norm = camera_object.location.normalized()

    # set the important objects' select value to True
    for obj in visual_objects_list:
        obj.select_set(True)

    # run a zealous zoom/crop to focus on the selected objects
    # NOTE: this seems to cause a slight change in camera vector
    bpy.ops.view3d.camera_to_view_selected()
    # deselect all the objects; for bookkeeping sake
    bpy.ops.object.select_all(action='DESELECT')

    if maintain_vector:
        # calculate the zealous zoom magnitude
        new_cam_mag = camera_object.location.magnitude

        # back away from the zealous zoom and correct any potential changes to
        # the original camera vector, scale the cam_loc_norm by 
        # scaling*new_cam_mag
        camera_object.location = scaling * new_cam_mag * cam_loc_norm

    else:
        # back away from the zealous zoom without worrying about recreating the
        # original view vector
        camera_object.location *= scaling

    return


