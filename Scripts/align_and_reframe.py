"""
Run this script when a structure is initially loaded into a blender visualization
"""

# PREAMBLE
import MDAnalysis

# to handle blender objects... this is gonna be the painful part
import bpy

def rotate_pcs(vertices):
    """
    remove CoM translation, rotate the vertices to align principal components of
    vertices along the xyz dimensions, add back CoM translation, return vertices

    use singular value decomposition to calculate the rotation matrix for 
    aligning the vertices to the xyz dimensions

    :parameter vertices: 2d array, assumed to be Nx3 associated with the atomic
        positions of a single structure. 
    :return new_vertices: 2d array, same shape as input vertices but now rotated.
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
    
    return new_vertices


def resize_camera(camera_object,visual_objects):
    """
    https://blender.stackexchange.com/questions/78048/full-script-for-fitting-the-camera-view-to-the-only-object-in-a-scene
    """


