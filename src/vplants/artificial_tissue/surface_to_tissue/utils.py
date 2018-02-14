# -*- coding: utf-8 -*-
# -*- python -*-
#
#       Artificial Tissue
#
#       Copyright 2017-2018 INRIA
#
#       File author(s): Hadrien Oliveri <hadrien.oliveri@inria.fr>
#
#       File contributor(s): Hadrien Oliveri <hadrien.oliveri@inria.fr>, Guillaume Cerutti <guillaume.cerutti@inria.fr>
#
#       Distributed under the Cecill-C License.
#       See accompanying file LICENSE.txt or copy at
#           http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.html
#
###############################################################################

import numpy as np
import logging
from itertools import product
from openalea.image.spatial_image import SpatialImage, empty_image_like
from scipy.spatial import cKDTree


def topomesh_to_binary_image(mesh, voxelsize=(.5, .5, .5), verbose=False, debug=False):
    """

    :param polydata:
    :param voxelsize:
    :param bounds:
    :param verbose:
    :param debug:
    :return:
    """

    from vtk import VTK_MAJOR_VERSION as vtk_version
    from vtk import vtkImageData, vtkPolyDataToImageStencil, vtkImageStencil, VTK_UNSIGNED_CHAR
    from vtk.util.numpy_support import vtk_to_numpy
    from openalea.cellcomplex.property_topomesh.triangular_mesh import topomesh_to_triangular_mesh

    logging.getLogger().setLevel(logging.INFO if verbose else logging.DEBUG if debug else logging.ERROR)

    mesh, matching = topomesh_to_triangular_mesh(mesh, degree=3, coef=1.0, property_name=None, property_degree=None)
    polydata = mesh._repr_vtk_()
    bounds = polydata.GetBounds()

    # Initializing image
    logging.info("Initializing image")
    white_image = vtkImageData()
    white_image.SetSpacing(voxelsize)
    dim = [int(np.ceil((bounds[i * 2 + 1] - bounds[i * 2]) / voxelsize[i])) for i in xrange(3)]

    # FIXME: does not work if dimensions are different
    dim = [np.amax(dim)] * 3

    logging.info("Image dimension: " + str(dim))
    white_image.SetDimensions(dim)
    white_image.SetExtent(-1, dim[0] - 1, -1, dim[1] - 1, -1, dim[2] - 1)
    origin = [bounds[i * 2] for i in xrange(3)]
    # origin = [.5 *(bounds[i * 2]+bounds[i * 2 + 1]) for i in xrange(3)]

    white_image.SetOrigin(origin)
    logging.info("Image origin: " + str(white_image.GetOrigin()))
    if vtk_version < 6:
        white_image.SetScalarTypeToUnsignedChar()
        white_image.SetNumberOfScalarComponents(1)
        white_image.AllocateScalars()
    else:
        white_image.AllocateScalars(VTK_UNSIGNED_CHAR, 1)

    count = white_image.GetNumberOfPoints()
    for i in xrange(count):
        white_image.GetPointData().GetScalars().SetTuple1(i, 255)

    pol2stenc = vtkPolyDataToImageStencil()
    if vtk_version < 6:
        pol2stenc.SetInput(polydata)
    else:
        pol2stenc.SetInputData(polydata)
    pol2stenc.SetOutputOrigin(origin)
    pol2stenc.SetOutputSpacing(voxelsize)
    pol2stenc.SetOutputWholeExtent(white_image.GetExtent())

    imgstenc = vtkImageStencil()
    if vtk_version < 6:
        imgstenc.SetInput(white_image)
        imgstenc.SetStencil(pol2stenc.GetOutput())
    else:
        imgstenc.SetInputData(white_image)
        imgstenc.SetStencilConnection(pol2stenc.GetOutputPort())
    imgstenc.ReverseStencilOff()
    imgstenc.SetBackgroundValue(0)
    imgstenc.Update()

    im = imgstenc.GetOutput()
    rows, cols, depth = im.GetDimensions()
    sc = im.GetPointData().GetScalars()
    a = vtk_to_numpy(sc)
    a = a.reshape(rows, cols, depth)

    assert a.shape == im.GetDimensions()

    return SpatialImage(a, vdim=1, dtype=np.uint16)


def topomesh_binary_image(topomesh, size=None, voxelsize=(1,1,1), origin=[0,0,0], verbose=False, debug=False):

    from vtk import VTK_MAJOR_VERSION as vtk_version
    from vtk import vtkImageData, vtkPolyDataToImageStencil, vtkImageStencil, VTK_UNSIGNED_CHAR
    from vtk.util.numpy_support import vtk_to_numpy
    from openalea.cellcomplex.property_topomesh.triangular_mesh import topomesh_to_triangular_mesh

    logging.getLogger().setLevel(logging.INFO if verbose else logging.DEBUG if debug else logging.ERROR)

    mesh, matching = topomesh_to_triangular_mesh(topomesh, degree=3, coef=1.0, mesh_center=[0,0,0])
    polydata = mesh._repr_vtk_()
    bounds = polydata.GetBounds()

    if origin is None:
        origin = [bounds[i*2] for i in xrange(3)]

    if size is None:
        size = [(bounds[i*2 + 1] - origin[i]) for i in xrange(3)]

    shape = [int(np.ceil(s/v)) for s,v in zip(size,voxelsize)]
    
    print "--> Bounds : ",bounds
    print "--> Size : ",size
    print "--> Shape : ",shape
    print "--> Origin : ",origin

    white_image = vtkImageData()
    white_image.SetSpacing(voxelsize)
    white_image.SetDimensions(shape)
    
    # white_image.SetExtent(origin[0], origin[0]+shape[0]*voxelsize[0], origin[1], origin[1]+shape[1]*voxelsize[1], origin[2], origin[2]+shape[2]*voxelsize[2])
    # white_image.SetExtent(0, shape[0]*voxelsize[0], 0, shape[1]*voxelsize[1], 0, shape[2]*voxelsize[2])
    white_image.SetExtent(0, shape[0]-1, 0, shape[1]-1, 0, shape[2]-1)
    # white_image.SetOrigin([bounds[i*2] for i in xrange(3)])
    # white_image.SetOrigin(origin)
    # white_image.SetOrigin(np.array(origin) - 20.*np.array([bounds[i*2] for i in xrange(3)]))
    if vtk_version < 6:
        white_image.SetScalarTypeToUnsignedChar()
        white_image.SetNumberOfScalarComponents(1)
        white_image.AllocateScalars()
    else:
        white_image.AllocateScalars(VTK_UNSIGNED_CHAR, 1)

    print "White image : ",white_image.GetDimensions()

    count = white_image.GetNumberOfPoints()
    for i in xrange(count):
        white_image.GetPointData().GetScalars().SetTuple1(i, 255)

    pol2stenc = vtkPolyDataToImageStencil()
    if vtk_version < 6:
        pol2stenc.SetInput(polydata)
    else:
        pol2stenc.SetInputData(polydata)
    # pol2stenc.SetOutputOrigin(np.array(origin) - 2.*np.array([bounds[i*2] for i in xrange(3)]))
    pol2stenc.SetOutputOrigin(origin)
    # pol2stenc.SetOutputOrigin([bounds[i*2] for i in xrange(3)])
    pol2stenc.SetOutputSpacing(voxelsize)
    # pol2stenc.SetOutputWholeExtent(bounds)
    pol2stenc.SetOutputWholeExtent(0, shape[0]-1, 0, shape[1]-1, 0, shape[2]-1)

    imgstenc = vtkImageStencil()
    if vtk_version < 6:
        imgstenc.SetInput(white_image)
        imgstenc.SetStencil(pol2stenc.GetOutput())
    else:   
        imgstenc.SetInputData(white_image)
        imgstenc.SetStencilConnection(pol2stenc.GetOutputPort())
    imgstenc.ReverseStencilOff()
    imgstenc.SetBackgroundValue(0)
    imgstenc.Update()

    binary_image = imgstenc.GetOutput()
    print "Mask image : ",binary_image.GetDimensions()
    rows, cols, depth = binary_image.GetDimensions()
    sc = binary_image.GetPointData().GetScalars()
    binary_img = vtk_to_numpy(sc)
    binary_img = np.transpose(binary_img.reshape(depth,cols,rows),(2,1,0))

    return SpatialImage(binary_img.astype(np.uint16), voxelsize=voxelsize, origin=origin)


def random_seeds(img, nb_seeds=100, density=None, points=None, replace=False, verbose=False, debug=False):
    """

    :param img:
    :param nb_seeds:
    :param density:
    :param points:
    :param replace:
    :param verbose:
    :param debug:
    :return:
    """

    assert nb_seeds > 0, "Illegal number of seeds: " + str(nb_seeds)
    logging.getLogger().setLevel(logging.INFO if verbose else logging.DEBUG if debug else logging.ERROR)

    if points is None:
        assert img is not None
        points = list(product(range(img.shape[0]), range(img.shape[1]), range(img.shape[2])))
        points = filter(lambda p: img[p], points)

    volume = len(points)
    logging.info("Generating " + str(nb_seeds) + " random seeds in " + str(volume) + " voxels.")

    res_img = empty_image_like(img)
    logging.info("Uniform density of seeds." if density is None else "Non-uniform density of seeds")

    rand_indices = np.random.choice(a=len(points), size=nb_seeds, replace=replace, p=density)
    seeds = np.array(points)[rand_indices]

    labels = np.uint16(np.arange(nb_seeds) + 1)
    res_img[tuple(np.array(seeds).T)] = labels

    return seeds, labels, res_img


def voronoi(seeds, labels, mask=None, points=None, verbose=False, debug=False):
    """

    :param seeds:
    :param labels:
    :param mask:
    :param points:
    :param verbose:
    :param debug:
    :return:
    """
    logging.getLogger().setLevel(logging.INFO if verbose else logging.DEBUG if debug else logging.ERROR)

    if points is None:
        assert mask is not None and mask.shape == mask.shape
        points = filter(lambda p: mask[p],
                        list(product(xrange(mask.shape[0]), xrange(mask.shape[1]), xrange(mask.shape[2]))))

    # KD tree for nearest point research
    nearest = cKDTree(data=seeds).query(points)[1]
    res_img = empty_image_like(mask)
    res_img[tuple(np.transpose(points))] = labels[nearest]

    return res_img


def voronoi_image(seeds, labels, mask, seed_tensors=None, verbose=False, debug=False):
    """
    """
    
    voxelsize = np.array(mask.voxelsize) if hasattr(mask,'voxelsize') else np.ones(3,float)
    
    if seed_tensors is None:
        seed_tensors = np.array([np.diag(np.ones(3)) for s in seeds])

    points = np.transpose(np.mgrid[0:mask.shape[0],0:mask.shape[1],0:mask.shape[2]],(1,2,3,0))*voxelsize
    seed_distances = np.array([np.power(np.einsum('...ij,...ij->...i', points-seed,np.einsum('...ij,...j->...i',np.linalg.inv(np.power(seed_tensor,2.)),points-seed)),0.5) for seed,seed_tensor in zip(seeds,seed_tensors)]) 
    voronoi_img = ((mask!=0)*np.array(labels)[np.argmin(seed_distances,axis=0)]).astype(np.uint16)

    return voronoi_img



def centroids(img, labels, bg=0, coords=None):
    """

    :param img: input labelled image
    :param labels:
    :param bg: background value
    :param coords: list of point coordinates
    :return:
    """

    if coords is None:
        assert np.ndim(img) == 3
        l0, l1, l2 = np.shape(img)
        coords = filter(lambda (x, y, z): img[x, y, z] != bg, product(xrange(l0), xrange(l1), xrange(l2)))

    labels_coords = {l: [] for l in labels if l != bg}

    for x, y, z in coords:
        labels_coords[img[x, y, z]] += [[x, y, z]]

    centroids = dict(zip(labels, [None] * len(labels)))
    for l in labels:
        pos = labels_coords[l]
        centroids[l] = np.sum(pos, axis=0).astype(float) / len(pos)

    return np.array(centroids.values())


def cvt(mask, seeds, labels, steps=1e3, voronoi_img=None, res_path=None, points=None,
        method="lloyd", intermediary_step=5e6, verbose=False, debug=False, save=False, threshold=1.e-9):
    """
    Centroidal Voronoi tessellation

    :param mask:
    :param seeds:
    :param labels:
    :param steps:
    :param voronoi_img:
    :param res_path:
    :param points:
    :param method:
    :param intermediary_step:
    :return:
    """

    assert method in ['lloyd', 'macqueen'], "Illegal method argument " + method
    assert steps >= 0, "Illegal number of steps"
    from openalea.image.serial.all import imsave

    logging.getLogger().setLevel(logging.INFO if verbose else logging.DEBUG if debug else logging.ERROR)
    logging.info("Centroidal Voronoi tessellation. Method: " + method)

    def residual(s1, s2):
        return np.max(np.linalg.norm((np.asarray(s1) - np.asarray(s2)) * np.array(mask.voxelsize), axis=1))\
               / np.power(len(points) * np.prod(mask.voxelsize), 1/3.)

    seeds = np.double(seeds)

    if method == "lloyd":

        assert np.ndim(mask) == 3
        assert voronoi_img is not None
        prev_seeds = seeds
        for step in xrange(int(steps)):
            seeds = centroids(voronoi_img, labels, coords=points)
            res = residual(seeds, prev_seeds)
            logging.info("CVT step " + str(step) + " --> Residual:" + str(res))
            voronoi_img = voronoi(seeds, labels, mask, points, verbose=verbose, debug=debug)
            prev_seeds = seeds

            # Intermediary saving.
            if res_path is not None and save:
                imsave(res_path, voronoi_img)

            if np.isclose(0, res, atol=threshold):  # TODO specify convergence threshold
                logging.info("Lloyd algorithm converged.")
                break

    elif method == "macqueen":

        # FIXME: does not converge
        raise NotImplementedError()

        points = np.double(points)
        nb_seeds = len(seeds)
        memory = dict(zip(range(nb_seeds), [0] * nb_seeds))
        recompute_voronoi = True
        for step in xrange(len(steps)):
            if not step % (intermediary_step / 5):
                logging.info("CVT step " + str(step))
            random_point = points[np.random.choice(len(points))]
            closest = np.argmin(np.linalg.norm(seeds - random_point, axis=1))
            memory[closest] += 1
            seeds[closest] = (random_point + memory[closest] * np.asarray(seeds[closest])) / (memory[closest] + 1)

            # Intermediary saving.
            if not step % intermediary_step and step:
                voronoi_img = voronoi(seeds, labels, mask, np.uint16(points))
                c = centroids(voronoi_img, labels)
                res = residual(seeds, c)
                logging.info("CVT step " + str(step) + " --> Residual:" + str(res))

                if res_path is not None and save:
                    imsave(res_path, voronoi_img)

                if np.isclose(0, res):  # TODO specify convergence threshold
                    recompute_voronoi = False
                    logging.info("McQueen algorithm converged.")
                    break

        if recompute_voronoi:
            voronoi_img = voronoi(seeds, labels, mask, points)
            if res_path is not None:
                imsave(res_path, voronoi_img)

    return seeds, voronoi_img
