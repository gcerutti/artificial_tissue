# -*- coding: utf-8 -*-
# -*- python -*-
#
#       Artificial Tissue
#
#       Copyright 2017-2018 INRIA
#
#       File author(s): Hadrien Oliveri <hadrien.oliveri@inria.fr>
#
#       File contributor(s): Hadrien Oliveri <hadrien.oliveri@inria.fr>
#
#       Distributed under the Cecill-C License.
#       See accompanying file LICENSE.txt or copy at
#           http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.html
#
###############################################################################

import argparse
from itertools import product
import os
from utils import *
from openalea.image.serial.all import imread, imsave
from openalea.mesh.property_topomesh_io import read_ply_property_topomesh

def mesh_to_cvt_image(input, output='.', voxelsize=(.5, .5, .5), max_step=1e9, nbcells=100, method="lloyd",
                      verbose=False, debug=False, save=False):
    """

    Args:
        input: Path to PLY mesh file.
        output: Path to directory where intermediary and final images are stored. Default is current directory.
        voxelsize: Voxel size. Default is (.5,.5,.5).
        max_step: Maximum step for the centroidal Voronoi tessellation. Default is 1e9
        nbcells: Number of cells. Default is 100. # XXX What happens if the number of cells is higher than the number of voxels encompassed by the mesh ?
        method: Method to compute the centroidal Voronoi tessellation: 'lloyd' or 'mcqueen'. Default is 'lloyd'. #FIXME McQueen does not seem to converge. Bug?
        verbose: Verbose mode. Default is False
        debug: Debug mode. Default is False
        save: Save images. Default is False

    Returns: CVT image

    """

    def save_img(path, img):
        logging.info("Saving image to " + path)
        imsave(path, voronoi_img)

    assert os.path.exists(input), "Input file: " + input + " does not exist."
    assert method in ["lloyd", "mcqueen"], "Wrong method."

    if output is not None and save and not os.path.exists(output):
        os.makedirs(output)

    logging.getLogger().setLevel(logging.INFO if verbose else logging.DEBUG if debug else logging.ERROR)
    logging.info("### Generating cellularized image from mesh. ###")

    # Loading mesh
    logging.info("Loading mesh")
    mesh = read_ply_property_topomesh(input)

    # Binarizing image
    logging.info("Rasterizing mesh to binary image.")
    mask = topomesh_to_binary_image(mesh=mesh, voxelsize=voxelsize, verbose=verbose, debug=debug)

    if save and output is not None: save_img(os.path.join(output, "bin.inr"), mask)

    points = filter(lambda p: mask[p], list(product(range(mask.shape[0]), range(mask.shape[1]), range(mask.shape[2]))))

    # Generating random seeds
    logging.info("Choosing initial seeds.")
    density = None  # TODO
    seeds, labels, seed_img = random_seeds(img=mask, nb_seeds=nbcells, density=density, points=points,
                                           verbose=verbose, debug=debug, replace=False)
    if save and output is not None: save_img(os.path.join(output, "seeds.inr"), seed_img)

    # Initial Voronoi partioning
    logging.info("Initializing Voronoi diagram.")
    voronoi_img = voronoi(seeds, labels, mask=mask, points=points, verbose=verbose, debug=debug)
    if save and output is not None: save_img(os.path.join(output, "voronoi.inr"), voronoi_img)


    # Computing a centroidal Voronoi tessellation.
    logging.info("Computing a centroidal Voronoi tessellation.")
    centroid_img, cvt_img = cvt(mask=mask, seeds=seeds, labels=labels, res_path=os.path.join(output, "cvt.inr"),
                            steps=max_step, points=points, method=method, voronoi_img=voronoi_img, save=save, verbose=verbose,
                                debug=debug)

    return cvt_img

def main():
    """

    Returns:

    """

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', help='Path to input mesh', required=True)
    parser.add_argument('-o', '--output', help='Path to output files directory', default="./output")
    parser.add_argument('-n', '--nbcells', help='Number of cells', default=100, type=int)
    parser.add_argument('--step', help='Maximal number of steps for CVT', default=1e9, type=int)
    parser.add_argument('-v', '--verbose', default=False, action='store_true', help='verbose')
    parser.add_argument('-d', '--debug', default=False, action='store_true', help='debug')
    parser.add_argument('-s', '--save', default=False, action='store_true', help='save output image to file')
    parser.add_argument('--voxelsize', help='Voxel size', default=[.025, .025, .025], nargs=3, type=float)
    parser.add_argument('-m', '--method', help='Method for CVT [\'lloyd\', \'mcqueen\']', default='lloyd')

    args = parser.parse_args()

    mesh_to_cvt_image(input=args.input, output=args.output, method=args.method, verbose=args.verbose, debug=args.debug,
        save=args.save, voxelsize=args.voxelsize, nbcells=args.nbcells, max_step=args.step)

if __name__ == "__main__":
    main()

