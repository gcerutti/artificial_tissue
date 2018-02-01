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
import os
from utils import *
from openalea.image.serial.all import imread, imsave
from openalea.mesh.property_topomesh_io import read_ply_property_topomesh


def mesh_to_cvt_image(input, output='.', voxelsize=(.5, .5, .5), max_step=1e9, nbcells=100, method="lloyd",
                      verbose=False, debug=False, save=False):
    """

    :param input:
    :param output:
    :param voxelsize:
    :param method:
    :param verbose:
    :param debug:
    :return:
    """

    assert os.path.exists(input), "Input file: " + input + " does not exist."
    assert method in ["lloyd", "mcqueen"], "Wrong method."

    if output is not None and save and not os.path.exists(output):
        os.makedirs(output)

    logging.getLogger().setLevel(logging.INFO if verbose else logging.DEBUG if debug else logging.ERROR)
    logging.info("Generating image from mesh.")

    # Loading mesh
    logging.info("Loading mesh")
    mesh = read_ply_property_topomesh(input)

    # Binarized image
    logging.info("Rasterizing mesh to binary image.")
    mask = topomesh_to_binary_image(mesh=mesh, voxelsize=voxelsize, verbose=verbose, debug=debug)
    if save and output is not None:
        logging.info("Saving image to " + output)
        bin_path = os.path.join(output, "bin.inr")
        imsave(bin_path, mask)
    from itertools import product
    points = filter(lambda p: mask[p], list(product(range(mask.shape[0]), range(mask.shape[1]), range(mask.shape[2]))))

    # Random seeds
    logging.info("Choosing initial seeds.")
    density = None  # TODO
    seeds, labels, seed_img = random_seeds(img=mask, nb_seeds=nbcells, density=density, points=points,
                                           verbose=verbose, debug=debug, replace=False)
    if save and output is not None:
        logging.info("Saving image to " + output)
        seed_path = os.path.join(output, "seeds.inr")
        imsave(seed_path, seed_img)

    # Initial Voronoi partioning
    logging.info("Initializing Voronoi diagram.")
    voronoi_img = voronoi(seeds, labels, mask=mask, points=points, verbose=verbose, debug=debug)
    if save and output is not None:
        logging.info("Saving image to " + output)
        voronoi_path = os.path.join(output, "voronoi.inr")
        imsave(voronoi_path, voronoi_img)

    # CVT
    logging.info("Computing CENTROIDAL VORONOI TESSELLATION")
    centroid_img, cvt_img = cvt(mask=mask, seeds=seeds, labels=labels, res_path=os.path.join(output, "cvt.inr"),
                            steps=max_step, points=points, method=method, voronoi_img=voronoi_img, save=save, verbose=verbose,
                                debug=debug)

    return cvt_img

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', help='path to input mesh', required=True)
    parser.add_argument('-o', '--output', help='path to output files directory', default="./output")
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

