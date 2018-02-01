# -*- coding: utf-8 -*-
# -*- python -*-
#
#       Artificial Tissue
#
#       Copyright 2017-2018 INRIA
#
#       File author(s): Guillaume Cerutti <guillaume.cerutti@inria.fr>
#
#       File contributor(s): Guillaume Cerutti <guillaume.cerutti@inria.fr>
#
#       Distributed under the Cecill-C License.
#       See accompanying file LICENSE.txt or copy at
#           http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.html
#
###############################################################################

import numpy as np


def test_1_layer_sphere():
    from vplants.artificial_tissue.layered_sphere import layered_sphere_tissue_image
    
    size = 60
    n_points = 10

    img = layered_sphere_tissue_image(size=size, n_points=n_points)

    assert img.shape == (size,size,size)
    assert len(np.unique(img)) == n_points+2


def test_2_layer_sphere():
    from vplants.artificial_tissue.layered_sphere import layered_sphere_tissue_image
    
    size = 60
    n_points = 10

    img = layered_sphere_tissue_image(size=size, n_points=n_points, n_layers=2)

    assert img.shape == (size,size,size)