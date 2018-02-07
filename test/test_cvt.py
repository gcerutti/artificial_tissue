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

def test_cvt_sphere():
    import vplants.artificial_tissue
    pkg_path = vplants.artificial_tissue.__path__[0]+"/../../.."
    from vplants.artificial_tissue.surface_to_tissue.mesh2cvt import mesh_to_cvt_image

    n_points = 10
    img = mesh_to_cvt_image(input=pkg_path+"/share/data/meshes/icosphere.ply",nbcells=n_points,voxelsize=(1,1,1))

    assert len(np.unique(img))==n_points+1
