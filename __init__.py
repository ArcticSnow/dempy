# -*- coding: utf-8 -*-
#
# Copyright © 2014 Simon Filhol
# Licensed under the terms of the MIT License. 

"""
DEM toolbox:
    - create artificial
    - process raster and DEM
    - analyze
    - save DEM as raster files
"""
import numpy as np
from dempy import fftdem
from dempy import linearDetrend
from dempy import MatrixGenerator
from dempy import rasterTool
from dempy import plotting
from dempy import morphometry
from dempy import smooth2d
from dempy import Qgis_tool