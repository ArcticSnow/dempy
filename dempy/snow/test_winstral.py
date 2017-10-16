#test_winstral


# library required for running the script in Python
import matplotlib.pyplot as plt
from matplotlib.colors import LightSource
import numpy as np
import winstral

import sys 
sys.path.append('/home/arcticsnow/github/dempy/dempy/decomposition/')
import ds


dem = ds.diamondSquare(100, 100, 100, .8)


ls = LightSource(azdeg=315, altdeg=45)

# extrac cellsize from raster
cellsize = 1        # extract pixel size from raster

# set the wind direction and the maximum search distance
in_wind = -90
dmax = 100 #m

# execute the algorithm
Sx_max = winstral.winstralSX(dem, cellsize, dmax, in_wind)

# Plot the results
plt.figure()
plt.subplot(121)
cmap = plt.cm.gist_earth
plt.imshow(ls.shade(dem,cmap=plt.cm.terrain, vert_exag=.1, blend_mode='soft'))
plt.title('Dem sample')

plt.subplot(122)
plt.imshow(Sx_max)
plt.colorbar()
plt.title('Sx_max, Dmax='+str(dmax))
plt.show()
