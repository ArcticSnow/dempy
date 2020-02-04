import pyximport; pyximport.install()
import diamondSquare as ds

from time import time

t0 = time()
a = ds.diamondSquare(8000, 8000, 1, 0.2)
t1 = time()
print('Processing time: ' + str(t1 - t0))


from matplotlib import pyplot as plt 
plt.imshow(a)
plt.show()


