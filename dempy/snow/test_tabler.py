import tabler

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

'''
Example script on how to use Tabler's function
'''


line = np.zeros(1000)

for i in range(1, line.__len__()):
        line[i] = line[i - 1] + np.random.uniform(-1, 1)

lineM = pd.Series(line).rolling(window=10).mean()
line = np.array(lineM)

dx= 1.2
xline = np.arange(0, line.__len__() * dx, dx)

tab = tabler.tabler_profile(line,dx=dx)

plt.figure()
plt.plot(xline, line)
plt.plot(tab[0,:],tab[1,:])
plt.show()





























