import numpy as np
import matplotlib.pyplot as plt

nu, cool1, cool2 = np.loadtxt(fname='_COOLING', skiprows=1, unpack=True)

fig, ax = plt.subplots()

ax.plot(nu, cool1, color='r')
ax.plot(nu, cool2, color='b')
plt.show()