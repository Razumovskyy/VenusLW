import numpy as np
import matplotlib.pyplot as plt

height, cool1, cool2 = np.loadtxt(fname='_COOLING', skiprows=1, unpack=True)
fig, ax = plt.subplots()

ax.plot(1013.25*cool1, height, color='r')
ax.plot(1013.25*cool2, height, color='b')
plt.xlabel('Cooling rate [K/day]')
plt.ylabel('Height [km]')
plt.grid(visible=True, which='major', color='g')
plt.xlim(-0.2, 0.2)
plt.ylim(0, 100)
plt.show()