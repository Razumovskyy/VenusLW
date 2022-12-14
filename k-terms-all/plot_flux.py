import numpy as np
import matplotlib.pyplot as plt

height, fup1, fup2 = np.loadtxt(fname='CO2_Fup_3_', skiprows=1, unpack=True)

fig, ax = plt.subplots()

ax.plot(fup1, height, color='r')
ax.plot(fup2, height, color='b')
plt.grid(visible=True, which='major', color='g')
plt.xlim(0, 100)
plt.ylim(40, 120)
plt.xlabel('Flux [W/m^2]')
plt.ylabel('Height [km]')
plt.show()