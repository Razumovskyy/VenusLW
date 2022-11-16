from turtle import color
import numpy as np
import matplotlib.pyplot as plt

height, fup1, fup2 = np.loadtxt(fname='_Fup_', skiprows=1, unpack=True)

fig, ax = plt.subplots()

ax.plot(fup1, height, color='r')
ax.plot(fup2, height, color='b')
plt.grid(which='major', axis='both', color='g')
plt.xlim(0, 50)
plt.ylim(0, 140)
plt.xlabel('Flux, W/m^2')
plt.ylabel('Height, km')
plt.show()