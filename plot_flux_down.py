from turtle import color
import numpy as np
import matplotlib.pyplot as plt

height, fdown1, fdown2 = np.loadtxt(fname='_Fdown_', skiprows=1, unpack=True)

fig, ax = plt.subplots()

ax.plot(fdown1, height, color='r')
ax.plot(fdown2, height, color='b')
plt.grid(which='major', color='g')
plt.xlim(0, 500)
plt.ylim(0, 140)
plt.show()