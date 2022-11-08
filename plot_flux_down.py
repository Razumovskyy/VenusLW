from turtle import color
import numpy as np
import matplotlib.pyplot as plt

nu, fdown1, fdown2 = np.loadtxt(fname='_Fdown_', skiprows=1, unpack=True)

fig, ax = plt.subplots()

ax.plot(nu, fdown1, color='r')
ax.plot(nu, fdown2, color='b')
plt.show()