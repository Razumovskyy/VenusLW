from turtle import color
import numpy as np
import matplotlib.pyplot as plt

nu, fup1, fup2 = np.loadtxt(fname='_Fup_', skiprows=1, unpack=True)

fig, ax = plt.subplots()

ax.plot(nu, fup1, color='r')
ax.plot(nu, fup2, color='b')
plt.show()