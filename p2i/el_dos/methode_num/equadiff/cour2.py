from numpy import pi, cos, sin
import numpy as np
import matplotlib.pyplot as plt

#2.1 Schéma d'Euler explicite
# Yi2 = Yi1 + dt * F
#2)
l=1; g=10; alpha=1; T=10; N=100; dt= T/N
X = np.linspace(0, T, N+1)
Y1 = np.zeros(N+1)
Y2 = np.zeros(N+1)

theta = alpha*pi/180 #deg --> rad
theta_p = 0
Y1[0], Y2[0] = theta, theta_p

for t in range(1, N+1):
    theta = theta + dt*theta_p
    theta_p = theta_p -dt*(g/l)*sin(theta)
    Y1[0], Y2[0] = theta, theta_p

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14,4)) #le 1,2 veut dire une ligne 2 colonne soit deux graphes

ax1.set_xlabel(r'temps séparé de dt=1/0')
ax1.set_ylabel(r'\theta(t)')
ax1.set_title(r'Evolution de \theta(t) au cours du temps.')

ax2.set_xlabel(r'temps séparé de dt=1/0')
ax2.set_ylabel(r'\theta^.(t)')
ax2.set_title(r'Evolution de \theta^.(t) au cours du temps.')

ax1.plot(X, Y1)
ax2.plot(X, Y2)


