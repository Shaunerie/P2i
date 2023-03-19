import numpy as np
import matplotlib.pyplot as plt
from numpy import cos, tan, sin, arctan, pi


#paramètre géo.
h = 0.1 
alpha = 15*pi/180

#mouvement d'entrée:
#psi point = omega = une cste
psi = np.linspace(0, 2*pi, 360) # l'incrément corrsepond à un degré

#les inconnus sont theta et x

#on veut plot la fonction analytique:
solution_analytique_1 = -h*cos(psi)*tan(alpha)
solution_analytique_2 = -tan(psi)*cos(alpha)
"""
#y'a des probleme avec la discontinuite donc calcul maitrisé:
sauter = []
for i in range(1, psi.shape[0]):
    if abs((-tan(psi[i])*cos(alpha)) -(-tan(psi[i-1])*cos(alpha))) > 2:
        sauter.append(i)
psi2 = np.array([0])
solution_analytique_2 = np.array([0])
for valeur in range(psi.shape[0]):
    if valeur not in sauter:
        psi2 = np.append(psi2, np.array([psi[valeur]]))
        solution_analytique_2 = np.append(solution_analytique_2, np.array([(-tan(psi[valeur])*cos(alpha))]))
fig, ax = plt.subplots(figsize=(8,6))
"""
#solution1
plt.xlabel(r'$ \psi $ (rd)')
plt.ylabel('x (m)')
plt.title("Position de l'axe de sortie(bleu).")
plt.plot(psi, solution_analytique_1, label='x analytique', color = 'blue')
plt.legend()
plt.show()

#solution2
plt.xlabel(r'$ \psi $ (rd)')
plt.ylabel('tan(theta) (SI)')
plt.title("Angle entre le boitier et l'arbre d'entrée(rouge).")
plt.plot(psi, solution_analytique_2,'-', label='theta analytique', color = 'red')
plt.legend()
plt.show()

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#de manière numérique:
#equation 1: x*cos(alpha) + h*cos(psi)*sin(alpha) = 0
#equation 2: x*sin(theta)*sin(alpha) − h*(cos(psi)*sin(theta)*cos(alpha) + sin(psi)*cos(theta)) = 0
psi_plot = np.linspace(0, 2*pi, 36) #ici on passe l'incrément à 10degré
x_plot = np.zeros(36)
theta_plot = np.zeros(36)
for i in range(36):
    #variables à présenter à la fin:

    #Newton
    omega = 1 #parce que je connais pas l'influence de psi point
    ecart = 1
    #les inconnues sont x, theta
    x = 0
    theta = 0
    z = np.array([x, theta])
    #le paramètre variable est psi:
    psi = psi_plot[i]
    t = 0
    while ecart > 1e-8 and t < 10000:
        f1 = x*cos(alpha) + h*cos(psi)*sin(alpha)
        f2 = x*tan(theta)*sin(alpha) - h*(cos(psi)*tan(theta)*cos(alpha) + sin(psi))
        f = np.array([f1, f2])
        #on pose que tan(theta est l'inconnu don on dérive en csq:)
        J1 = [cos(alpha), 0]
        J2 = [-sin(theta)*sin(alpha), x*sin(alpha) - h*(cos(psi)*cos(alpha))]
        J = np.array([J1, J2])

        dz = -np.linalg.inv(J)@f
        z = z + dz
        #on actualise l'écart et les valeurs:
        ecart = np.linalg.norm(dz)
        x, theta = z[0], z[1]
        t += 1
    #a la fin de la boucle on a trouver x et theta assez proche alors:
    theta_plot[i] = z[1]
    x_plot[i] = z[0]

#on peut ploter x et theta:
psi = np.linspace(0, 2*pi, 360)
#x
fig_phi, ax = plt.subplots(figsize=(8, 6))
ax.plot(psi, solution_analytique_1, label = 'x analytique')
ax.plot(psi_plot, x_plot, 'r*', label = 'x numérique')
ax.set_xlabel(r'$ \psi $  (rd)')
ax.set_ylabel(r'x (m)')
ax.set_title("Position de l'axe de sortie.")
ax.legend()
#theta mais attention il y a des discontinuités:
fig_phi, ax = plt.subplots(figsize=(8, 6))
ax.plot(psi, solution_analytique_2, label = r'$\theta$ analytique')
ax.plot(psi_plot, theta_plot, 'r*', label = r'$\theta$ numérique')
ax.set_xlabel(r'$ \psi $  (rd)')
ax.set_ylabel(r'$\theta$ (rd)')
ax.set_title("Orientation boitier/arbre entrée.")
ax.legend()

