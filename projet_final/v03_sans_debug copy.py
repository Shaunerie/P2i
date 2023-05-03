import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.patches as patches
from matplotlib.patches import Rectangle
from matplotlib import animation
import time
from numba import njit


"""
L'algorithme utilisé provient du pdf du professeur 'Mecanique_des_Materiaux_Granulaires_16_17.pdf'

Hypothèses persos:
- Les grains sont des sphères de même masse et de même rayon
- On est dans le vide
"""

def trajectoire(POSITION, nb_grains, paroiGauche, paroiDroite):
    """
    Affiche la trajectoire des grains dans un graphe matplotlib.

    Paramètres
    ==========
    POSITION : np.array, tableau des positions
    hauteur_silo : float, hauteur du silo
    largeur_silo : float, largeur du silo
    nb_grains : int, nombre de grains

    Retour
    ======
    rien
    """
    print("Affichage de la trajectoire...")

    fig, ax = plt.subplots()
    ax.set_aspect('equal')

    #dessin du silo dans le tableau graphique matplot
    X1 = np.linspace(-0.5, -0.02, 100)
    X2 = np.linspace(0.02, 0.5, 100)
    plt.plot(X1, paroiGauche(X1), color='black')
    plt.plot(X2, paroiDroite(X2), color='black')
    
    for grain in range(nb_grains):
        ax.plot(POSITION[:, grain, 0], POSITION[:, grain, 1], label="grain {}".format(grain))

    plt.grid()
    plt.legend()
    plt.show()

def grain_anime(POSITION, nb_grains, rayon, paroiGauche, paroiDroite):
    """
    Fait une animation de la chute des grains dans le silo.

    Paramètres
    ==========
    POSITION : np.array, tableau des positions
    hauteur_silo : float, hauteur du silo
    largeur_silo : float, largeur du silo
    nb_grains : int, nombre de grains
    rayon : float, rayon des grains

    Retour
    ======
    rien    
    """
    print("Animation en cours...")

    fig, ax = plt.subplots()
    ax.set_aspect('equal')
    
    #dessin du silo dans le tableau graphique matplot
    X1 = np.linspace(-0.5, -0.02, 100)
    X2 = np.linspace(0.02, 0.5, 100)
    plt.plot(X1, paroiGauche(X1), color='black')
    plt.plot(X2, paroiDroite(X2), color='black')
    
    #dessin des grains dans le tableau graphique matplot
    grains = []
    texts = []
    for grain in range(nb_grains):
        grains.append(ax.add_patch(patches.Circle((POSITION[0, grain, 0], POSITION[0, grain, 1]), radius=rayon, fill=True, color='black')))
        texts.append(ax.text(POSITION[0, grain, 0], POSITION[0, grain, 1], str(grain), ha='center', va='center', fontsize=8, color='white'))
    
    time_text = ax.text(0.05, 0.99, '', transform=ax.transAxes, verticalalignment='top', fontsize=12)
    def animate(i):
        #Affiche l'indice du temps en haut a gauche de l'écran
        for grain in range(nb_grains):
            grains[grain].center = (POSITION[i, grain, 0], POSITION[i, grain, 1])
            texts[grain].set_position((POSITION[i, grain, 0], POSITION[i, grain, 1]))
        return grains + texts + [time_text]
    
    ani = animation.FuncAnimation(fig, animate, frames=int(POSITION.shape[0]/1e0), interval=1, blit=True)
    plt.show()



def calcul_allongement_normal(vitesse_i, vitesse_j, position_i, position_j, rayon_i, rayon_j):
    """
    Calcul de l'allongement/distance normal à partir de l'équation
    Paramètres
    ==========
    vitesse_i : np.array, vitesse du grain i
    vitesse_j : np.array, vitesse du grain j
    position_i : np.array, position du grain i
    position_j : np.array, position du grain j
    rayon_i : float, rayon du grain i
    rayon_j : float, rayon du grain j

    Retour
    ======
    allongement_normal : float, allongement normal entre les grains i et j
    """
    return np.linalg.norm(position_i - position_j) - (rayon_i + rayon_j)

def calcul_allongement_tangentiel(i, j, indice_temps, pas_de_temps, VITESSE, POSITION):
    """
    Calcul de l'allongement tangentielle à partir de l'équation

    Paramètres
    ==========
    i : int, indice du grain i
    j : int, indice du grain j
    indice_temps : int, indice du temps actuel
    pas_de_temps : float, pas de temps
    VITESSE : np.array, tableau des vitesses
    POSITION : np.array, tableau des positions

    Retour
    ======
    allongement_tangentiel : float, allongement tangentiel entre les grains i et j
    """
    #print("Contact!")
    #(4.8): d/dt(delta_t) = (vitesse_i - vitesse_j)*tangente 
    # C'est le produit vectorielle de la différence de vitesse sur la tangente du contact entre les deux grains i et j
    #C'est cette formule (4.8) qu'il faut intégrer temporellement en utilisant les données des vitesses(tableaux numpy)
    #On calcul ensuite l'allongement tangentielle via la formule (4.8)
    #Pour ca on cherche le moment d'impact:
    #On cherche le moment d'impact en cherchant le moment ou la distance normal est nulle
    k = 0 # car on sait que c'est déjà le cas
    while calcul_allongement_normal(vitesse_i=VITESSE[indice_temps-k, i],
                          vitesse_j=VITESSE[indice_temps-k, j],
                          position_i=POSITION[indice_temps-k, i],
                          position_j=POSITION[indice_temps-k, j],
                          rayon_i=rayon,
                          rayon_j=rayon) < 0 and k < indice_temps:
        k += 1
    #On a donc trouvé le moment d'impact avec moment = indice_temps-k, NB: si k = 0 alors impact a lieu à indice_temps
    impact = indice_temps - k
    #On peut calculer le tableau des vitesses tangentielles:
    produit_scalaire_array = np.zeros(indice_temps+1) #tableau 1D car projection sur la tangente, 
    #indice_temps+1, parce que plus lisible pour les indices,
    #autrement dit on aura que des zeros jusqu'a l'indice indice_temps-impact   
    for tps_actuel in range(impact, indice_temps+1):
        vecteur_normal = (POSITION[tps_actuel,i] - POSITION[tps_actuel,j])/np.linalg.norm(POSITION[tps_actuel,i] - POSITION[tps_actuel,j])
        #On veut connaître le sens du vecteur tangent qui doit être dans la même direction que la vitesse absolue
        vecteur_tangent = np.array([-vecteur_normal[1], vecteur_normal[0]])
        if np.dot(vecteur_tangent, VITESSE[tps_actuel,i]) < 0:
            vecteur_tangent = -vecteur_tangent
        vitesse_relative = VITESSE[tps_actuel,i] - VITESSE[tps_actuel,j]
        produit_scalaire = np.dot(vitesse_relative, vecteur_tangent)
        produit_scalaire_array[tps_actuel] = produit_scalaire
        
    #On doit maintenant intégrer la valeur de la vitesse tangentielle sur le temps
    allongement_tangentiel = 0
    for i in range(len(produit_scalaire_array)):
        allongement_tangentiel += produit_scalaire_array[i]
    allongement_tangentiel *= pas_de_temps
    return allongement_tangentiel

def application_efforts_distance(masse):
    """
    Application des efforts à distance (par exemple la pesanteur).

    Paramètres
    ==========
    masse : float, masse du grain

    Retour
    ======
    forces : np.array, forces appliquée au grain
    """
    return np.array([0, -masse*9.81])
    

def voisinage(i, POSITION, GRILLE, c, indice_temps):
    """
    Détermine la liste de voisinage du grain i

    Paramètres
    ==========
    i : int, indice du grain i
    POSITION : np.array, tableau des positions
    GRILLE : dict, grille de voisinage
    c : float, pas d'espace
    indice_temps : int, indice du temps

    Retour
    ======
    grains_voisins : np.array, tableau des indices des grains en contact avec le grain i
    """
    x = int(POSITION[indice_temps,i,0]/c)
    y = int(POSITION[indice_temps,i,1]/c)
    voisinage = []
    for voisin in GRILLE[x, y]:
        if voisin != i:
            voisinage.append(voisin)

    for j in range(-1, 2):
        for k in range(-1, 2):
            if not(j == 0 and k == 0):
                case_voisine_x, case_voisine_y = int(POSITION[indice_temps,i,0]/c)+j, int(POSITION[indice_temps,i,1]/c)+k
                if (case_voisine_x, case_voisine_y) in GRILLE.keys():
                    voisinage += GRILLE[case_voisine_x, case_voisine_y]

    return voisinage










def algoprincipal(POSITION, VITESSE, ACCELERATION, VITESSE_DEMI_PAS, nb_grains, raideur_normale, nb_temps, pas_de_temps, rayon, masse, raideur_tangentielle, indice_temps, temps, hauteur_silo_haut, hauteur_silo_bas, largeur_silo_droite, largeur_silo_gauche, raideur_mur, paroiDroite, paroiGauche, coefficient_de_frottement):
    """
    Algorithme principal
    
    Paramètres
    ==========
    POSITION : np.array, tableau des positions
    VITESSE : np.array, tableau des vitesses
    ACCELERATION : np.array, tableau des accélérations
    VITESSE_DEMI_PAS : np.array, tableau des vitesses à demi pas
    nb_grains : int, nombre de grains
    nb_temps : int, nombre de pas de temps
    pas_de_temps : float, pas de temps
    rayon : float, rayon des grainsX1 = np.linspace(-0.5, -0.1, 100)
    X2 = np.linspace(0.1, 0.5, 100)
    masse : float, masse des grains
    raideur_tangentielle : float, raideur tangentielle
    indice_temps : int, indice du temps
    temps : float, temps
    hauteur_silo : float, hauteur du silo
    largeur_silo : float, largeur du silo
    raideur_mur : float, raideur des murs
    paroiDroite : lambda, fonction de la paroi droite
    paroiGauche : lambda, fonction de la paroi gauche

    Retour
    ======
    rien
    """


    start_time = time.time()

    while indice_temps < nb_temps-1:
        #Actualisation du temps
        temps += pas_de_temps
        indice_temps += 1
        if indice_temps % 1000 == 0:
            print("indice_temps ", indice_temps, "/", nb_temps)
            print("temps:", temps)
        
        # Mise a jour des tableaux de position et vitesse à temps k
        for grain in range(nb_grains):
            #Calcul de la nouvelle position du grain
            POSITION[indice_temps, grain, :] = POSITION[indice_temps-1, grain, :] + VITESSE_DEMI_PAS[indice_temps-1, grain, :]*pas_de_temps
    
            
            #Calcul de la vitesse du grain
            VITESSE[indice_temps, grain, :] = VITESSE_DEMI_PAS[indice_temps-1, grain, :] + ACCELERATION[indice_temps-1, grain, :]*pas_de_temps/2

            #On définit une grille pour discrétiser l'espace selon le pas d'espace c, a chaque case on met la liste des grains qui sont dans cette case
            c = 3*rayon #taille de la case
            GRILLE = {(i,j):[] for i in range(int(largeur_silo_gauche/c)-1, int(largeur_silo_droite/c)+2) for j in range(int(hauteur_silo_bas/c)-1, int(hauteur_silo_haut/c)+2)}


            #On associe à chaque case de la grille les grains qui sont dans cette case
            for i in range(nb_grains):
                try:
                    pos_case = (int(POSITION[indice_temps,i,0]/c), int(POSITION[indice_temps,i,1]/c))
                    GRILLE[pos_case[0], pos_case[1]].append(i)
                except:
                    print("problème avec le grain", i)
                    print("position:", POSITION[indice_temps,i,0], POSITION[indice_temps,i,1])
                    print("pos_case:", pos_case)
                    print("indice_temps:", indice_temps)
                    trajectoire(POSITION, nb_grains, paroiDroite=paroiDroite, paroiGauche=paroiGauche)
                    grain_anime(POSITION, nb_grains, rayon, paroiDroite=paroiDroite, paroiGauche=paroiGauche)
                    exit()
            condition_arret = POSITION[indice_temps, grain, 1] < 0.0
            if condition_arret:
                POSITION[indice_temps, grain, 0] = -0.05*grain
                POSITION[indice_temps, grain, 1] = -0.1
                VITESSE[indice_temps, grain, 1] = 0.0
                VITESSE[indice_temps, grain, 0] = 0.0
                ACCELERATION[indice_temps, grain, 1] = 0.0
                ACCELERATION[indice_temps, grain, 0] = 0.0
                VITESSE_DEMI_PAS[indice_temps, grain, 1] = 0.0
                VITESSE_DEMI_PAS[indice_temps, grain, 0] = 0.0

                            
        # Calcul des efforts de contact pour mise à jour des vitesses à temps k+1/2 et accélérations à temps k
        for grain1 in range(nb_grains):
            #Force à distance = gravité
            force_resultante = np.array([0.0, 0.0])
            force_resultante += application_efforts_distance(masse)

            #Rencontre avec une paroi du silo ?
            vecteur_directeur_paroi_gauche = np.array([1.0, -1/0.34])
            vecteur_orthogonal_paroi_gauche = np.array([1/0.34, 1.0]) / np.sqrt(1 + (1/0.34)**2)

            vecteur_directeur_paroi_droite = np.array([1.0, 1/0.34])
            vecteur_orthogonal_paroi_droite = np.array([-1/0.34, 1.0]) / np.sqrt(1 + (1/0.34)**2)

            A = 1/0.34
            B = 1
            C = -0.5
            distance_a_la_gauche = abs(A * POSITION[indice_temps, grain1, 0] + B * POSITION[indice_temps, grain1, 1] + C) / np.sqrt(A**2 + B**2)

            A = -1/0.34
            B = 1
            C = -0.5
            distance_a_la_droite = abs(A * POSITION[indice_temps, grain1, 0] + B * POSITION[indice_temps, grain1, 1] + C) / np.sqrt(A**2 + B**2)

            penetration_gauche = distance_a_la_gauche - rayon
            penetration_droite = distance_a_la_droite - rayon

            if penetration_gauche < 0:
                force_resultante += -raideur_mur * penetration_gauche * vecteur_orthogonal_paroi_gauche

            elif penetration_droite < 0:
                force_resultante += -raideur_mur * penetration_droite * vecteur_orthogonal_paroi_droite

            voisins = voisinage(i=grain1, POSITION=POSITION, c=c, GRILLE=GRILLE, indice_temps=indice_temps)
            for grain2 in voisins:
                if grain1 != grain2:
                    #On définit la force de contact entre les deux grains:
                    force_contact = 0

                    allongement_normal = calcul_allongement_normal(vitesse_i=VITESSE[indice_temps, grain1,:],
                                                                   vitesse_j=VITESSE[indice_temps, grain2, :],
                                                                   position_i=POSITION[indice_temps, grain1, :],
                                                                   position_j=POSITION[indice_temps, grain2, :],
                                                                   rayon_i=rayon,
                                                                   rayon_j=rayon)
                    #Effort normal
                    if allongement_normal < 0:
                        vecteur_normal = (POSITION[indice_temps, grain1, :] - POSITION[indice_temps, grain2, :])/np.linalg.norm(POSITION[indice_temps, grain1, :] - POSITION[indice_temps, grain2, :])
                        force_contact += -raideur_normale * allongement_normal * vecteur_normal

                        #Effort tangentiel
                        allongement_tangentiel = calcul_allongement_tangentiel(grain1, grain2, indice_temps, pas_de_temps, VITESSE, POSITION)
                        vecteur_normal = (POSITION[indice_temps, grain1, :] - POSITION[indice_temps, grain2, :])/np.linalg.norm(POSITION[indice_temps, grain1, :] - POSITION[indice_temps, grain2, :])
                        vecteur_tangentiel = np.array([-vecteur_normal[1], vecteur_normal[0]])
                        if np.dot(vecteur_tangentiel, VITESSE[indice_temps,grain1]) < 0:
                            vecteur_tangentiel = -vecteur_tangentiel
                        force_contact += -raideur_tangentielle * allongement_tangentiel * vecteur_tangentiel
                
                    # 17. Mise à jour de la résultante des forces sur grain1
                    force_resultante += force_contact
            
            frotemment = -coefficient_de_frottement * VITESSE[indice_temps, grain1, :]
            force_resultante += frotemment

            # Calcul de l'accélération du grain à partir de l'équation (4.1)
            ACCELERATION[indice_temps][grain1] = force_resultante / masse
        
            # Calcul de la vitesse de demi-pas à k+1/2 à partir de l'équation (4.19)
            VITESSE_DEMI_PAS[indice_temps][grain1] = VITESSE_DEMI_PAS[indice_temps-1][grain1] + ACCELERATION[indice_temps][grain1] * pas_de_temps / 2


    # Fin de la boucle principale
    print("Fin de la simulation")
    print("Temps de calcul: ", time.time() - start_time, "secondes")

    #Affichage:
    trajectoire(POSITION, nb_grains, paroiGauche, paroiDroite)
    grain_anime(POSITION, nb_grains, rayon, paroiGauche, paroiDroite)




if __name__ == "__main__":
    #Définition grain
    nb_grains = 10
    masse = 1e-3 #kg    
    rayon = 1e-2 #m
    raideur_normale = 100 #N/m
    raideur_tangentielle = 50 #N/m
    coefficient_de_frottement = 0.0005 #N/m
    #Pour le roulement trop compliqué, on utilise leur rotation

    #Définition du silo
    hauteur_silo_bas = -0.5  #m
    hauteur_silo_haut = 2.5 #m
    largeur_silo_gauche = -0.5 #m
    largeur_silo_droite = 0.5 #m
    raideur_mur = 100 #N/m
    paroiGauche = lambda x : -1/0.34*x+0.5
    paroiDroite = lambda x : 1/0.34*x+0.5

    #Définition du temps
    temps = 0
    indice_temps = 0
    pas_de_temps = 1e-3 #s
    duree_simulation = 5
    nb_temps = int(duree_simulation/pas_de_temps)

    #ON PLACE LE REPERE EN BAS A GAUCHE (0,0) DU SILO COIN BAS GAUCHE Y VERS LE HAUT.
    #TABLEAUX NUMPY
    POSITION = np.zeros((nb_temps, nb_grains, 2))   
    for grain in range(nb_grains):
        POSITION[0, grain, 0] = rayon*grain + rayon
        POSITION[0, grain, 1] = 1.5
    print(POSITION[0])
    VITESSE = np.zeros((nb_temps, nb_grains, 2))
    VITESSE[0,:,:] = 0 # pas défini au début, on commece à 1 pour la vitesse et à 0 pour la vitessse de demi pas
    VITESSE_DEMI_PAS = np.zeros((nb_temps, nb_grains, 2))
    VITESSE_DEMI_PAS[0] = np.random.uniform(low=-0.005, high=0.005, size=(nb_grains, 2)) #RAPPEL: Le silo fait un metre par un metre...
    ACCELERATION = np.zeros((nb_temps, nb_grains, 2))
    ACCELERATION[0,:,:] = 0



    start = int(input("Start ou pas? entrez 1 ou 0"))
    if start:
        algoprincipal(POSITION=POSITION, VITESSE=VITESSE, VITESSE_DEMI_PAS=VITESSE_DEMI_PAS, ACCELERATION=ACCELERATION, nb_grains=nb_grains, rayon=rayon, masse=masse, raideur_normale=raideur_normale, raideur_tangentielle=raideur_tangentielle, coefficient_de_frottement=coefficient_de_frottement, paroiGauche=paroiGauche, paroiDroite=paroiDroite, raideur_mur=raideur_mur, pas_de_temps=pas_de_temps, nb_temps=nb_temps, indice_temps=indice_temps, hauteur_silo_bas=hauteur_silo_bas, hauteur_silo_haut=hauteur_silo_haut, largeur_silo_gauche=largeur_silo_gauche, largeur_silo_droite=largeur_silo_droite, temps=temps)