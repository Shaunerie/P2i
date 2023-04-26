import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

"""
L'algorithme utilisé provient du pdf du professeur 'Mecanique_des_Materiaux_Granulaires_16_17.pdf'

Hypothèses persos:
- Les grains sont des sphères de même masse et de même rayon
- On est dans le vide
"""
def calcul_distance_normal(vitesse_i, vitesse_j, position_i, position_j, rayon_i, rayon_j):
    """
    Calcul de l'allongement/distance normal à partir de l'équation
    in: vitesse_i, vitesse_j, position_i, position_j (tableaux numpy), rayon_i, rayon_j (float)
    out: distance/allongement normal (float)
    """
    return np.linalg.norm(position_i - position_j) - (rayon_i + rayon_j)

def calcul_allongement_tangentiel(i, j, indice_temps, pas_de_temps, VITESSE, POSITION):
    """
    Calcul de l'allongement tangentielle à partir de l'équation
    in: VITESSE, POSITION (tableaux numpy 3D), rayon_i, rayon_j (float), indice_temps (int), pas_de_temps (float)
    out: allongement tangentielle (float)
    """

    #Si il ne sont pas en contact on ne fait rien:
    if calcul_distance_normal(vitesse_i=VITESSE[indice_temps][i][:],
                              vitesse_j=VITESSE[indice_temps][j][:],
                              position_i=POSITION[indice_temps][i][:],
                              position_j=POSITION[indice_temps][j][:],
                              rayon_i=rayon,
                              rayon_j=rayon) > 0:
        return 0

    #Sinon on continue
    else:
        print("Contact!")
        #(4.8): d/dt(delta_t) = (vitesse_i - vitesse_j)*tangente 
        # C'est le produit vectorielle de la différence de vitesse sur la tangente du contact entre les deux grains i et j
        #C'est cette formule (4.8) qu'il faut intégrer temporellement en utilisant les données des vitesses(tableaux numpy)
        #On calcul ensuite l'allongement tangentielle via la formule (4.8)
        #Pour ca on cherche le moment d'impact:
        #On cherche le moment d'impact en cherchant le moment ou la distance normal est nulle
        j = 2 # car on sait que c'est déjà le cas et que c'était déjà le cas à l'indice_temps-1
        while calcul_distance_normal(vitesse_i=VITESSE[indice_temps-j][i][:],
                              vitesse_j=VITESSE[indice_temps-j][j][:],
                              position_i=POSITION[indice_temps-j][i][:],
                              position_j=POSITION[indice_temps-j][j][:],
                              rayon_i=rayon,
                              rayon_j=rayon) > 0:
            j += 1
        #On a donc trouvé le moment d'impact avec moment = indice_temps-j

        #On peut calculer le tableau des vitesses tangentielles:
        norme_tangentielle_array = np.zeros(j+1) #car on prend en compte le temps actuelle + la différence représenté par j, tableau 1D car projection sur la tangente


        for tps_actuel in range(indice_temps-j, indice_temps+1):
            vecteur_normal = (POSITION[tps_actuel][i][:] - POSITION[tps_actuel][j][:])/np.linalg.norm(POSITION[tps_actuel][i][:] - POSITION[tps_actuel][j][:])
            vitesse_relative = VITESSE[tps_actuel][i][:] - VITESSE[tps_actuel][j][:]
            vitesse_tangentielle = vitesse_relative - np.dot(vitesse_relative, vecteur_normal)*vecteur_normal
            norme_vitesse_tangentielle = np.linalg.norm(vitesse_tangentielle)
            norme_tangentielle_array[tps_actuel] = norme_vitesse_tangentielle
            
        #On doit maintenant intégrer la valeur de la vitesse tangentielle sur le temps
        allongement_tangentiel = 0
        for i in range(len(norme_tangentielle_array)):
            allongement_tangentiel += norme_tangentielle_array[i]
        allongement_tangentiel *= pas_de_temps
        return allongement_tangentiel
    


def application_efforts_distance(masse):
    """
    Application des efforts à distance (par exemple la pesanteur)
    in: masse (float)
    out: effort à distance (tableau numpy)
    """
    return np.array([0, -masse*9.81])
    


def algoprincipal(POSITION, VITESSE, ACCELERATION, VITESSE_DEMI_PAS, ALLONGEMENT_TANGENTIEL, nb_grains, nb_temps, pas_de_temps, rayon, masse, raideur_tangentielle, indice_temps, temps, hauteur_silo, largeur_silo, raideur_mur):
    """
    Algorithme principal
    in: POSITION, VITESSE, ACCELERATION, VITESSE_DEMI_PAS (tableaux numpy), nb_grains, nb_temps, pas_de_temps, rayon, masse, raideur_tangentielle (float)
    out: None
    """

    while indice_temps < nb_temps-1:
        #Actualisation du temps
        temps += pas_de_temps
        indice_temps += 1
        
        # Mise a jour des tableaux de position et vitesse à temps k
        for grain in range(nb_grains):
            #Calcul de la nouvelle position du grain
            POSITION[indice_temps][grain] = POSITION[indice_temps-1][grain] + VITESSE_DEMI_PAS[indice_temps-1][grain] * pas_de_temps
            
            #Calcul de la vitesse du grain
            VITESSE[indice_temps][grain] = VITESSE_DEMI_PAS[indice_temps-1][grain] + ACCELERATION[indice_temps-1][grain] * pas_de_temps/2
        
        # Calcul des efforts de contact pour mise à jour des vitesses à temps k+1/2 et accélérations à temps k
        for grain1 in range(nb_grains):
            #Force à distance = gravité
            force_resultante = application_efforts_distance(grain1)

            #Rencontre avec un mur ?
            condition_droite = POSITION[indice_temps][grain1][0] + rayon > largeur_silo
            condition_gauche = POSITION[indice_temps][grain1][0] - rayon < 0
            condition_haut = POSITION[indice_temps][grain1][1] + rayon > hauteur_silo
            condition_bas = POSITION[indice_temps][grain1][1] - rayon < 0

            if...

            
            for grain2 in range(nb_grains):
                if grain1 != grain2:
                    allongement = calcul_allongement_tangentiel(grain1, grain2, indice_temps, pas_de_temps, VITESSE, POSITION)
                    ALLONGEMENT_TANGENTIEL[grain1][grain2] = calcul_allongement_tangentiel(grain1, grain2, indice_temps, pas_de_temps, VITESSE, POSITION)
                    #En cas de contact entre grain1 et grain2, calcul de l'effort de contact
                    force_contact = raideur_tangentielle * ALLONGEMENT_TANGENTIEL[grain1][grain2]
                    
                    # 17. Mise à jour de la résultante des forces sur grain1
                    force_resultante += force_contact

            # Calcul de l'accélération du grain à partir de l'équation (4.1)
            ACCELERATION[indice_temps][grain1] = force_resultante / masse
            
            # Calcul de la vitesse de demi-pas à k+1/2 à partir de l'équation (4.19)
            VITESSE_DEMI_PAS[indice_temps][grain1] = VITESSE_DEMI_PAS[indice_temps-1][grain1] + ACCELERATION[indice_temps][grain1] * pas_de_temps / 2

        # 22. Sauvegarde, Tracé, Monitoring, etc.
        

    # 23. Fin de la boucle principale
    print("Fin de la simulation")
    # 24. Affichage des résultats
    fig, ax = plt.subplots()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    
    for grain in range(nb_grains):
        ax.plot(POSITION[:, grain, 0], POSITION[:, grain, 1], label="grain {}".format(grain))
    
    #dessin du silo



    plt.grid()
    plt.legend()
    plt.show()





if __name__ == "__main__":
    #Définition grain
    nb_grains = 10
    masse = 1e-3 #kg
    rayon = 1e-2 #m
    raideur_normale = 1000 #N/m
    raideur_tangentielle = 1000 #N/m
    coefficient_de_frottement = 0.5 #N/m
    #Pour le roulement trop compliqué, on utilise leur rotation
    #Définition du silo
    hauteur_silo = 0.5 #m
    largeur_silo = 0.5 #m
    raideur_mur = 1000 #N/m

    temps = 0
    indice_temps = 0
    pas_de_temps = 0.01
    duree_simulation = 10
    nb_temps = int(duree_simulation/pas_de_temps)
    #Définition du silo
    hauteur_silo = 1 #m
    largeur_silo = 1 #m
    raideur_mur = 1000 #N/m

    #ON PLACE LE REPERE EN BAS A GAUCHE (0,0) DU SILO COIN BAS GAUCHE Y VERS LE HAUT.
    #TABLEAUX NUMPY
    POSITION = np.zeros((nb_temps, nb_grains, 2))
    POSITION[0] = np.random.uniform(low=0, high=1, size=(nb_grains, 2))
    VITESSE = np.zeros((nb_temps, nb_grains, 2))
    VITESSE[0,:,:] = 0 # pas défini au début, on commece à 1 pour la vitesse et à 0 pour la vitessse de demi pas
    VITESSE_DEMI_PAS = np.zeros((nb_temps, nb_grains, 2))
    VITESSE_DEMI_PAS[0] = np.random.uniform(low=-1, high=1, size=(nb_grains, 2))
    ACCELERATION = np.zeros((nb_temps, nb_grains, 2))
    ACCELERATION[0,:,:] = np.random.uniform(low=-1, high=1, size=(nb_grains, 2))

    ALLONGEMENT_TANGENTIEL = np.zeros((nb_grains, nb_grains)) #On stocke les allongements tangentiels pour chaque couple de grains en contact

    algoprincipal(POSITION, VITESSE, ACCELERATION, VITESSE_DEMI_PAS, ALLONGEMENT_TANGENTIEL, nb_grains, nb_temps, pas_de_temps, rayon, masse, raideur_tangentielle, indice_temps, temps, hauteur_silo, largeur_silo, raideur_mur)