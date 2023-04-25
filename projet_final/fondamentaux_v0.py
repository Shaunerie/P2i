import numpy as np
"""
L'algorithme utilisé provient du pdf du professeur 'Mecanique_des_Materiaux_Granulaires_16_17.pdf'

Hypothèses persos:
- Les grains sont des sphères de même masse et de même rayon
- On est dans le vide
"""


def calcul_nouvelle_position(position, vitesse, acceleration, pas_de_temps):
    """
    Calcul de la nouvelle position à partir de l'équation (4.18)
    in: position, vitesse, acceleration, pas_de_temps (tableaux numpy)
    out: nouvelle position (tableau numpy)
    """
    return position + vitesse * pas_de_temps

def calcul_vitesse(vitesse, acceleration, pas_de_temps):
    """
    Calcul de la vitesse de demi-pas de temps à partir de l'équation (4.20)
    in: vitesse, acceleration, pas_de_temps (tableaux numpy)
    out: vitesse de demi-pas de temps (tableau numpy)
    """
    return vitesse + acceleration * pas_de_temps/2

def calcul_disance_normal(vitesse_i, vitesse_j, position_i, position_j, rayon_i, rayon_j):
    """
    Calcul de l'allongement/distance normal à partir de l'équation
    in: vitesse_i, vitesse_j, position_i, position_j (tableaux numpy), rayon_i, rayon_j (float)
    out: distance/allongement normal (float)
    """
    return np.linalg.norm(position_i - position_j) - (rayon_i + rayon_j)

def calcul_allongement_tangentiel(vitesse_i, vitesse_j, position_i, position_j, rayon_i, rayon_j):
    """
    Calcul de l'allongement tangentielle à partir de l'équation
    in: vitesse_i, vitesse_j, position_i, position_j (tableaux numpy), rayon_i, rayon_j (float)
    out: allongement tangentielle (float)
    """
    #(4.8): d/dt(delta_t) = (vitesse_i - vitesse_j)*tangente 
    # C'est le produit vectorielle de la différence de vitesse sur la tangente du contact entre les deux grains i et j
    #C'est cette formule (4.8) qu'il faut intégrer temporellement en utilisant les données des vitesses(tableaux numpy)
    #On calcul ensuite l'allongement tangentielle via la formule (4.8)
    vecteur_normal = (position_i - position_j)/np.linalg.norm(position_i - position_j)
    vitesse_relative = vitesse_i - vitesse_j
    vitesse_tangentielle = vitesse_relative - np.dot(vitesse_relative, vecteur_normal)*vecteur_normal
    #On doit maintenant intégrer la valeur de la vitesse tangentielle sur le temps
    #En fait c pas facile parce qu'on doit utiliser le cumsum avec le pas de temps mais a partir 
    #d'un temps précis, celui du contact entre les deux grains... ou ca correspond au repos du ressort
    #tangentielle (allongement_tangetielle = 0)


def application_efforts_distance(masse):
    """
    Application des efforts à distance (par exemple la pesanteur)
    in: masse (float)
    out: effort à distance (tableau numpy)
    """
    # Remplacez cette ligne par l'application des efforts à distance (par exemple la pesanteur)
    return numpy.array([0, -masse*9.81])
    

def calcul_effort_contact(grain1, grain2):
    # Remplacez cette ligne par le calcul de l'effort de contact entre grain1 et grain2
    pass

def mise_a_jour_resultante(grain, force_contact):
    # Remplacez cette ligne par la mise à jour de la résultante des forces sur grain
    pass


if __name__ == "__main__":
    # 0. Définition des paramètres physiques du problème
    masse = 1e-3 #kg
    rayon = 1e-2 #m
    raideur_normale = 1 #N/m
    raideur_tangentielle = 1 #N/m
    coefficient_de_frottement = 0.5 #N/m
    nb_grains = 1000
    # 2. Définition de la valeur du pas de temps
    pas_de_temps = 0.01
    # 3. Définition de la durée souhaitée de la simulation
    duree_simulation = 10
    # 4. Initialisation du temps
    temps = 0
    #Pour le roulement trop compliqué, on utilise leur rotation
    nb_temps = duree_simulation//pas_de_temps
    #Définition du silo
    hauteur_silo = 1 #m
    largeur_silo = 1 #m

    # 1. Définition de l'état initial pour chaque particule composant l'échantillon
        ???
    # 5. Tant que temps < duree_simulation, boucle principale
    while temps < duree_simulation:
        # 6. Avancée dans le temps
        temps += pas_de_temps
        
        # 7. Boucle sur tous les grains
        for grain in etat_initial:
            # 8. Calcul de la nouvelle position du grain
            calcul_nouvelle_position(grain)
            
            # 9. Calcul de la vitesse du grain
            calcul_vitesse(grain)
            
            # 10. Mise à jour de l'allongement des ressorts tangentiels
            calcul_allongement(grain)
        
        # 12. Boucle sur tous les grains
        for grain1 in etat_initial:
            # 13. Initialisation des résultantes d'efforts sur le grain
            force_resultante = np.zeros(3)
            
            # 14. Application des efforts à distance (par exemple la pesanteur)
            application_efforts_distance(grain1)
            
            # 15. Boucle sur tous les grains
            for grain2 in etat_initial:
                if grain1 != grain2:
                    # 16. En cas de contact entre grain1 et grain2, calcul de l'effort de contact
                    force_contact = calcul_effort_contact(grain1, grain2)
                    
                    # 17. Mise à jour de la résultante des forces sur grain1
                    mise_a_jour_resultante(grain1, force_contact)
            
            # 19. Calcul de l'accélération du grain à partir de l'équation (4.1)
            # (Remplacez cette ligne par le calcul de l'accélération à partir de l'équation (4.1))
            
            # 20. Calcul de la vitesse de demi-pas de temps à partir de l'équation (4.19)
            # (Remplacez cette ligne par le calcul de la vitesse de demi-pas à partir de l'équation (4.19))

        # 22. Sauvegarde, Tracé, Monitoring, etc.
        # (Ajoutez ici le code pour sauvegarder, tracer et/ou surveiller les données de la simulation)

    # 23. Fin de la boucle principale

