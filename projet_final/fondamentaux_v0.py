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

def calcul_allongement(vitesse_i, vitesse_j, position_i, position_j):
    """
    Calcul de l'allongement des ressorts tangentiels à partir de l'équation (4.8)
    in: vitesse_i, vitesse_j, position_i, position_j (tableaux numpy)
    out: allongement (tableau numpy)
    """
    #Raideur tangentielle = allongement
    #(4.8): d\delta/dt = (vitesse_i - vitesse_j)*tangente 
    # C'est le produit vectorielle de la différence de vitesse sur la tangente du contact entre les deux grains i et j
    #C'est cette formule (4.8) qu'il faut intégrer temporellement en utilisant les données des vitesses(tableaux numpy)
    tangente = ?
    derivee = ?
    integration = ?
    return integration

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

# 1. Définition de l'état initial pour chaque particule composant l'échantillon
# (À adapter en fonction de votre situation)
etat_initial = []

# 2. Définition de la valeur du pas de temps
pas_de_temps = 0.01

# 3. Définition de la durée souhaitée de la simulation
duree_simulation = 10

# 4. Initialisation du temps
temps = 0

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

