import numpy as np

def verlet(X, t):
    if t == 0:
        return
    pass

if __name__ == '__main__':
    #X(t + dt) = 2X(t) - X(t - dt) - a(t)*(dt)^2

    K = 10 ; uma= 1.66*10**-27 ; n = 55.55 ; a = 156*10**-12 ; X0 = 2.3*a ; V_trac = 10**-3 ; N = 10 ; 
    N_iter = 10**5 ; dt = 10**-15; dx = V_trac*dt ; l0=(N+1)*X0 ; 


    #Initiatialisation:
    X = np.zeros((N+1, N_iter+1))


    #Réccurence
    t = 1
    for atome in range(N):
        X[atome] = verlet(X, t)
