import numpy as np
import scipy as scp
import matplotlib.pylab as plt
from mpl_toolkits.mplot3d import Axes3D


G = 39.398 #En U.A, masas solares y años
N = 1
h = 1.0/N
#########
n = 3  ## Numero de particulas
def RsideV(g,M,Rrel,NormRrel):
	return -1*g*M*Rrel*NormRrel**(-3)
m = np.genfromtxt("Masa_Prueba.txt")
XPas = np.genfromtxt("CI_Prueba.txt", delimiter=",", usecols=(0,1,2))  ##Posicion global del sistema de n particulas en el instante inicial
VPas = np.genfromtxt("CI_Prueba.txt", delimiter=",", usecols=(3,4,5))	##Las primeras 3(n) columnas serán la posicion y las otras 3(n) serán la velocidad
#########
XPres = np.zeros((n,3))  ##Posicion global del sistema de n particulas en un instante posterior
VPres = np.zeros((n,3))

for i in range(n):
	for k in range(n): 
		if i!=k:
			VPres[i,:] = VPas[i,:] + 0.5*h*np.sum(RsideV(G, m[k], XPas[k,:] - XPas[i,:], np.linalg.norm(XPas[i,:] - XPas[k,:])))  # La suma debe ser sobre k
