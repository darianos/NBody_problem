import numpy as np
import scipy as scp
import matplotlib.pylab as plt
from mpl_toolkits.mplot3d import Axes3D


G = 39.398 #En U.A, masas solares y años
N = 365
h = 1.0/N
#########
n = 3  ## Numero de particulas
m = np.genfromtxt("Masa_Prueba.txt")
XPas = np.genfromtxt("CI_Prueba.txt", delimiter=",", usecols=(0,1,2))  ##Posicion global del sistema de n particulas en el instante inicial
VPas = np.genfromtxt("CI_Prueba.txt", delimiter=",", usecols=(3,4,5))	##Las primeras 3(n) columnas serán la posicion y las otras 3(n) serán la velocidad
AceleracionINI = np.zeros((n,3))
#########
XPres = np.zeros((n,3))  ##Posicion global del sistema de n particulas en un instante posterior
VPres = np.zeros((n,3))
Aceleracion = np.zeros((n,3))
#########
for i in range(n):
	B = 0
	for k in range(n):
		if i!=k:
			AceleracionINI[i,:] = G*m[k]*(XPas[k,:] - XPas[i,:])*np.linalg.norm(XPas[k,:] - XPas[i,:])**(-3)
			B = B + AceleracionINI[i,:]
		VPres[i,:] = VPas[i,:] + 0.5*h*B
##########
plt.figure()
ax = plt.axes(projection="3d")
ax.set_xlabel("Posicion en x (U.A)")
ax.set_ylabel("Posicion en y (U.A)")
ax.set_zlabel("Posicion en z (U.A)")

for T in range(N): #For de tiempo
	for i in range(n): #For de selección de particulas
		A = 0
		for k in range(n): #For sobre todas las particulas
			if i!=k: #Selección de particulas, diferentes a la i-ésima
				Aceleracion[i, :] = G*m[k]*(XPas[k, :] - XPas[i, :])*(np.linalg.norm(XPas[k, :] - XPas[i, :]))**(-3)
				A = A + Aceleracion[i, :]  # Aceleración sobre el cuerpo i-ésimo
			XPres[i, :] = XPas[i, :] + h*VPas[i, :]
			VPres[i, :] = VPas[i, :] + h*A
			ax.scatter(XPres[0, 0], XPres[0, 1], XPres[0, 2], marker="_",c="red")  # el segundo indice hace referncia a la parte coordenada, x=0,y=1,z=2
			ax.scatter(XPres[1, 0], XPres[1, 1], XPres[1, 2], marker="_", c="green")
			ax.scatter(XPres[2, 0], XPres[2, 1], XPres[2, 2], marker="_", c="blue")


	for f in range(n):
		XPas[f, :] = XPres[f, :]
		VPas[f, :] = VPres[f, :]

plt.grid()
plt.show()

