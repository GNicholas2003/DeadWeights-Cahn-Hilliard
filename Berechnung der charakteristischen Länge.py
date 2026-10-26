import numpy as np
import matplotlib.pyplot as plt
from ch_code import *

def DFT(t):
    matrix = conc_field[t,:,:]
    DFT = np.fft.fft2(matrix)
    return DFT

def DFT_Zeigen(t):
    plt.figure(t*10)
    plt.imshow(conc_field[t,:,:])
    plt.show()

DFT_Zeigen(90)
DFT_Zeigen(100)

def S(t):
    s = 0
    k = 0
    center_y, center_x = (DFT(t).shape[0] - 1) / 2, (DFT(t).shape[1] - 1) / 2
    y, x = np.ogrid[:DFT(t).shape[0], :DFT(t).shape[1]] # Nummeriert die Matrix durch
    distance_matrix = np.sqrt((y - center_y)**2 + (x - center_x)**2) #Berechnet den Abstand zum Zentrum

    unique_distances = np.unique(np.round(distance_matrix, 10))
    a = []
    for dist in unique_distances:
        
        group = DFT(t)[np.isclose(distance_matrix, dist)] #Alle Elemente mit ähnlichem Abstand
        group = np.append(group, dist.real)
        average = np.mean(group[:-2])
        structur_factor = average.real**2+average.imag**2

        a.append(structur_factor)

        if structur_factor > s:
            s = structur_factor
            k = dist  
    return s, k

def L(t):
    return 2*np.pi/S(t)[1]




A = []
B = [n for n in range(0, 1000)]
for i in range(0, 1000):
    A.append(L(i))


plt.figure(2)
plt.plot(B, A)
plt.xscale('log')
#plt.yscale('log')
#plt.savefig("L_vs_t_random")
plt.show()


A = L(0)
a = np.log(L(50)/A)/np.log(50)
print('--------------------------------------------------')
print(A, a)



