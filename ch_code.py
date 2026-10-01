import numpy as np
import matplotlib.pyplot as plt
#Task 1: Samuel Simulating concentration profile over a NxN grid

N=64 #grid_size
t=10000 #total_passed_time/steps
conc_field=np.zeros((t,N,N)) #initialize concentration field
mobility_coeff=1 #constant for C-H-Problem
kappa=1 #constant for C-H-Problem
stp_sz=0.01 #individual time step, important for Discrete Laplacian
grd_spc=1 #spacing between grids, important for Discrete Laplacian
set_cbar=True #legend, should only be plotted once
def conc_fnc(c): #concentration function, with double well
    return (c**2-1)**2/4
def conc_fnc_der(c): #derivative of concentration function
    return (c**2-1)*c

c=np.random.uniform(low=-1, high=1, size=(N,N)) #initialize random concentrations with same size as grid
const_c=np.full((N,N),-0.5)#constant_concentration
node_c=np.zeros((N,N))
node_c[10,10]=0.5
conc_field[0,:,:]=node_c #assign concentration to grid

#Test command:
# print(conc_field[0,:,:])
# print(conc_field[0,1,1])
# print(conc_field[0,11%N,11%N])
#Discretized Laplacian operator function:
def Lap(f,dx):
    #Consider periodic boundary conditions
    return (
        np.roll(f, 1, axis=0) + # move elements along specified axis and position, e.g. axis 0 is x axis, 1 means moved by one position
        np.roll(f, -1, axis=0) + 
        np.roll(f, 1, axis=1) + 
        np.roll(f, -1, axis=1) - 
        4 * f
    ) / (dx**2)


for i in range(1,t):
    
    #calculate partial derivatives and add to former concentration profile
    prev_c=conc_field[i-1,:,:] #use previous concentration field as temporal variable
    mu=conc_fnc_der(prev_c)-kappa*Lap(prev_c,grd_spc)#calculate mu(c)
    d_conc_d_t=mobility_coeff*Lap(mu,grd_spc)#calculate dc/dt
    conc_field[i,:,:]=prev_c+(d_conc_d_t*stp_sz)

    #plot every 500 steps in a graph
    if (i-1)%500==0:
                im= plt.imshow(
                    conc_field[i,:,:],
                    cmap='RdBu',            # Red-White-Blue colormap (Red = -1, White = 0, Blue = +1)
                    origin='lower',         # Places (0, 0) at the bottom-left corner
                    vmin=-1.0,              # Fixes the minimum colorbar limit
                    vmax=1.0,               # Fixes the maximum colorbar limit
                    interpolation='bicubic' # Smooths pixel transitions (useful for continuous fields)
                    
                )
                if set_cbar: #set legend
                    cbar=plt.colorbar(im)
                    cbar.set_label('Concentration Field ',fontsize=11)
                    cbar.set_ticks([-1.0,0.0,1.0])
                    cbar.set_ticklabels(['-1.0 (Phase A)', '0.0 (Interface)', '+1.0 (Phase B)'])
                    set_cbar=False
                plt.title(f'Step {i-1}')
                plt.tight_layout()
                plt.savefig(f'Node_Conc_field_step{i-1}')



# #Task 2: Nicholas 
# #use formulas to find m(t) and F(t) and plot them versus time
# #lalalalalalalala lalalalal



# #Task 3: Nikola

# #Task 4: Moritz+Luca