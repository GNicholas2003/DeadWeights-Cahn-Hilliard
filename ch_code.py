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
def double_well(c): #concentration function, with double well
    return (c**2-1)**2/4
def double_well_der(c): #derivative of concentration function
    return (c**2-1)*c
def asym_double_well(c):
     return (c**2-1)**2/4+c#*alpha
def asym_double_well_der(c): #derivative of concentration function
    return (c**2-1)*c+500
def huggins(c):
    return c*np.log(c)+(1-c)*np.log(1-c) + 2.5 * c*(1-c)
def huggins_der(c):
    return np.log(c / (1-c)) + 2.5 * (1 - 2*c)
def triple_well(c):
    return (c+1)**2 * c**2 * (c-1)**2
def triple_well_der(c):
    return 2 * c * (c**2 - 1) * (3*c**2 - 1)
def non_convex(c):
    return (c**2 - 1)**2 / 4 - c
def non_convex_der(c):
    return c**3 - c - 500
def periodic(c):
    return 1-np.cos(c)
def periodic_der(c):
    return np.sin(c)

c=np.random.uniform(low=-1, high=1, size=(N,N)) #initialize random concentrations with same size as grid
const_c=np.full((N,N),-0.5)#constant_concentration
node_c=np.zeros((N,N))
node_c[10,10]=0.5
conc_field[0,:,:]=c #assign concentration to grid
#HUGGINS MODE
#conc_field[0, :, :] = (conc_field[0, :, :] + 1) / 2

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
    mu=double_well_der(prev_c)-kappa*Lap(prev_c,grd_spc)#calculate mu(c)
    d_conc_d_t=mobility_coeff*Lap(mu,grd_spc)#calculate dc/dt
    conc_field[i,:,:]=prev_c+(d_conc_d_t*stp_sz)

    #plot every 500 steps in a graph
    if (i-1)%500==0:
                im= plt.imshow(
                    conc_field[i,:,:],# * 2 - 1, last step, to convert to -1 to 1 range for Huggins model
                    cmap='RdBu',            # Red-White-Blue colormap (Red = -1, White = 0, Blue = +1)
                    origin='lower',         # Places (0, 0) at the bottom-left corner
                    #vmin=-1.0,              # Fixes the minimum colorbar limit
                    #vmax=1.0,               # Fixes the maximum colorbar limit
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
                #plt.savefig(f'Conc_field_step{i-1}-double_well.png')



# #Task 2: Nicholas 
# #use formulas to find m(t) and F(t) and plot them versus time

# #Task 3: Nikola

# #Task 4: Moritz+Luca
