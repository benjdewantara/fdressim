# I need to style this control file,
# so that one can just type "from core import *"
# and simply build a model then run a simulation

from core import *
from differentiator import runSimulation


# Firstly, we can specify the dimension of the cartesian grid
nz, ny, nx = 1, 1, 5
dims = (nz, ny, nx)
g = Grid(dims)

# Then, we can specify the whole reservoir dimension
Lz, Ly, Lx = 75, 1000, 5000
resDimension = (Lz, Ly, Lx)

# Build the fluid and rock model
# See to this later!
f = Fluid(refRho=62.428, refPres=14.7, compress=3.5*1e-6, mu=10)
r = Rock(refPoro=0.18, refPres=14.7, compress=0, perm=0.015)

# rho is in lbm/ft^3
# refPres is in psi
# compress is in psi^-1
# mu is in cP
# perm is in D (Darcy)



# We contain all information in a Reservoir object
res = Reservoir(grid=g, fluid=f, rock=r, resDim=resDimension)

# By default, the moment we declare a Node object, a no-flow Neumann
# condition has already been imposed if the Node is a boundary Node.

# Set the initial pressure array
res.setInitPressure(6000)


# Set a source/sink in coordinate (0, 0, 3)
res.grid.nodes[np.ravel_multi_index((0, 0, 3), res.grid.dims)].setSrc(-150)


# Finally, run the simulation!
runSimulation(res, dt=15, nTime=24)

#runSimulation2(res, dt=15, nTime=30)




