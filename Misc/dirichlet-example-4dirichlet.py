# I need to style this control file,
# so that one can just type "from core import *"
# and simply build a model then run a simulation

from core import *
from differentiator import runSimulation
from visualizeUsingMRST import generateMFile


# Firstly, we can specify the dimension of the cartesian grid
nz, ny, nx = 1, 51, 51
dims = (nz, ny, nx)
g = Grid(dims)

# Then, we can specify the whole reservoir dimension
Lz, Ly, Lx = 75, 5000, 5000
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



# We contain all these informations in a Reservoir object
res = Reservoir(grid=g, fluid=f, rock=r, resDim=resDimension)

# By default, the moment we declare a Node object, a no-flow Neumann
# condition has already been imposed if the Node is a boundary Node.
# But we can specify another condition with another value as follows
#bc = BoundaryCondition("d", 6000)
#res.addBoundaryCondition(bc, x='before')
#res.addBoundaryCondition(bc, x='after')
#res.addBoundaryCondition(bc, y='before')
#res.addBoundaryCondition(bc, y='after')

# Set the initial pressure array
res.setInitPressure(6000)


# Set a source/sink in coordinate (0, 0, 3)
res.grid.nodes[np.ravel_multi_index((0, 25, 25), res.grid.dims)].setSrc(-150)


# Finally, run the simulation!
runSimulation(res, dt=0.25, nTime=5)



# Temporarily use MRST's plotCellData
#generateMFile('example-4', res.grid.dims, res.resDim, 23)



