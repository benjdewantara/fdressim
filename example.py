# I need to style this dashboard file,
# so that one can just type "from core import *"
# and simply build a model then run a simulation

from core import *

# Firstly, we can specify the dimension of the cartesian grid
nz, ny, nx = 1, 5, 5
dims = (nz, ny, nx)
g = Grid(dims)

# Then, we can specify the whole reservoir dimension
Lz, Ly, Lx = nz*75, ny*1000, nx*1000
resDimension = (Lz, Ly, Lx)

# Build the fluid and rock model
# Implement this soon!
f = Fluid("Some fluid")
r = Rock("Some rock")

# We contain all these informations in a Reservoir object
res = Reservoir(grid=g, fluid=f, rock=r, resDim=resDimension)

# By default, the moment we declare a Node object, a no-flow Neumann
# condition has already been imposed if the Node is a boundary Node.
# But we can specify another condition with another value as follows
bc = BoundaryCondition()
res.addBoundaryCondition(bc, x='before')

# Set the initial pressure array
res.setInitPressure(1000)



