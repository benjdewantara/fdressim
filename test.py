import numpy as np
from scipy import linalg

sle = np.array([[-0.6878, 0.1268, 0, 0, 0],
                [0.1268, -0.8146, 0.1268, 0, 0],
                [0, 0.1268, -0.8146, 0.1268, 0],
                [0, 0, 0.1268, -0.8146, 0.1268],
                [0, 0, 0, 0.1268, -0.6878]])

known = np.array([-3365.98, -3365.98, -3365.98, -3215.98, -3365.98])

nextPres = linalg.solve(sle, known)