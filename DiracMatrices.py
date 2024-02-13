from utils import FourVector, Matrix


class DiracMatrices(FourVector):
    """
    This class provides the Dirac matrices. Note that this class
    inherits from the 'FourVector' class. This is because the Dirac
    matrices also transform under the Minkowski metric, just like
    standard four-vectors.
    """
    def __init__(self, v0 = None, v1 = None, v2 = None, v3 = None):
        """
        Initialize the Dirac matrices. Ideally this would not be mutable.
        """
        g0 = Matrix([0.0, 0.0, 1.0, 0.0], [0.0, 0.0, 0.0, 1.0],
                    [1.0, 0.0, 0.0, 0.0], [0.0, 1.0, 0.0, 0.0])
        g1 = Matrix([0.0, 0.0, 0.0, 1.0], [0.0, 0.0, 1.0, 0.0],
                    [0.0, -1.0, 0.0, 0.0], [-1.0, 0.0, 0.0, 0.0])
        g2 = Matrix([0.0, 0.0, 0.0, -1.0j], [0.0, 0.0, 1.0j, 0.0],
                    [0.0, 1.0j, 0.0, 0.0], [-1.0j, 0.0, 0.0, 0.0])
        g3 = Matrix([0.0, 0.0, 1.0, 0.0], [0.0, 0.0, 0.0, -1.0],
                    [-1.0, 0.0, 0.0, 0.0], [0.0, 1.0, 0.0, 0.0])
        FourVector.__init__(self, g0, g1, g2, g3)







