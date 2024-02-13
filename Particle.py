from utils import Vector


class Particle:
    """
    This class represents a particle.
    """
    def __init__(self, data, p, h):
        """
        Initialize the 'Particle' class, given 'data' of type
        'ParticleData' for that particle type, the momentum
        four-vector 'p', and the helicity 'h'.
        """
        from math import sqrt
        self.data = data
        self.p = +p
        if self.p[0] < 0:
            self.p[0] = sqrt(sum([pj**2 for pj in p[1:]]) + data.mass**2)
        self.h = float(h)
    def w(self):
        """
        Return the Dirac spinor for this particle.
        """
        from math import sqrt
        if self.data.spin != 2: return None
        # Check if particle or anti-particle.
        s = -1 if self.data.pid < 0 else 1
        p = sqrt(sum([pj**2 for pj in self.p[1:]]))
        # Handle if |p| == p[3].
        if p + self.p[3] == 0:
            xi = 1.0
            if s*self.h == 1: kappa = [0, 1]
            elif s*self.h == -1: kappa = [-1, 0]
            else: kappa = [0, 0]
        # Handle otherwise.
        else:
            xi = 1.0/sqrt(2.0*p*(p + self.p[3]))
            if s*self.h == 1:
                kappa = [p + self.p[3], self.p[2]*1.0j + self.p[1]]
            elif s*self.h == -1:
                kappa = [self.p[2]*1.0j - self.p[1], p + self.p[3]]
            else:
                kappa = [0, 0]
        hp = xi*sqrt(self.p[0] + self.h*p)
        hm = xi*sqrt(self.p[0] - self.h*p)
        # Return the anti-particle spinor.
        if s == -1:
            return Vector(-self.h*kappa[0]*hp, -self.h*kappa[1]*hp,
                           self.h*kappa[0]*hm,  self.h*kappa[1]*hm)
        # Return the particle spinor.
        else:
            return Vector(kappa[0]*hm, kappa[1]*hm, kappa[0]*hp, kappa[1]*hp)
    def wbar(self):
        """
        Return the bar Dirac spinor for this particle.
        """
        w = ~self.w()
        w[0], w[1], w[2], w[3] = w[2], w[3], w[0], w[1]
        return w













