from DiracMatrices import DiracMatrices

class Annihilate:
    """
    This class defines the cross-section function needed to calculate
    the integrated cross-section of e+ e- -> mu+ mu-.
    """
    def __init__(self, p1, p2, p3, p4):
        """
        Initialize the
        """
        from math import pi
        self.p1 = p1
        self.p2 = p2
        self.p3 = p3
        self.p4 = p4
        self.dmu = DiracMatrices()
        self.dml = ~self.dmu
        # Calculate the cross-section prefactor ((hbar c)/(8 pi))^2 in
        # units m^2 GeV^2.
        self.xspre = (1.97326979e-16/(8*pi))**2
        # Calculate the matrix-element prefactor (-4 pi alpha).
        self.mepre = -4*pi/137.0
    def me(self):
        """
        Return the matrix element given the state of the internally
        represented particles.
        """
        p0 = self.p1.p + self.p2.p
        return self.mepre/p0**2*sum([
            (self.p3.wbar()*self.dmu[mu]*self.p4.w())*
            (self.p2.wbar()*self.dml[mu]*self.p1.w())
            for mu in range(0, 4)])
    def xs(self, phi, theta):
        """
        Return the cross-section in m^2 for a given phi and theta.
        """
        from math import sqrt, cos, sin
        ct = cos(theta)
        st = sin(theta)
        q = sqrt(self.p1.p[0]**2 - self.p3.data.mass**2)
        p = sqrt(sum([self.p1.p[mu]**2 for mu in range(1, 4)]))
        self.p3.p[0] = self.p1.p[0]
        self.p3.p[1] = q*st*cos(phi)
        self.p3.p[2] = q*st*sin(phi)
        self.p3.p[3] = q*ct
        self.p4.p = ~self.p3.p
        me = self.me()
        try: me2 = me.real**2 + me.imag**2
        except: me2 = me**2
        return self.xspre*me2*q/p*st/(self.p1.p[0] + self.p2.p[0])**2

