from utils import Matrix, FourVector, Vector

from ParticleDataBase import ParticleData, ParticleDatabase
from DiracMatrices import DiracMatrices
from CrossSection import Annihilate
from Integrator import Integrator
from Particle import Particle



def circle(x, y):
    """
    Return 1 if 'x' and 'y' in a unit circle, 0 otherwise.
    """
    from math import sqrt
    f = sqrt(1 - x**2)
    return 0 if abs(y) > f else 1

def show(expr):
    """
    Print and evaluate an expression.
    """
    print("%s\n%10s:\n%s\n%s\n" % ("-"*20, expr, "-"*20, eval(expr)))


if __name__== "__main__":
    
        # ParticleDatabase class.
    pdb = ParticleDatabase()
    show("pdb['e+']")
    show("pdb['e-']")

    # Check the anti-particles.
    pds = sorted([pd for pd in pdb.items() if type(pd[0]) == int])
    for pid, pd in pds:
        if pid > 0 and pid < 30 and not pd.anti: show("pd.name")

    # DiracMatrices class.
    dm = DiracMatrices()
    show("dm[0]")
    show("dm[1]")
    show("dm[2]")
    show("dm[3]")

    # Particle class.
    pp = FourVector(-1, 22.0, 150.0, 400)
    ph = 1
    p = Particle(pdb["e-"], pp, ph)
    show("p.w()")
    show("p.wbar()")
    p.h = -1
    show("p.w()")
    show("p.wbar()")

    # Diract momentum-space equation.
    m = Matrix([1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0],
                         [0, 0, 0, 1])*abs(p.p)
    w = p.wbar()
    t = None
    for i in range(0, 4):
        if not t: t = (dm[i]*(~p.p)[i] - m)*w
        else: t += (dm[i]*(~p.p)[i] - m)*w
    show("sum([ti for ti in t])")

    # Circle integration.
    i = Integrator(circle, -1, 1, -1, 1)
    show("i.mc()")

    # Helicity cross-sections.
    from math import pi
    p1 = FourVector(-1.0, 0.0, 0.0, 100)
    p2 = FourVector(-1.0, 0.0, 0.0, -100)
    p3 = FourVector(0.0, 0.0, 0.0, 0.0)
    p4 = FourVector(0.0, 0.0, 0.0, 0.0)
    pp1 = Particle(pdb["e-"], p1, 1)
    pp2 = Particle(pdb["e+"], p2, 1)
    pp3 = Particle(pdb["mu-"], p3, 1)
    pp4 = Particle(pdb["mu+"], p4, 1)
    a = Annihilate(pp1, pp2, pp3, pp4)
    i = Integrator(a.xs, 0.0, 2*pi, 0, pi)
    for h1 in [-1, 1]:
        pp1.h = h1
        for h2 in [-1, 1]:
            pp2.h = h2
            for h3 in [-1, 1]:
                pp3.h = h3
                for h4 in [-1, 1]:
                    pp4.h = h4
                    print("%2i %2i %2i %2i %8.1e" % (
                        h1, h2, h3, h4, i.mc(1000)/1e-31))

    # Cross-section integration.
    try: import pythia8; py = pythia8.Pythia("", False)
    except: py = None

    for p in [5.0, 10.0, 50.0, 100.0, 1000.0]:

        # Problem set cross-section.
        from math import pi
        p1 = FourVector(-1.0, 0.0, 0.0, p)
        p2 = FourVector(-1.0, 0.0, 0.0, -p)
        p3 = FourVector(0.0, 0.0, 0.0, 0.0)
        p4 = FourVector(0.0, 0.0, 0.0, 0.0)
        pp1 = Particle(pdb["e-"], p1, 1)
        pp2 = Particle(pdb["e+"], p2, 1)
        pp3 = Particle(pdb["mu-"], p3, 1)
        pp4 = Particle(pdb["mu+"], p4, 1)
        a = Annihilate(pp1, pp2, pp3, pp4)
        i = Integrator(a.xs, 0.0, 2*pi, 0, pi)
        xs0 = 0
        for h3 in [-1, 1]:
            pp3.h = h3
            for h4 in [-1, 1]:
                pp4.h = h4
                xs = []
                for h1 in [-1, 1]:
                    pp1.h = h1
                    for h2 in [-1, 1]:
                        pp2.h = h2
                        xs += [i.mc(1000)/1e-31]
                xs0 += sum(xs)/len(xs)
        show("xs0")

        # Pythia 8 cross-section.
        if py:
            py.readString("Print:quiet = on")
            py.readString("Beams:idA = 11")
            py.readString("Beams:idB = -11")
            py.readString("Beams:frameType = 3")
            py.readString("Beams:pzA = %r" % p)
            py.readString("Beams:pzB = -%r" % p)
            py.readString("PDF:lepton = off")
            py.readString("PartonLevel:all = off")
            py.readString("SigmaProcess:alphaEMorder = 0")
            py.readString("StandardModel:alphaEM0 = %r" % (1.0/137.0))
            py.readString("WeakSingleBoson:ffbar2ffbar(s:gm) = on")
            py.init()
            acc = 0
            for i in range(0, 10000):
                py.next()
                if py.process[5].idAbs() == 13: acc += 1
            xs1 = float(acc)/py.info.nAccepted()*py.info.sigmaGen()
        else: xs1 = 1.0
        show("p, xs1/xs0")





