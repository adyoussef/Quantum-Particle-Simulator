

class Integrator:
    """
    This class integrates a two variable function.
    """
    def __init__(self, f, xmin, xmax, ymin, ymax):
        """
        Initialize the integrator, given a function 'f', a minimum x
        'xmin', a maximum x 'xmax', a minumum y 'ymin', and a
        maxumimum y 'ymax.
        """
        self.f = f
        self.xmin = xmin
        self.xmax = xmax
        self.ymin = ymin
        self.ymax = ymax
        self.xdif = xmax - xmin
        self.ydif = ymax - ymin
    def mc(self, n = 1000):
        """
        Perform MC integration for given number of sampling points 'n'.
        """
        import random
        t = 0
        for i in range(0, n):
            x = self.xmin + random.random()*self.xdif
            y = self.ymin + random.random()*self.ydif
            t += self.f(x, y)
        return t/float(n)*self.xdif*self.ydif







