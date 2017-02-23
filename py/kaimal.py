import numpy as np


class Kaimal(object):

    def __init__(self, f=np.arange(1e-2, 10, 1e-2)):
        self.freq = f

    def Suu(self, U=1, z=1, ustar2=1):
        f0 = U / z
        fstar = self.freq / f0
        return (105 / (1 + 33 * fstar) ** (5. / 3.)) * (ustar2 / f0)

    def Sww(self, U=1, z=1, ustar2=1):
        f0 = U / z
        fstar = self.freq / f0
        return (2 / (1 + 5.3 * (fstar) ** (5. / 3.))) * (ustar2 / f0)

    def Suw(self, U=1, z=1, ustar2=1):
        f0 = U / z
        fstar = self.freq / f0
        return (14.0 / (1 + 9.6 * fstar) ** (2.4)) * (ustar2 / f0)

    def __getitem__(self, idx):
        if idx == 0:
            return self.Suu()
        elif idx == 2:
            return self.Sww()
        elif idx == 1:
            return None
        raise IndexError("Invalid index.")
