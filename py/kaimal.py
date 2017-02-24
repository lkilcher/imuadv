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


def nd_cospec(dat, inds):
    z = dat.z
    Uhor = np.abs(dat.U[inds])[:, None]
    f0 = Uhor / z
    ustar2 = np.sqrt(dat.upwp_[inds] ** 2 +
                     dat.vpwp_[inds] ** 2)[:, None]
    # ustar2 = np.abs(dat.upwp_[inds])[:, None]
    fstar = dat.freq / f0
    Cdata = (-np.sign(dat.u[inds][:, None]) *
             dat['Cspec_vel'][1][inds].real *
             f0 / ustar2)
    return Cdata, fstar

# def nd_cospec(dat, inds):
#     z = dat.z
#     Uhor = np.abs(dat.U[inds])[:, None]
#     fstar = dat.freq * z / Uhor
#     #ustar2 = (kappa * z) ** 2 * dat.S2[inds, None]
#     ustar2 = (kappa * z * dat.dudz[inds, None]) ** 2
#     Cdata = (-np.sign(dat.u[inds][:, None]) *
#              dat['Cspec_vel'][1][inds].real *
#              Uhor / (ustar2 * z))
#     return Cdata, fstar

# def nd_cospec(dat, inds):
#     z = dat.z
#     Uhor = np.abs(dat.U[inds])[:, None]
#     fstar = dat.freq * z / Uhor
#     # The 0.0659 is from a fit of dudz to u.
#     ustar2 = dat.u[inds, None] ** 2 / 77400.
#     # ustar2 = (0.005 * dat.u[inds, None]) ** 2
#     ustar2 *= (kappa * z) ** 2
#     Cdata = (-np.sign(dat.u[inds][:, None]) *
#              dat['Cspec_vel'][1][inds].real *
#              Uhor / (ustar2 * z))
#     return Cdata, fstar


def bin_cospec(Cdin, fsin, fbins):
    nb = len(fbins)
    df = np.diff(fbins)
    fbe = np.empty(nb + 1)
    fbe[1:-1] = fbins[:-1] + df / 2
    fbe[0] = fbins[0] - df[0] / 2
    fbe[-1] = fbins[-1] + df[-1] / 2
    Cs = np.empty_like(fbins)
    for idx in range(nb):
        itmp = (fbe[idx] < fsin) & (fsin < fbe[idx + 1])
        Cs[idx] = Cdin[itmp].mean()
    return Cs
