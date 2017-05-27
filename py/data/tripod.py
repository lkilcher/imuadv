import tools as tbx
import dolfyn.adv.api as avm

fname = 'tripod_data'


def pull():
    print("Retrieving btest data...")
    tbx.retrieve(
        'https://www.dropbox.com/s/qppgv82xwxrcs0u/TTT_Vector_Feb2011.vec?dl=1',
        tbx.datdir + fname + '.VEC',
        hash='d95c230aa67f53ae'
    )


def process():

    # Something strange happens with the measurements after 11800000
    dat = avm.read_nortek(tbx.datdir + fname + '.VEC', npings=11800000)
    dat.props['declination'] = 16.95  # Degrees East.
    dat.props['heading_offset'] = 9  # Degrees CCW
    # Yes, this is correct, when the head points up, the ADV is 'down'!?
    dat.config.orientation = 'down'

    avm.rotate._inst2earth(dat, use_mean_rotation=True)

    avm.clean.GN2002(dat)

    avm.rotate.earth2principal(dat)

    binner = avm.TurbBinner(n_bin=5 * 60 * dat.fs, fs=dat.fs)

    bd = binner(dat)

    bd.add_data('Cspec_vel', binner.calc_vel_cpsd(dat.vel), 'spec')

    bd.save(tbx.datdir + fname + '_pax_b5m.h5')


def load():
    return avm.load(tbx.datdir + fname + '_pax_b5m.h5')
