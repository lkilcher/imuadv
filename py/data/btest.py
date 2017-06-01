import dolfyn.adv.api as avm
try:
    from . import tools as tbx
except ValueError:
    import tools as tbx

fnm = 'btest-C'


def pull():
    print("Retrieving btest data...")
    tbx.retrieve(
        'http://mhkdr.openei.org/files/223/Bench_test_data.VEC',
        tbx.datdir + fnm + '.VEC',
        hash='d4a1e3555a40e13b'
    )


def process():
    print("Processing btest data...")
    outfile = tbx.datdir + fnm + '.h5'
    if tbx.checkhash(outfile, 'f42d5ca36cc324b5'):
        print("   hash check passed; skipping processing.")
        return
    # Skip the last 600 points or so because there is a NaN in the IMU
    # data there for some reason.
    dat = avm.read_nortek(tbx.datdir + fnm + '.VEC', npings=3348490)
    dat.save(outfile)
