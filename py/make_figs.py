import ptools as pt

# This is set to ion() in ptools, but when running this 'make all'
# script we don't want to display all of the figures. Just generate
# them.
pt.plt.ioff()


def execfile(fname, global_vars=None, local_vars=None):
    with open(fname) as f:
        code = compile(f.read(), fname, 'exec')
        exec(code, global_vars, local_vars)

###
# This is the 'base file' for the map (before annotations were added
# in inkscape)

# execfile('make_map04.py')
