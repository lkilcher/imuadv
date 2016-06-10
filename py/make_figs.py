
def execfile(fname, global_vars=None, local_vars=None):
    with open(file) as f:
        code = compile(f.read(), fname, 'exec')
        exec(code, global_vars, local_vars)

# execfile('imu_uncertainty.py')
