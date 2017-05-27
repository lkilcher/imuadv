"""
This script generates all of the figures utilized in the publication.

The following figures are used by the document:

+ 01 map04_annot.pdf
- 02 TTM_Simple
- 03 TTM_image01_annot
- 04 StableMoor_Composite
- 05 Torpedo_Image01_annot
+ 06 stationary_noise04
+ 07 SpecFig04_filtering2
+ 08 TimeFig02
+ 09 SpecFig02_TTM02B-top
+ 10 SpecFig02_SMnose
+ 11 SpecFig03_TTT
+ 12 StressSpec_TTM_04vp
+ 13 CoSpecND03_TTM-both
+ 14 TurbTime_TTM_01
+ 15 EpsVProd01
+ 16 EpsVU_03

Above:
  '+' indicates there is a script here that produces a file
  '-' indicates that there is not (these are images)

"""

import ptools as pt
from subprocess import call
import os
import sys

# This is set to ion() in ptools, but when running this 'make all'
# script we don't want to display all of the figures. Just generate
# them.
pt.plt.ioff()


def cp_env():
    e = os.environ.copy()
    path = ''
    for p in sys.path:
        if p == '':
            continue
        path += p + ':'
    path = path[:-1]  # Drop the last ":"
    e['PYTHONPATH'] = path
    return e


def run(script_name):
    call(['python', script_name],
         env=cp_env())

print("""
FIGURE 1: map04.pdf
This is the starting file for the map; annotations were added using
inkscape.
""")
run('make_map04.py')

###############################
# FIGURES 2-5 are images, with inkscape annotations.
###

print("""
FIGURE 6: stationary_noise04.pdf
The IMU noise-level figure.
""")
run('imu_uncertainty.py')


print("""
FIGURE 7: SpecFig04_filtering2.pdf
This figure shows how spectra change with Accel high-pass filter frequency
""")
run('spec04.py')


print("""
FIGURE 8: TimeFig02.pdf
The mean velocity figure (shows agreement between AWAC and ADV)
""")
run('TTM_TimeFigs.py')


print("""
FIGURE 9: SpecFig02_TTM02B-top.pdf
The TTM spectra figure.
""")
run('spec02.py')


print("""
FIGURE 10: SpecFig02_SMnose.pdf
The SMB spectra figure.
""")
run('spec_SM01.py')


print("""
FIGURE 11: SpecFig03_TTT.pdf
The Torpedo spectra figure.
""")
run('spec03.py')


print("""
FIGURE 12: StressSpec_TTM_04vp.pdf
FIGURE 13: CoSpecND03_TTM-both.pdf
The cross-spectra ('stress-spec') figures
""")
run('stresses.py')


print("""
FIGURE 14: TurbTime_TTM_01.pdf
The turbulence time-series figure.
FIGURE 15: EpsVProd01.pdf
The epsilon-production scatter plot
""")
run('TTM_TKEfigs01.py')


print("""
FIGURE 16: EpsVU_03.pdf
The epsilon-U^3 scatter plots
""")
run('epsVU.py')
