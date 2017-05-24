"""
The following figures are used by the document:

01+ map04_annot.pdf
02- TTM_Simple
03- TTM_image01_annot
04- StableMoor_Composite
05- Torpedo_Image01_annot
06+ stationary_noise04
07+ SpecFig04_filtering2
08+ TimeFig02
09+ SpecFig02_TTM02B-top
10+ SpecFig02_SMnose
11+ SpecFig03_TTT
12+ StressSpec_TTM_04vp
13+ CoSpecND03_TTM-both
14+ TurbTime_TTM_01
15+ EpsVProd01
16+ EpsVU_03

Above:
  '+' indicates there is a script here that produces a file
  '-' indicates that there is not (these are images)

"""

import ptools as pt

# This is set to ion() in ptools, but when running this 'make all'
# script we don't want to display all of the figures. Just generate
# them.
pt.plt.ioff()


############
# # FIGURE 1
# # This is the 'base file' for the map (before annotations were added
# # in inkscape)

# execfile('make_map04.py')


############
# FIGURES 2-5 are images, with inkscape annotations.
###

############
# # FIGURE 6
# # This is the 'base file' for the map (before annotations were added
# # in inkscape)

execfile('imu_uncertainty.py')
