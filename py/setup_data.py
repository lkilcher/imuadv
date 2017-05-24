"""
The following figures are used by the document:

stationary_noise04
SpecFig04_filtering2
TimeFig02
SpecFig02_TTM02B-top
SpecFig02_SMnose
SpecFig03_TTT
StressSpec_TTM_04vp
CoSpecND03_TTM-both
TurbTime_TTM_01
EpsVProd01
EpsVU_03


The following data files are still needed:

# June2012
AWAC
TTM / NREL:raw,bin both pax

# June2014
ttm02b top + bot
ttm01b top + bot

# SM2015
Nose mode files. Others?
!!! Currently I'm using the SM2015 TTM for imu_uncertainty fig. Need to fix this.

"""

import data.bathy
import data.btest

print("""
This script downloads and pre-processes data files for making the
figures. This can take some time (several hours, depending on your
internet connection), so be ready to be patient.
""")


###
# This is the bathymetry data for Figure 1
data.bathy.pull()
data.bathy.process()

###
# This is the bench-test data, for Figure 6
data.btest.pull()
data.btest.process()

