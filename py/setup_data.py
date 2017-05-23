"""
The following figures are used by the document:

X map04_annot.pdf
X TTM_Simple
X TTM_image01_annot
X StableMoor_Composite
X Torpedo_Image01_annot
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


The following data files are needed:

# bench test
C

# June2012
AWAC
TTM / NREL:raw,bin both pax

# June2014
ttm02b top + bot
ttm01b top + bot

# SM2015
Nose mode files. Others?

"""

import data.bathy
data.bathy.pull_data()
data.bathy.process_data()

