"""
The following figures are used by the document:

01+ map04_annot.pdf
  :map04.py
  . bathy
02- TTM_Simple
03- TTM_image01_annot
04- StableMoor_Composite
05- Torpedo_Image01_annot
06+ stationary_noise04
  :imu_uncertainty.py
  . btest-C
  . !!! sm15-TTM-VELMOOR (switch to J14?)
07+ SpecFig04_filtering2
  :spec04.py
  . j14.ttm02b-top.raw
08+ TimeFig02
  :TTM_TimeFigs.py
  . j12.awac.bnd
  ? j12.awac.raw
  . j12.nrel-pax.bnd
  ? j12.nrel-pax.raw
09+ SpecFig02_TTM02B-top
  :spec02.py
  . dat = j14.load('ttm02b-top', 'pax',
                   bin=True)
10+ SpecFig02_SMnose
  :spec_SM01.py
  . SMN data (reorganize)
11+ SpecFig03_TTT
  :spec03.py
  . sm15.ttt_davit_b-pax.bnd
12+ StressSpec_TTM_04vp
  stresses.py
  . j14.load('ttm02b-top', 'pax', bin=True)
  ? j14.load('ttm02b-bottom', 'pax', bin=True)
  ? smdat.load('SMN-5s', bindat=True)
  . TTT avm.load('/Users/lkilcher/data/pnnl/TTT_Vector_Feb2011_pax_b5m.h5')
13+ CoSpecND03_TTM-both
  :stresses.py
  ... Same as 12
14+ TurbTime_TTM_01
  :TTM_TKEfigs01.py
  . dat = j14.load('ttm02b-top', 'pax',
                   bin=True)
  ? dat2 = j14.load('ttm02b-bot', 'pax',
                    bin=True)
  ? dat3 = j14.load('ttm01b-top', 'pax',
                    bin=True).subset(slice(3, -2))
  ? dat4 = j14.load('ttm01b-bot', 'pax',
                    bin=True).subset(slice(3, -2))
15+ EpsVProd01
  :TTM_TKEfigs01.py
  . dat = j14.load('ttm02b-top', 'pax',
                   bin=True)
  . dat2 = j14.load('ttm02b-bot', 'pax',
                    bin=True)
  ? dat3 = j14.load('ttm01b-top', 'pax',
                    bin=True).subset(slice(3, -2))
  ? dat4 = j14.load('ttm01b-bot', 'pax',
                    bin=True).subset(slice(3, -2))
16+ EpsVU_03
  :epsVU.py
  . dat['ttm'] = j14.load('ttm02b-top', 'pax',
                          bin=True)
  . dat['sm'] = sm.load('SMN-5s', bindat=True)

The following data files are still needed:

# June2012
AWAC
TTM / NREL:raw,bnd both pax

# June2014
ttm02b top + bot
ttm01b top + bot

# SM2015
Nose mode files. Others?
!!! Currently I'm using the SM2015 TTM for imu_uncertainty fig. Need to fix this.

"""

import data.bathy
import data.btest
import data.ttm_june2012.setup_data as j12setup
import data.ttm_june2014.setup_data as j14setup
import data.smb_may2015.setup_data as m15setup

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

###
# This is the June 2012 TTM dataset
j12setup.pull()
j12setup.process()

###
# This is the June 2014 TTM dataset
j14setup.pull()
j14setup.process()
j14setup.ttmlean.process()

###
# This is the May 2015 SMB dataset
m15setup.pull()
m15setup.process_SMN()
m15setup.process_torpedo()
