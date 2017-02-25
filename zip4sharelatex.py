import zipfile
import glob

z = zipfile.ZipFile('Kilcher_etal_IMU-ADV.zip', 'w', compression=8)

figfiles = ['map04_annot',
            'TTM_Simple',
            'TTM_image01_annot',
            'StableMoor_Composite',
            'Torpedo_Image01_annot',
            'stationary_noise04',
            'SpecFig04_filtering2',
            'TimeFig02',
            'SpecFig02_TTM02B-top',
            'SpecFig02_SMnose',
            'SpecFig03_TTT',
            'StressSpec_TTM_04vp',
            'CoSpecND03_TTM-both',
            'TurbTime_TTM_01',
            'EpsVProd01',
            'EpsVU_03',
]

files = glob.glob('*.tex')
files += ['fig/' + fnm + '.pdf' for fnm in figfiles]
files += ['ametsoc.cls',
          'defs.sty',
          'thisbib.bib']

for fl in files:
    z.write(fl)

z.close()
