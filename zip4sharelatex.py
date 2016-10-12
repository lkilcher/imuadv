import zipfile
import glob

z = zipfile.ZipFile('Kilcher_etal_IMU-ADV.zip', 'w', compression=8)

files = glob.glob('*.tex')
files += glob.glob('fig/*')
files += ['ametsoc.cls',
          'defs.sty',
          'all.bib']

for fl in files:
    z.write(fl)

z.close()
