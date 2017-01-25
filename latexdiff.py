#!/usr/bin/python
import glob
import sys
from subprocess import call
import os
from shutil import copyfile

mainfile = 'Kilcher_etal_IMU-ADV'

if sys.argv < 2:
    raise Exception("You must specify the root we're comparing against.")

indir = './diff/'

nm1 = sys.argv[1]
dir1 = indir + nm1 + '/'

if len(sys.argv) > 2:
    nm2 = sys.argv[2]
    dir2 = indir + nm2 + '/'
else:
    nm2 = 'now'
    now_dir = './'

if nm1.startswith('d_') or nm2.startswith('d_'):
    raise Exception("It looks like you're trying to latexdiff a latexdiff!?!!")

outdir = indir + '/diff_' + nm1 + '-' + nm2 + '/'

try:
    os.stat(outdir)
except:
    os.mkdir(outdir)
print outdir

for fl in glob.glob(dir1 + '*.tex'):
    fnm = fl.rsplit('/', 1)[-1]
    outfile = open(outdir + fnm, "w")
    tmp = call(['latexdiff', fl, now_dir + fnm], stdout=outfile)
    outfile.close()
    print fl, fnm

for fl in ['fig', 'ametsoc.cls', 'ametsoc2014.bst', 'thisbib.bib', 'defs.sty']:
    try:
        os.remove(outdir + fl)
    except:
        pass
    os.symlink(os.path.abspath('./' + fl), outdir + fl)

os.chdir(outdir)

call(['pdflatex', '--interaction', 'nonstopmode', mainfile + '.tex'])
call(['bibtex', mainfile + '.aux'])
call(['pdflatex', '--interaction', 'nonstopmode', mainfile + '.tex'])
call(['pdflatex', '--interaction', 'nonstopmode', mainfile + '.tex'])
copyfile(mainfile + '.pdf', '../diff_' + nm1 + '-' + nm2 + '.pdf')
os.chdir('../../')
