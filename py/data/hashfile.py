import tools as tbx
import sys

sys.argv.pop(0)

for fnm in sys.argv:

    h = tbx.sha(fnm)

    print("Hash of file: '{}'".format(fnm))
    print(h)
    print("")
