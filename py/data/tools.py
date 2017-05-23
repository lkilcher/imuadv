import shutil
import os.path as path
try:
    from urllib.request import urlopen
except:
    from urllib import urlopen

datdir = path.dirname(path.realpath(__file__))


def retrieve(url, fname):
    response = urlopen(url)
    with open(fname, 'wb') as f:
        shutil.copyfileobj(response, f)
