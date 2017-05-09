import shutil
try:
    from urllib.request import urlopen
except:
    from urllib import urlopen


def retrieve(url, fname):
    response = urlopen(url)
    with open(fname, 'wb') as f:
        shutil.copyfileobj(response, f)
