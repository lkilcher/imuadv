import shutil
import os.path as path
import hashlib
try:
    from urllib.request import urlopen
except:
    from urllib import urlopen

datdir = path.dirname(path.realpath(__file__)) + '/'


def sha(fname):
    h = hashlib.sha256()
    with open(fname, "rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()[:16]


def checkhash(fname, hash=None):
    if not path.isfile(fname) or hash is None:
        return None
    return sha(fname) == hash


def retrieve(url, fname, hash=None):
    val = checkhash(fname, hash)
    if val is None:
        pass
    elif checkhash(fname, hash):
        print("   Hash test passed; skipping download.")
        return
    else:
        print("   Hash test failed; "
              "overwriting existing file...")
    response = urlopen(url)
    with open(fname, 'wb') as f:
        shutil.copyfileobj(response, f)
