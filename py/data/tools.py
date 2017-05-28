import shutil
import os.path as path
import hashlib
import ssl
try:
    from urllib.request import urlopen
except:
    from urllib import urlopen

datdir = path.dirname(path.realpath(__file__)) + '/'

insecure_ctx = ssl.create_default_context()
insecure_ctx.check_hostname = False
insecure_ctx.verify_mode = ssl.CERT_NONE


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


def retrieve(url, fname, hash=None, verify=True):
    url_kwargs = {}
    if verify:
        url_kwargs['context'] = insecure_ctx
    val = checkhash(fname, hash)
    if val is None:
        pass
    elif checkhash(fname, hash):
        print("   Hash test passed; skipping download.")
        return
    else:
        print("   Hash test failed; "
              "overwriting existing file...")
    response = urlopen(url, **url_kwargs)
    with open(fname, 'wb') as f:
        shutil.copyfileobj(response, f)
