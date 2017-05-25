from base import tbx, datdir


def pull():
    print("Retrieving SMB nose-mode data...")
    tbx.retrieve(
        '<NEED LINK HERE!>',
        datdir + '<NEED FNAME>.VEC',
        hash='<NEED HASH HERE>'
    )
    print("Retrieving SMB Bottom-Track data...")
    tbx.retrieve(
        '<NEED LINK HERE!>',
        datdir + '<NEED FNAME>.mat',
        hash='<NEED HASH HERE>'
    )
