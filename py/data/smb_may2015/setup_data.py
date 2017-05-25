from base import tbx, datdir
import process_data as pd


def pull():
    print("Retrieving SMB nose-mode data...")
    tbx.retrieve(
        '<NEED LINK HERE!>... dropbox?',
        datdir + 'SMN.VEC',
        hash='707ae7fa06c55d20'
    )
    print("Retrieving SMB Bottom-Track data...")
    tbx.retrieve(
        '<NEED LINK HERE!>... dropbox?',
        datdir + 'BT_IMU_sync_data.mat',
        hash='109bcde2ef771de9'
    )


def process():
    pd.pre_process_adv()
    pd.merge_adv_bt({'5s': 0.2})
