from base import tbx, datdir
import process_data as pd

mhkdr = 'http://mhkdr.openei.org/files/'


def pull():
    print("Retrieving SMB nose-mode data...")
    tbx.retrieve(
        mhkdr + '51/StableMoor_Nose_NREL02_May2015.VEC',
        datdir + 'SMN.VEC',
        hash='707ae7fa06c55d20'
    )
    print("Retrieving SMB Bottom-Track data...")
    tbx.retrieve(
        mhkdr + '223/BT_IMU_sync_data.mat',
        datdir + 'BT_IMU_sync_data.mat',
        hash='109bcde2ef771de9'
    )
    print("Retrieving Davit Torpedo data from May2015 (SMB deployment) ...")
    tbx.retrieve(
        mhkdr + '223/Turbulence_torpedo_data.VEC',
        datdir + 'TTT_Davit_B_May2015.VEC',
        hash='422831f6f5f1873f'
    )


def process_SMN():
    pd.pre_process_adv('SMN')
    pd.merge_adv_bt('SMN', filt_freqs={'5s': 0.2})


def process_torpedo():
    pd.pre_process_adv('TTT_Davit_B_May2015')
    pd.process_basic('TTT_Davit_B_May2015')
