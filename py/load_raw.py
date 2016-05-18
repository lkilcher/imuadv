import ttm.sm2015 as sm15

print("Loading raw data...")
smn = sm15.load('SM_Nose', bin=False)
smp = sm15.load('SM_wing_Port', bin=False)
sms = sm15.load('SM_wing_Star', bin=False)
print("Done.")
