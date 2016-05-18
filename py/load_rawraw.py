import ttm.sm2015 as sm15

print("Loading raw RAW data...")
smn = sm15.load('SM_Nose', coordsys='raw', bin=False)
smp = sm15.load('SM_wing_Port', coordsys='raw', bin=False)
sms = sm15.load('SM_wing_Star', coordsys='raw', bin=False)
print("Done.")
