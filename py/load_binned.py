import ttm.sm2015 as sm15

print("Loading binned data...")
smn = sm15.load('SM_Nose', bin=True)
smp = sm15.load('SM_wing_Port', bin=True)
sms = sm15.load('SM_wing_Star', bin=True)
print("Done.")
