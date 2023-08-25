import os
import numpy as np

ellList=np.arange(40,380,20)
ellList=np.append(np.array([50]),ellList, axis=0)
ellList = [100] #,380]
#ellList = [40, 60, 140, 160, 200, 220, 240, 260, 280, 300, 320, 340] #][20,30,50,80,120,150,180] #,180]
for ell in ellList:
    cmd = 'python pythonFiles/addRecDimAtEll_smth.py --ell={0:d}'.format(ell)
    print(cmd)
    os.system(cmd)

