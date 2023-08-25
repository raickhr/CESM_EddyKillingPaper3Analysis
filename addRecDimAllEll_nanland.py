import os
import numpy as np

ellList=np.arange(40,380,20)
ellList=np.append(np.array([50]),ellList, axis=0)
#ellList = [300,380]
ellList = [100] #,180]
for ell in ellList:
    cmd = 'python pythonFiles/addRecDimAtEll_nanland.py --ell={0:d}'.format(ell)
    print(cmd)
    os.system(cmd)

