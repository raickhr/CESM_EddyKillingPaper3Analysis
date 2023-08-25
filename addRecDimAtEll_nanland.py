import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument(
    '--ell', help="FilterLength", type=int, default=40)

args = parser.parse_args()
ell = args.ell
#.nc_100._Filtered.nc
#tavg_365_371.nc_0040_Filtered.nc
for i in range(10,12):#366):
    timeval = 49*365 + 7*i + 3
    dimInfo = 'defdim("time",1);time[time]={0:f};time@long_name="Time";time@units="days since 0000-01-01 00:00:00"'.format(timeval)
    fileName = '{1:d}km_nanland/{0:03d}.nc_{1:04d}_Filtered.nc'.format(i,ell)
    cmd = "ncap2 -s '{0:s}' -O {1:s} {1:s}".format(dimInfo, fileName)
    cmd2 = 'ncks -O --mk_rec_dmn time {0:s} {0:s}'.format(fileName)
    cmd4 = "ncecat -O -u time {0:s} {0:s}".format(fileName)
    cmd3 = "ncwa -O -a time {0:s} {0:s}".format(fileName)
    
    print(cmd)
    os.system(cmd)

    print(cmd2)
    os.system(cmd2)
    
    print(cmd3)
    os.system(cmd3)

    print(cmd4)
    os.system(cmd4)




