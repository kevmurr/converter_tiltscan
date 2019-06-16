#set these values first!
detector_dist=0.3
scanmotor="yML2"
coord_of_db=(480,680) #(478,800) (482,681)
target="Mo"
appendix=None #Put None if none needed
###################################################
cut_directbeam=True
db_window=10

import numpy as np
import sys
import os
import h5py as h5


motor_list=["yML1","pML1","yML2","pML2"]
target_list={
    "Mo": 17.48,
    "Cu": 8,
    "Ru": 20.2,
}

flip_dict={
    "yML1":False,
    "pML1":False,
    "yML2":True, #tested
    "pML2":True, #tested
}
energy=target_list[target]
if scanmotor not in motor_list:
    IOError("Motorname not recognized")
flip_bool=flip_dict[scanmotor]
if scanmotor.startswith("y"):
    orientation="h"
    roi = (coord_of_db[0] - 7, coord_of_db[0] + 7, 0, 1555)
else:
    orientation="v"
    roi = (0, 515, coord_of_db[1] - 7, coord_of_db[1] + 7)

example_fullname=str(sys.argv[1])
save_dir=str(sys.argv[2])+"/"
datfile=str(sys.argv[3])
example_directory=os.path.split(example_fullname)[0]+"/"
example_filename=os.path.split(example_fullname)[1]
scanname=example_filename[:7]
print("Scanname: {}".format(scanname))
allfiles=os.listdir(example_directory)

scanfiles=[]
for x in allfiles:
    y=os.path.split(x)[1]
    if y.startswith(scanname):
        scanfiles.append(x)
scanfiles.sort()
if orientation=="h":
    for i in range(0,len(scanfiles),1):
        filename_now=scanfiles[i]
        print("Loading "+example_directory+filename_now)
        f_now=h5.File(example_directory+filename_now,"r")
        data_now=f_now["/entry/instrument/detector/data"][0,:,:]
        data_roi=data_now[roi[0]:roi[1],roi[2]:roi[3]]
        if i==0:
            tilt_array=np.zeros((len(scanfiles),data_roi.shape[1]))
        sumarray=np.sum(data_roi,axis=0)
        if cut_directbeam==True:
            sumarray[coord_of_db[1]-db_window:coord_of_db[1]+db_window]=np.zeros((2*db_window))
        tilt_array[i,:]=sumarray
        f_now.close()
    thetafac=360*55*10**-6/(detector_dist*2*np.pi)
    thetamax=thetafac*(coord_of_db[1]-roi[2])
    thetamin=thetafac*(coord_of_db[1]-roi[3])
    thetaarr=np.linspace(thetamax,thetamin,tilt_array.shape[1])

if orientation=="v":
    for i in range(0,len(scanfiles),1):
        filename_now=scanfiles[i]
        print("Loading "+example_directory+filename_now)
        f_now=h5.File(example_directory+filename_now,"r")
        data_now=f_now["/entry/instrument/detector/data"][0,:,:]
        data_roi=data_now[roi[0]:roi[1],roi[2]:roi[3]]
        if i==0:
            tilt_array=np.zeros((len(scanfiles),data_roi.shape[0]))
        sumarray=np.sum(data_roi,axis=1)
        if cut_directbeam==True:
            sumarray[coord_of_db[0]-db_window:coord_of_db[0]+db_window]=np.zeros((2*db_window))
        tilt_array[i,:]=sumarray
    thetafac=360*55*10**-6/(detector_dist*2*np.pi)
    thetamax=-thetafac*(coord_of_db[0]-roi[0])
    thetamin=-thetafac*(coord_of_db[0]-roi[1])
    thetaarr=np.linspace(thetamax,thetamin,tilt_array.shape[1])

translation_file=np.loadtxt(datfile,skiprows=24)
translation_vec=np.zeros((translation_file.shape[0],1))
translation_vec=translation_file[:,1]
#if flip==True:
#    translation_vec=-translation_vec
#    print("Flipping")


if appendix!=None:
    f=h5.File(save_dir+scanname+"_energy_scan_"+str(appendix)+".h5","w")
else:
    f=h5.File(save_dir+scanname+"_energy_scan.h5","w")

print("Saving as %s"%(save_dir+scanname+"_tilt_scan.h5"))
if flip_bool==True:
    tilt_array=np.flip(tilt_array,axis=0)
data=f.create_dataset("Data",data=tilt_array)
positions=f.create_dataset("Omega",data=translation_vec)
theta=f.create_dataset("2Theta",data=thetaarr)
energy=f.create_dataset("Energy",data=energy)
orientation_file=f.create_dataset("Orientation",data=str(orientation))

f.close()        
print("Done")