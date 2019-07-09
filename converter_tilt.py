


import numpy as np
import sys
import os
import h5py as h5
##################################################
#Parameters
dbeam_cut_range=20
roi_size=30
year="2019"
##################################################

logdir = "/gpfs/cfel/cxi/labs/MLL-Sigray/scan-logs/"
scanssdir = "/gpfs/cfel/cxi/labs/MLL-Sigray/scan-frames/"
save_dir= "/gpfs/cfel/cxi/labs/MLL-Sigray/Processed/{0}/".format(year)
mask="/gpfs/cfel/cxi/labs/MLL-Sigray/mask/lambda_mask1.h5"

scanno=str(input("What is the scan number? (e.g. 10)"))
scan_name="Scan_{0}".format(scanno)
ldir=scanssdir+"{0}/".format(scan_name)

datfile=open("{0}{1}.log".format(logdir,scan_name))
i=0

for line in datfile:
    if line.startswith("# Points count"):
        N_points=int(line.split(":")[1])
        positions = np.zeros((N_points)).astype("float")
        useframes = np.ones((N_points)).astype(np.bool)
    if line.startswith("# Device:") and line.endswith("Scanner")==False and line.endswith("Lambda\n")==False:
        scanmotor=str(line.split(":")[1][1:])
    if line.startswith("#")==False:
        current_entry=line.split(";")[2]
        if i==0:
            unit=current_entry[-5:]
            print(unit)
        current_pos=float(current_entry[:-5])
        positions[i]=current_pos
        i+=1
useframes[i:]=np.zeros_like(useframes[i:]).astype(np.bool)

target=str(input("What target was used? (Mo,Cu,Rh)"))
targets=["Mo","Cu","Rh"]
if target not in targets:
    print("Error! Invalid target selected! Has to be Mo, Cu or Rh.")
    quit()


detector_dist=float(input("What was the detector distance (meters)?"))

target_list={
    "Mo": 17.48,
    "Cu": 8,
    "Rh": 20.2,
}

flip_dict={
    "Yaw-LENSE-UP\n":False,
    "Pitch-LENSE-UP\n":False,
    "Yaw-LENSE-DOWN\n":True, #tested
    "Pitch-LENSE-DOWN\n":True, #tested
}
energy=target_list[target]

flip_bool=flip_dict[scanmotor]

if scanmotor.startswith("Yaw"):
    orientation="h"
else:
    orientation="v"

maskpath="/data"
maskfile=h5.File(mask,"r")
mask=maskfile[maskpath][()].astype(np.int)


startframe=0
for i in range(0,N_points):
    try:
        if i==startframe:
            file_now="{0}_".format(scan_name)+"{:05.0f}_Lambda.nxs".format(i)
            example_f=h5.File(ldir+"/"+file_now,"r")
            example_data=example_f["/entry/instrument/detector/data"][0,:,:]
            print(np.sum(example_data))
            if np.sum(example_data)<10:
                startframe+=1
                print("First frame was only zeros")
                useframes[i]=False
            else:
                db_coord=np.unravel_index(np.argmax(np.multiply(mask,example_data)),example_data.shape)
                print("Direct beam pixel coordinate is {0}.".format(db_coord))
                roi=(db_coord[0]-int(roi_size/2),db_coord[0]+int(roi_size/2),db_coord[1]-int(roi_size/2),db_coord[1]+int(roi_size/2))
                print("Startframe is {0}. Creating full data array.".format(i))
                if orientation=="h":
                    data_full=np.zeros((N_points,example_data.shape[1]))
                else:
                    data_full = np.zeros((N_points, example_data.shape[0]))
        else:
            file_now = "/{0}_".format(scan_name) + "{:05.0f}_Lambda.nxs".format(i)
            f_now=h5.File(ldir+"/"+file_now,"r")
            if orientation=="h":
                data_now=np.multiply(f_now["/entry/instrument/detector/data"][0,:,:],mask)
                line_now=np.sum(data_now[roi[0]:roi[1],:],axis=0)
                line_now[(db_coord[1]-int(dbeam_cut_range/2)):(db_coord[1]+int(dbeam_cut_range/2))]=np.zeros_like(line_now[(db_coord[1]-int(dbeam_cut_range/2)):(db_coord[1]+int(dbeam_cut_range/2))])
            else:
                data_now = np.muliply(f_now["/entry/instrument/detector/data"][0, :, :], mask)
                line_now = np.sum(line_now[:, roi[2]:roi[3]],axis=1)
                line_now[(db_coord[0] - int(dbeam_cut_range / 2)):(db_coord[0] + int(dbeam_cut_range / 2))] = np.zeros_like(line_now[(db_coord[0] - int(dbeam_cut_range / 2)):(db_coord[0] + int(dbeam_cut_range / 2))])
            data_full[i,:]=line_now
    except (KeyError,OSError):
        if i==startframe:
            startframe+=1
        print("Didnt find file {} or proper data in file.".format(file_now))
        useframes[i]=False


print(db_coord)

if orientation=="h":
    thetafac=360*55*10**-6/(detector_dist*2*np.pi)
    thetamax=thetafac*(db_coord[1])
    thetamin=thetafac*(db_coord[1]-data_full.shape[1])
    thetaarr=np.linspace(thetamax,thetamin,data_full.shape[1])

if orientation=="v":
    thetafac=360*55*10**-6/(detector_dist*2*np.pi)
    thetamax=-thetafac*(db_coord[0])
    thetamin=-thetafac*(db_coord[0]-data_full.shape[1])
    thetaarr=np.linspace(thetamax,thetamin,data_full.shape[1])

translation_vec=positions
#if flip==True:
#    translation_vec=-translation_vec
#    print("Flipping")

if os.path.isdir(save_dir+scan_name)==False:
    os.mkdir(save_dir+scan_name)
f=h5.File("{0}{1}/{1}_tiltscan.h5".format(save_dir,scan_name),"w")

print("Saving as {0}{1}/{1}_tiltscan.h5".format(save_dir,scan_name))
if flip_bool==True:
    data_full=np.flip(data_full,axis=0)
f.create_dataset("Data", data=data_full)
f.create_dataset("Omega", data=translation_vec)
f.create_dataset("2Theta", data=thetaarr)
f.create_dataset("Energy", data=energy)
f.create_dataset("Orientation", data=str(orientation))

f.close()        
print("Done")