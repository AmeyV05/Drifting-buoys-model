# processing the vel and pos data from simulations. 
# This is a standalone script. 
import numpy as np
import os
import shutil
import pandas as pd
import matplotlib.pyplot as plt
np.set_printoptions(threshold=np.inf)
params = {'legend.fontsize': '11',
         'axes.labelsize': 'medium',
         'axes.titlesize':'x-large',}
plt.rcParams.update(params)


# reading the sim_data file from the folder

# buoy number input.

Bnum=input("Please input the Buoy number from the list: [02, 03, 07, 09, 12, 13, 14, 16]. And press enter: ")
print("Running for: Buoy_"+Bnum)

bname='BUOY_'+Bnum
loc='../../generated_data/'+bname+'/h1f1A1T1O1Po1Pt1/'

#savefile location.
test_loc=loc+'experiments/tidal+PG+cor/GTSM/'

#processing for filtered velocity data to create quiver plots.
file='simdata_'+bname+'.xlsx'
df=pd.read_excel(loc+file)
Xsfil=np.array(df['Unnamed: 2'][2:], dtype=float)
nnan=~np.isnan(Xsfil.astype(float))
Xsfil=Xsfil[nnan]
Ysfil=np.array(df['Unnamed: 3'][2:], dtype=float)
Ysfil=Ysfil[nnan]
Xbfil=np.array(df['Unnamed: 4'][2:], dtype=float)
Xbfil=Xbfil[nnan]
Ybfil=np.array(df['Unnamed: 5'][2:], dtype=float)
Ybfil=Ybfil[nnan]
Ubfil=np.array(df['Unnamed: 6'][2:], dtype=float)
Ubfil=Ubfil[nnan]
Vbfil=np.array(df['Unnamed: 7'][2:], dtype=float)
Vbfil=Vbfil[nnan]
Usfil=np.array(df['Unnamed: 8'][2:], dtype=float)
Usfil=Usfil[nnan];Usfil=Usfil[1:]
Vsfil=np.array(df['Unnamed: 9'][2:], dtype=float)
Vsfil=Vsfil[nnan];Vsfil=Vsfil[1:]
#normalisation of velocity fields. 
Ns=np.sqrt(Usfil**2+Vsfil**2)
UsfilN=Usfil/Ns;VsfilN=Vsfil/Ns
Nb=np.sqrt(Ubfil**2+Vbfil**2)
UbfilN=Ubfil/Nb;VbfilN=Vbfil/Nb
# Checking for dimension mismatch.
print("Length of the filtered array datsets [sim,obs]: ["+str(len(Usfil))+","+str(len(Ubfil))+"]")

print("Reading tidal dataset from main excel file.")
# Reading data from main pos_vel data. 
tloc='../../generated_data/'+bname+'/Pos_Vel_data.xlsx'
dft=pd.read_excel(tloc)

# numtaps to index dataset.
numtaps=24*2*2+1
edge=int(numtaps/2)
#indexing the data as the fitering process removes the start and end of the datasets. 
Tib=np.array(dft['Date(GMT)'][edge:-1*edge-1])  #-1 to remove one last extra element
Ut=np.array(dft['Ut'][edge:-1*edge-1])
Vt=np.array(dft['Vt'][edge:-1*edge-1])
#GTSM
Pgxt=np.array(dft['Pgxt'])[edge:-1*edge-1]
Pgyt=np.array(dft['Pgyt'])[edge:-1*edge-1]
#FES2014
Fpgx=np.array(dft['Fpgx'])[edge:-1*edge-1]
Fpgy=np.array(dft['Fpgy'])[edge:-1*edge-1]

print("Length of the filtered array datsets [Ut,Pgt]: ["+str(len(Ut))+","+str(len(Pgxt))+"]")



#plotting of  residual vel. magnitudes.

fig=plt.figure(figsize=(8,10),frameon=True)
plt.plot(Ns,color='blue',label='sim')
plt.plot(Nb,color='red',label='obs')
plt.ylabel('Residual velocity mag')
plt.legend()
plt.savefig(test_loc+'velmag.jpg',format='jpg')
plt.close(fig)



def pltvecttraj(llim,ulim,trajprt):
	Y=np.zeros(len(Usfil))
	Tib1=Tib[llim:ulim]
	fig=plt.figure(figsize=(8,10),frameon=True)
	plt.quiver(Tib1,Y[llim:ulim],Usfil[llim:ulim],Vsfil[llim:ulim],scale=2,scale_units='inches',width=0.002,color='blue',label='sim')
	plt.quiver(Tib1,Y[llim:ulim],Ubfil[llim:ulim],Vbfil[llim:ulim],scale=2,scale_units='inches',width=0.002,color='red',label='obs')
	plt.quiver(Tib1,Y[llim:ulim],Ut[llim:ulim],Vt[llim:ulim],scale=2,scale_units='inches',width=0.002,color='green',label='tidal_cur')
	# plt.quiver(Tib1,-Pgyt[llim:ulim],+Pgxt[llim:ulim],scale=2,scale_units='inches',color='magenta',label='-pgtid')
	plt.legend()
	plt.xlabel('Time')
	plt.ylabel('Residual velocity vectors')
	labels=Tib1[::25]
	plt.xticks(Tib1[::25], labels, rotation='vertical')

	plt.savefig(test_loc+'model-obs-tid'+trajprt+'.jpg',format='jpg')
	plt.close(fig)


# plotting of residual vel vectors. 
trajprtvec=['1','2','3']
for i in trajprtvec:
	if i=='1':
		llim=0;ulim=250
		pltvecttraj(llim,ulim,i)
	elif i=='3':
		llim=-250;ulim=-1
		pltvecttraj(llim,ulim,i)
	else:
		llim=2000;ulim=2250
		pltvecttraj(llim,ulim,i)


## copying ice drift vel and lon lat files from the main folder.

shutil.copy(loc+'ice_drift.jpg',test_loc+'ice_drift.jpg') 
shutil.copy(loc+'velplot.jpg',test_loc+'velplot.jpg') 
shutil.copy(loc+'Longitude.jpg',test_loc+'Longitude.jpg') 
shutil.copy(loc+'Latitude.jpg',test_loc+'Latitude.jpg') 
shutil.copy(loc+file,test_loc+file) 
shutil.copy(loc+'out.log',test_loc+'out.log') 