## 
import numpy as np
import os
import shutil
import pandas as pd
import generalfunc as gf
import settings


# reading the sim_data file from the folder

# buoy number input.

def readsimdata(Bnum,folname):
	bname='BUOY_'+Bnum
	loc='../../generated_data/'+bname+'/'+folname+'/'
	file='simdata_'+bname+'.xlsx'
	df=pd.read_excel(loc+file,'Model_Data')
	Xis=np.array(df['Xis']);Yis=np.array(df['Yis'])
	Xib=np.array(df['Xib']);Yib=np.array(df['Yib'])
	Uis=np.array(df['Uis']);Vis=np.array(df['Vis'])
	Uib=np.array(df['Uib']);Vib=np.array(df['Vib'])
	Xsfil=np.array(df['Xsfil'],dtype=float);
	nnan=~np.isnan(Xsfil.astype(float))
	Xsfil=Xsfil[nnan]
	Ysfil=np.array(df['Ysfil']);Ysfil=Ysfil[nnan]
	Xbfil=np.array(df['Xbfil']);Xbfil=Xbfil[nnan]
	Ybfil=np.array(df['Ybfil']);Ybfil=Ybfil[nnan]
	Usfil=np.array(df['Usfil']);Usfil=Usfil[nnan]
	Vsfil=np.array(df['Vsfil']);Vsfil=Vsfil[nnan]
	Ubfil=np.array(df['Ubfil']);Ubfil=Ubfil[nnan]
	Vbfil=np.array(df['Vbfil']);Vbfil=Vbfil[nnan]
	T=np.array(df['Date(GMT)']);hi=np.array(df['Ice thickness'])
	# residual creation
	numtaps=2*24*2+1
	fedge=int(numtaps/2)
	Xsres=Xis[fedge:-fedge]-Xsfil;Ysres=Yis[fedge:-fedge]-Ysfil
	Xbres=Xib[fedge:-fedge]-Xbfil;Ybres=Yib[fedge:-fedge]-Ybfil
	Usres=Uis[fedge:-fedge]-Usfil;Vsres=Vis[fedge:-fedge]-Vsfil
	Ubres=Uib[fedge:-fedge]-Ubfil;Vbres=Vib[fedge:-fedge]-Vbfil
	## Reading from another sheet.
	df=pd.read_excel(loc+file,'FT_Data')
	SD={'Xis':Xis,'Yis':Yis,'Xib':Xib,'Yib':Yib,
		'Xsfil':Xsfil,'Ysfil':Ysfil,'Xbfil':Xbfil,'Ybfil':Ybfil,
		'Xsres':Xsres,'Ysres':Ysres,'Xbres':Xbres,'Ybres':Ybres,
		'Uis':Uis,'Vis':Vis,'Uib':Uib,'Vib':Vib,
		'Usfil':Usfil,'Vsfil':Vsfil,'Ubfil':Ubfil,'Vbfil':Vbfil,
		'Usres':Usres,'Vsres':Vsres,'Ubres':Ubres,'Vbres':Vbres,
		'Xsam':np.array(df['Xsam']),'Ysam':np.array(df['Ysam']),
		'Xbam':np.array(df['Xbam']),'Ybam':np.array(df['Ybam']),
		'Xsph':np.array(df['Xsph']),'Ysph':np.array(df['Ysph']),
		'Xbph':np.array(df['Xbph']),'Ybph':np.array(df['Ybph']),
		'Usam':np.array(df['Usam']),'Vsam':np.array(df['Vsam']),
		'Ubam':np.array(df['Ubam']),'Vbam':np.array(df['Vbam']),
		'Usph':np.array(df['Usph']),'Vsph':np.array(df['Vsph']),
		'Ubph':np.array(df['Ubph']),'Vbph':np.array(df['Vbph']),'Tft':np.array(df['Tft']),
		'BAmpLon':np.array(df['BAmpLon'][0:6]),'SAmpLon':np.array(df['SAmpLon'][0:6]),
		'BPhLon':np.array(df['BPhLon'][0:6]),'SPhLon':np.array(df['SPhLon'][0:6]),
		'BAmpLat':np.array(df['BAmpLat'][0:6]),'SAmpLat':np.array(df['SAmpLat'][0:6]),
		'BPhLat':np.array(df['BPhLat'][0:6]),'SPhLat':np.array(df['SPhLat'][0:6]),
		'tidearg':np.array(df['tideargLon'][0:6],dtype=int),'corarg':np.array(df['corarglon'][0:2],dtype=int),
		'BAmpU':np.array(df['BAmpU'][0:6]),'SAmpU':np.array(df['SAmpU'][0:6]),
		'BPhU':np.array(df['BPhU'][0:6]),'SPhU':np.array(df['SPhU'][0:6]),
		'BAmpV':np.array(df['BAmpV'][0:6]),'SAmpV':np.array(df['SAmpV'][0:6]),
		'BPhV':np.array(df['BPhV'][0:6]),'SPhV':np.array(df['SPhV'][0:6]),
		'tidearg':np.array(df['tideargU'][0:6],dtype=int),'corarg':np.array(df['corargU'][0:2],dtype=int),
		'T':T,'hi':hi}
	return(SD)


def readallbuoy(folder):
    BD={}
    for f in os.listdir(folder):
        if f!='BUOY_06.csv':
            strfile=folder+'/'+f
            buoy=os.path.basename(strfile)
            bname=buoy.split('.')[0]
            xbname=bname.split('_')[1]+'_x'
            ybname=bname.split('_')[1]+'_y'
            Tbname=bname.split('_')[1]+'_t'
            D=pd.read_csv(strfile)
            Xib=np.array(D['Lon'])[96:]
            Yib=np.array(D['Lat'])[96:]
            Tib=np.array(D['Date(GMT)'])[96:]
            BD[xbname]= Xib; BD[ybname]=Yib
            BD[Tbname]=Tib
    return(BD)