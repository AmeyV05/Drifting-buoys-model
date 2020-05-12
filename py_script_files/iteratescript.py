# a script for iterating over Enkf
import numpy as np 
import enkfmain as enkf 
import optCw as oCw
import generalfunc as gf
import shutil
import logging
import pandas as pd 
import openpyxl

s=enkf.enkfsettings()

def main(Bnum,indexing,Cor):
 # main program for Cw iteration. 
 prefix="BUOY_"
 bname=prefix+Bnum
 path = "../generated_data/"+bname
 global s
 M=28420
 pathmain=path+"/exp"+str(M)
 gf.mkdir_p(pathmain)
 # hvec=[0.3,0.5,1.0]
 # Cavec=[1.2e-04,0.5e-03,1.2e-03,2.1e-03]
 # Cwvec=[1.2e-04,3.5e-04,5.5e-04,6.5e-04,9.e-04,1.2e-03,3.5e-03,5.5e-03,6.5e-03,7.5e-03,8.5e-03,1.0e-02]
 # Cwvec=[1.2e-04,5.5e-04,1.2e-03,5.5e-03,6.5e-03,8.5e-03]
 Cwvec=[5.5e-04,5.5e-03]
 # Cwvec=[8.5e-03,1.0e-02]
 hvec=[0.5]
 Cavec=[1.2e-04]
 i=0
 iterno=1
 Jtensor=np.zeros((len(hvec),len(Cavec),len(Cwvec),iterno))
 for h in hvec:
  j=0
  pathi=pathmain+"/"+str(hvec[i])
  gf.mkdir_p(pathi)
  for ca in Cavec:
   pathj=pathi+"/"+str(Cavec[j])
   gf.mkdir_p(pathj)
   k=0
   JD={}
   writer = pd.ExcelWriter(pathj+'/Jac_data.xlsx', engine='xlsxwriter')
   dfcw = pd.DataFrame({'Cwvec': Cwvec})
   dfcw.to_excel(writer, startcol=0,startrow=0,index=False)
   workbook  = writer.book
   worksheet = writer.sheets['Sheet1']
   merge_format = workbook.add_format({
    'bold': 1,
    'border': 1,
    'align': 'center',
    'valign': 'vcenter',
     })
   workbook.close() 
   for cw in Cwvec:
    pathk=pathj+"/"+str(Cwvec[k])
    gf.mkdir_p(pathk)
    Jvec=[]
    for l in range(iterno):
     pathl=pathk+"/"+str(l+1)
     gf.mkdir_p(pathl)
     s=enkf.enkfsettings()
     s['h']=hvec[i] ; s['iCa']=Cavec[j] ; s['iCw']=Cwvec[k] 
     logging.info("Starting simulation with following values:")
     logging.info("Ice thickness: " + str(s['h']))
     logging.info("Ice-Air drag coeff: " + str(s['iCa']))
     logging.info("Ice-Water drag coeff: "+str(s['iCw']))
     logging.info("Running iteration number "+str(l+1)+".")
     logging.info("Model simulation started.")
     # [Xis,Yis,results,zeta_ktot,times,PD]=enkf.simulate(s,Bnum,indexing,Cor)
     # oCw.proEnkfdata(s,Xis,Yis,results,zeta_ktot,times,PD,pathl,l)  
     # logging.info("Computing cost function.")	
     # Jval=oCw.CFval(s,times,PD,results)
     # logging.info("Cost function value is: "+str(Jval))
     [resmat,wk1mat,wk2mat,Jval,times]=enkf.CFevalsim(s,Bnum,indexing,Cor)
     enkf.proEnkfCFdata(resmat,wk1mat,wk2mat,Jval,times,pathl)
     logging.info(str(l+1)+" Iteration Done")
     Jtensor[i,j,k,l]=Jval
     Jten=Jtensor[i,j,:,l];jval2excel(Jten,l,pathj)
    # Jtensor[i,j,k,:]=Jvec
    k+=1
   # JD={'Cwvec':Cwvec,'Jac1':Jtensor[i,j,:,0],'Jac2':Jtensor[i,j,:,1],
  	#   'Jac3':Jtensor[i,j,:,2],'Jac4':Jtensor[i,j,:,3],'Jac5':Jtensor[i,j,:,4]}
   # jval2excel(JD,pathj)
   logging.info("Excel file created with Jacobian values.")
   j+=1
  i+=1  	
 # srfile='/out.log'
 # shutil.copy(path+srfile,pathmain+srfile)
 # logging.info("Done")

def jval2excel(Jten,l,pathj):
 workbook = openpyxl.load_workbook(pathj+'/Jac_data.xlsx')
 writer = pd.ExcelWriter(pathj+'/Jac_data.xlsx', engine='openpyxl') 
 writer.book = workbook
 writer.sheets = dict((ws.title, ws) for ws in workbook.worksheets)
 dfjten= pd.DataFrame({str(l): Jten})
 dfjten.to_excel(writer,'Sheet1', startcol=l+1,startrow=0,index=False)
 writer.save()




  
