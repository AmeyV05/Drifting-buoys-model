# main script file.

import numpy as np
import body_program as BP

# buoy number decides the indexing
# '03' : 5856

def bindexing(Bnum,switcher):
 return switcher.get(Bnum,"Invalide Buoy number")

#
def readingalldata():


	return()
switcher={
  '02':3152,
  '03':5952,
  '09':3160,
  '07':1192,
  '12':3024,
  '13':2524,
  '14':2136,
  '16':5952
  }
#Bnum="15"
print("Do you wish to run the code for all the buoys or a particular buoy? ")
print("Type 1, if you want to run for all the buoys else any other number for a particular buoy data.")
count=int(input())
Cor=int(input("Do you want Coriolis force? If yes, type 1 else any other number: "))

if (count==1):
 for i in switcher:
  Bnum=i
  indexing=switcher[i]-96
  BP.main(Bnum,indexing,Cor)
else:
 print("Please input the Buoy number from the list: [02, 03, 07, 09, 12, 13, 14, 16]. And press enter.")
 Bnum=input()
 indexing=switcher.get(Bnum,"Invalid Buoy number")-96
 BP.main(Bnum,indexing,Cor)


