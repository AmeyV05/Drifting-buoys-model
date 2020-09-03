# This is an additional script with runs through all the forcing for a particular buoy .
import numpy as np 
import logging
import os
import sys
import importlib
import main_script
import settings

# Bnum = sys.argv[1]

# function for editing the settings file for forcevec.
def forcevecedit(fvec):
	with open('settings.py', 'r') as f:
		lines = f.readlines()
	# print(lines[1])
	lines=lines[2:]
	with open("settings.py", "w") as f:
		lines.insert(0, "global forcevec\n")
		lines.insert(1, "forcevec="+fvec+"\n")
		f.write("".join(lines))
		# print(lines[0:2])




def main2(Bnum):
	fvec='[1,1,1,1,1,1,1,1,1,1,1,1]'
	forcevecedit(fvec)
	importlib.reload(settings)
	main_script.main(Bnum,0)
	fveclist=['1','1','1','1','1','1','1','1','1','1','1','1']
	for i in range(len(fveclist)):
		fveclist=['1','1','1','1','1','1','1','1','1','1','1','1']
		if (i%2==0 and i!=0):
			fveclist[i]='0'
			fveclist[i+1]='0'
			fvecn='['+','.join(fveclist)+']'
			print(fvecn)
			forcevecedit(fvecn)
			importlib.reload(settings)
			main_script.main(Bnum,0)
		elif (i==0):
			fveclist[i]='0'
			fvecn='['+','.join(fveclist)+']'
			print(fvecn)
			forcevecedit(fvecn)
			importlib.reload(settings)
			main_script.main(Bnum,0)
			
	fvec='[1,1,1,1,1,1,1,1,1,1,1,1]'
	forcevecedit(fvec)

main2(sys.argv[1])