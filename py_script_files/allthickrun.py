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
def forcevecedit(fvec,hi):
	with open('settings.py', 'r') as f:
		lines = f.readlines()
	# print(lines[1])
	lines=lines[4:]
	with open("settings.py", "w") as f:
		lines.insert(0, "global forcevec\n")
		lines.insert(1, "forcevec="+fvec+"\n")
		lines.insert(2, "global h\n")
		lines.insert(3, "h="+str(hi)+"\n")
		f.write("".join(lines))
		# print(lines[0:2])


def main2(Bnum):
	# fvec='[1,1,1,1,1,1,1,1,1,1,1,1,0,0]'
	# hi=1.
	# forcevecedit(fvec,hi)
	# importlib.reload(settings)
	# 0 is for preprocess off and 1 is for on.
	fveclist=['1','1','1','1','1','1','1','1','1','1','1','1','0','0']
	hivec=[0.5,1.0,1.5,2.0,3.0]  ## h0- 0.5, h1=1 h2-1.5, h3-2.0, h4-3.0
	for i in range(len(hivec)):
		hi=hivec[i]
		fveclist=['1','1','1','1','1','1','1','1','1','1','1','1','0','0']
		fveclist[1]=str(i)
		fvecn='['+','.join(fveclist)+']'
		print(fvecn)
		forcevecedit(fvecn,hi)
		importlib.reload(settings)
		# 0 is for preprocess off and 1 is for on.
		main_script.main(Bnum,'0')

	fvec='[1,1,1,1,1,1,1,1,1,1,1,1,0,0]'
	hi=1.
	forcevecedit(fvec,hi)

main2(sys.argv[1])