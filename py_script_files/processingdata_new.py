## This function processes the data from all the models.
# It particularly reads the buoy data and other corresponding model data 
# then interpolates the data to buoy locations and outputs those vectors to an excel file. 
import pandas as pd
import numpy as np
import scipy.io
from scipy.interpolate import interp2d
from scipy.interpolate import griddata
import pykrige.kriging_tools as kt
from pykrige.ok import OrdinaryKriging
import generalfunc as gf
import logging
import os


def processall(Bnum,indexing,FD,sdate):
    prefix="BUOY_"
    bname=prefix+Bnum

        ##Buoy data processing
    BD=FD['BD'] #dictionary where buoy_data is saved.
    Xib=BD['Xib'][sdate:]; Yib=BD['Yib'][sdate:]
    Uib=BD['Uib'][sdate:]; Vib=BD['Vib'][sdate:]
    Tib=BD['Tib'][sdate:]