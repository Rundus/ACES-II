# --- data_viewer.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Set of functions useful for quickly looking at the contents of a single CDF file


# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"
__date__ = "2022-08-22"
__version__ = "1.0.0"
# -------------------


from cdflib import cdfread



# ------
# Inputs
# ------

from glob import glob

#Toggle
wFile = 3

inputFile = glob(r'D:\Data\ACESII\L0\high\*.cdf')
print(inputFile)
input_path = inputFile[wFile]
print(input_path)


cdfFile = cdfread.CDF(input_path)
info = cdfFile.cdf_info()
zvars = info['zVariables']

# toggles
which_zvar = 0

# Functions
def single_var_inquire(variable): # just pass in something from the zvars list
    print(cdfFile.varinq(variable))
def single_var_attributes(variable):
    print(cdfFile.varattsget(variable,expand=True))
    return cdfFile.varattsget(variable,expand=True)
def single_var_data(variable):
    print(cdfFile.varget(variable,expand=True))
    return cdfFile.varget(variable,expand=True)
def all_var_attributes():
    for i, zvar in enumerate(zvars):
        print(f'----- {zvar} VARIABLE PARAMETERS-----')
        varinq = cdfFile.varinq(zvar)
        for thing in varinq:
            print(f' {thing} :  {varinq[thing]}')

        print('\n')
        print(f'----- {zvar} Attribute PARAMETERS-----')
        varatt = cdfFile.varattsget(zvar)
        for thing in varatt:
            print(f' {thing} :  {varatt[thing]}')
        print('\n' * 2)
def global_file_attributes():
    print(cdfFile.globalattsget())
def all_var_attributes_mod():
    for i,zvar in enumerate(zvars):
        print(f'----- {zvar} VARIABLE PARAMETERS-----',end='')
        varinq = cdfFile.varinq(zvar)

        print(f'----- {zvar} Attribute PARAMETERS-----', end='')
        varatt = cdfFile.varattsget(zvar)

# EXECUTE CODE HERE
from ACESII_code.class_var_func import setupPYCDF
setupPYCDF()
from spacepy import pycdf

count_interval = single_var_data('Count_Interval')['Data']

print(len(count_interval))


