#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: David Krach 
         david.krach@mib.uni-stuttgart.de

Copyright 2024 David Krach, Felix Weinhardt

Permission is hereby granted, free of charge, to any person obtaining a 
copy of this software and associated documentation files (the “Software”), 
to deal in the Software without restriction, including without limitation 
the rights to use, copy, modify, merge, publish, distribute, sublicense, 
and/or sell copies of the Software, and to permit persons to whom the 
Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included 
in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS 
OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL 
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER 
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, 
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN 
THE SOFTWARE.

"""

### HEADER ------------------------------------------------------------------------

import numpy as np
import os, sys, glob, shutil, h5py 
import time
import csv
import argparse
import fileinput
import stat

import matplotlib.pyplot as plt
import matplotlib.image as img
from natsort import natsorted
from functools import partial
import multiprocessing

import localdrag as ld
from localdrag import *

### -------------------------------------------------------------------------------

verbose = False

argParser = argparse.ArgumentParser()
argParser.add_argument( "-dir",         "--workingDirectory",                                                     help="Working directory with .pgm files."                              )
argParser.add_argument( "-vs",          "--voxelsize",         type=float,                                        help="Defines the voxel size of the image in [m]."                     )
argParser.add_argument( "-height",      "--heightOfDomain",    type=float,                                        help="Defines the voxel height of the domain in [m]."                  )
argParser.add_argument( "-dumux",       "--dumux_path",        type=str,                                          help="Path to my dumux build directory, with executable."              )
argParser.add_argument( "-exec",        "--dumux_exec",        type=str,                                          help="Name of dumux executable."                                       )
argParser.add_argument( "-deltaP",      "--deltaPressure",     type=float, default=1.4e-3,                        help="Defines the pressure difference from inlet to outlet in [Pa]."   )
argParser.add_argument( "-get_phi",     "--get_porosity",                  default=False,  action="store_true",   help="True if porosities should be computed and saved."                )
argParser.add_argument( "-no_run",      "--no_run",                        default=False,  action="store_true",   help="True if files should be copied and simulation be run."           )
argParser.add_argument( "-no_postproc", "--no_postprocessing",             default=False,  action="store_true",   help="True permeability and porosity Data should be harvested -> csv." )
argParser.add_argument( "-no_lambda",   "--no_lambda_given",               default=False,  action="store_true",   help="Prefactor files are given."                                      )
argParser.add_argument( "-vtu",         "--vtuOutput",                     default=False,  action="store_true",   help="Results should be written in vtu files, default is True."        )
args = argParser.parse_args()

lambda_given       = not args.no_lambda_given
datadir            = f'{args.workingDirectory}'
voxelsize          = float(args.voxelsize)
height             = float(args.heightOfDomain)
dumuxpath          = str(args.dumux_path)
executable         = str(args.dumux_exec)
pressure           = float(args.deltaPressure)
generic_input_file = f'genericInput.input'

# This has to be adapted based on folder structure 
# output_file_name   = f'perm_{datadir.split("/")[1]}_{datadir.split("/")[2]}.csv'
output_file_name   = f'perm_single_prec.csv'

get_porosity       = args.get_porosity
copy_run           = not args.no_run
postprocessing     = not args.no_postprocessing
parallel           = False
vtuOutput          = args.vtuOutput

metadata = {}
metadata['lambda_given']       = lambda_given
metadata['datadir']            = datadir
metadata['voxelsize']          = voxelsize
metadata['height']             = height
metadata['pressure']           = pressure
metadata['dumuxpath']          = dumuxpath
metadata['executable']         = executable
metadata['generic_input_file'] = generic_input_file
metadata['output_file_name']   = output_file_name
metadata['get_porosity']       = get_porosity
metadata['copy_run']           = copy_run
metadata['postprocessing']     = postprocessing
metadata['verbose']            = verbose
metadata['parallel']           = parallel
metadata['vtuOutput']          = vtuOutput

ld.run_dumux.get_executeable(metadata)


filelist = glob.glob(f'{datadir}/*.pgm')
filelist = natsorted(filelist)

if get_porosity:
    ld.postprocessing_sweep.all_porosities(filelist, datadir)

if copy_run:
    print(f'Starting sequential computation.')
    for i in range(len(filelist)):
        # Run actual dumux simulation
        ld.run_dumux.run(metadata, filelist[i])
        
        print(f'Progress: {((i+1)/len(filelist) * 100):2.2f} % done, domain {i+1} / {len(filelist)}\n')

if postprocessing:
    ld.postprocessing_sweep.crawl_output_benchmark(filelist, metadata)



