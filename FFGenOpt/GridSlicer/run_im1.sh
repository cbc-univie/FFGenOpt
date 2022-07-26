#!/bin/bash
#Params: configFile, mdinp, mdout, qmout, paramfilename, outputDir, num of focus points, samples per dimension, (maxdims)
python3 GridSlicer.py ../FF_GenOpt_conf/FF_GenOpt_im1.cfg im1/im1.inp /dev/shm/im1_slicingmdout.log im1/im1.log im1/temp_parameters.str im1/ 50 21
