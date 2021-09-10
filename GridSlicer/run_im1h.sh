#!/bin/bash
#Params: configFile, mdinp, mdout, qmout, paramfilename, outputDir, num of focus points, samples per dimension, (maxdims)
python3 GridSlicer.py ../FF_GenOpt_conf/FF_GenOpt_im1h.cfg im1h/im1h.inp /dev/shm/im1h_slicingmdout.log im1h/im1h.log im1h/temp_parameters.str im1h/ 50 21

