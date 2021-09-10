#!/bin/bash
#Params: configFile, mdinp, mdout, qmout, paramfilename, outputDir, num of focus points, samples per dimension, (maxdims)
python3 GridSlicer.py ../FF_GenOpt_conf/FF_GenOpt_im1h_clap.cfg im1h_clap/im1h.inp /dev/shm/slicingmdout.log im1h_clap/im1h_init.log im1h_clap/temp_parameters.str im1h_clap/ 50 21

