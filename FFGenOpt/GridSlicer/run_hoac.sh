#!/bin/bash
#Params: configFile, mdinp, mdout, qmout, paramfilename, outputDir, num of focus points, samples per dimension, (maxdims)
python3 GridSlicer.py ../FF_GenOpt_conf/FF_GenOpt_hoac.cfg hoac/hoac.inp /dev/shm/hoacslicingmdout.log hoac/hoac.log hoac/temp_parameters.str hoac/ 50 21
