#!/bin/bash
#Params: configFile, mdinp, mdout, qmout, paramfilename, outputDir, num of focus points, samples per dimension, (maxdims)
python3 GridSlicer.py ../FF_GenOpt_conf/FF_GenOpt_oac.cfg oac/oac.inp /dev/shm/oacslicingmdout.log oac/oac.log oac/temp_parameters.str oac/ 50 21
