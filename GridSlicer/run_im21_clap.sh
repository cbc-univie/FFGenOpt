#!/bin/bash
#Params: configFile, mdinp, mdout, qmout, paramfilename, outputDir, num of focus points, samples per dimension, (maxdims)
python3 GridSlicer.py ../FF_GenOpt_conf/FF_GenOpt_im21_clap.cfg im21_clap/im21.inp /dev/shm/im21_clap_slicingmdout.log im21_clap/im21_freq.log im21_clap/temp_parameters.str im21_clap/ 50 21
