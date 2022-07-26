#!/bin/bash
#Params: configFile, mdinp, mdout, qmout, paramfilename, outputDir, num of focus points, samples per dimension, (maxdims)
python3 GridSlicer.py ../FF_GenOpt_conf/FF_GenOpt_otf_clap.cfg otf_clap/otf.inp /dev/shm/otf_clap_slicingmdout.log otf_clap/otf.log otf_clap/temp_parameters.str otf_clap/ 50 21
