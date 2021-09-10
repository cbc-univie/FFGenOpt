#!/bin/bash
#Params: configFile, mdinp, mdout, qmout, paramfilename, outputDir, num of focus points, samples per dimension, (maxdims)
python3 GridSlicer.py ../FF_GenOpt_conf/FF_GenOpt_oac_clap.cfg oac_clap/oac.inp /dev/shm/oac_clap_slicingmdout.log oac_clap/oac.log oac_clap/temp_parameters.str oac_clap/ 50 21
