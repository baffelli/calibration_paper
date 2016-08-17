import os
import glob
import pyrat.fileutils.gpri_files as gpf
import pyrat
import pyrat.gpri_utils
from types import SimpleNamespace
import re
import json






#def reflector_list_from_config():
#	reflectors = config['list_of_reflectors']
#	azidx = [ref[0][1] for ref in reflectors]
#	azidx_dec = [int(ref[0][1]/int(config['range_compression']['dec'])) for ref in reflectors]
#	ridx = [ref[0][0] for ref in reflectors]
#	undec_list = [ (r, az) for r, az in zip(ridx, azidx)]
#	dec_list = [ (r, az) for r, az in zip(ridx, azidx_dec)]
#	return undec_list, dec_list
#Heading
config['geocoding']= {}
config['geocoding']['DEM'] = 'geo/swissALTI3D_2016_Clip.dem'
config['geocoding']['ref_mli_par'] = 'slc/20140321_123707_AAAl.mli.par'
config['geocoding']['scan_heading'] = -69
config['geocoding']['table_name'] = 'Chutze'
include: '../Code/geocoding.snake'
#Set the working directory to be the current directory
workdir: './'



#Set the rule for geocoding


rule all:
    input:
        "slc/20140321_123707_AAAl.mli_gc.tif",
        "geo/Chutze.inc.tif",
        "geo/Chutze.sim_sar.tif",
        'geo/Chutze.ls_map.tif',
        'geo/Chutze.sh_map_gc.tif',
        'geo/Chutze.sim_sar_fgc.tif',
        'geo/Chutze.inc.tif',
        'geo/Chutze.u.tif',


rule cleanup:
    input:
    output:
    shell:
        """
            rm geo/Chutze*
            rm slc/*gc*
        """

rule mli:
    input:
        slc = "{slcname}.slc",
        slc_par = "{slcname}.slc.par",
    output:
        mli = "Chutze.mli",
        mli_par = "Chutze.mli.par",
    shell:
        "multi_look {input.slc} {input.slc_par} {output.mli} {output.mli_par} 2 4 0.5 0.35"



#rule to_bmp:
#    input:
#        mli = "{slcname}.mli",
#        mli_par = "{slcname}.mli.par",
#    output:
#        bmp = '{slcname}.mli.bmp'
#    run:
#        mli_par = gpf.par_to_dict(mli_par)
#        width = mli_par['width']
#        bmp_cmd = "raspwr {input.mli} " + str(width) + " - - - - - - - {output.bmp}"









