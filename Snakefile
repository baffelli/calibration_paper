import os
import glob
import pyrat.fileutils.gpri_files as gpf
import pyrat
import pyrat.gpri_utils
from types import SimpleNamespace
import re
import json





subworkflow old_data:
    workdir:    './processed'
    snakefile:  pyrat.rules['slc_to_calibrated_c']
    configfile: './calibration_configuration_20160222.json'

subworkflow new_data:
    workdir:    './processed'
    snakefile:  pyrat.rules['slc_to_calibrated_c']
    configfile: './calibration_configuration_chutze.json'

rule all:
    input:
        'fig/figure_1.pdf'
#        old_data('outputs/img/HV_gain.pdf'),
##        old_data('outputs/img/HV_loss.pdf'),
#        new_data('analysis.done'),
##        'doc/calibration_paper.pdf'



###############################################################################

rule fig1:
    output:
        'fig/figure_1.pdf'
    input:
        HH = new_data('slc_chan/20160914_145059_AAAl.slc'),
        VV = new_data('slc_chan/20160914_145059_BBBl.slc'),
        HH_desq = new_data('slc_desq/20160914_145059_AAAl.slc'),
        VV_desq = new_data('slc_desq/20160914_145059_BBBl.slc'),
        HH_corr = new_data('slc_corr/20160914_145059_AAAl.slc'),
        VV_corr = new_data('slc_corr/20160914_145059_BBBl.slc'),
        style = 'paper_style.rc'
    params:
        ridx = 1331,
        azidx = 1490
    script:
        'scripts/figure_1.py'








###############################################################################
### Remake the paper when sections or figures change
rule paper:
    input:
        sections = glob.glob('doc/sections/*.tex'),
        main_paper = 'doc/calibration_paper.tex',
        figures = glob.glob('outputs/img/*.pdf')
    output:
        paper_pdf = 'doc/calibration_paper.pdf'
    shell:
        """
            cd doc
            pdflatex $(basename {input.main_paper})
            bibtex $(basename {input.main_paper} .tex)
            pdflatex $(basename {input.main_paper})
        """


