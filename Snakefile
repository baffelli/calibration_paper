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
        'fig/figure_1.pdf',
        'fig/figure_2.pdf',
        'fig/figure_3.pdf'
#        old_data('outputs/img/HV_gain.pdf'),
##        old_data('outputs/img/HV_loss.pdf'),
#        new_data('analysis.done'),
##        'doc/calibration_paper.pdf'



#Define list of reflectors
list_of_reflectors = [
    [32, 5192, "t"],
    [1331, 1490, "t"],
    [830, 3503, "t"],
    [1032, 5703, "t"],
    [1051, 5617, "t"],
    [3518, 3545, "t"]
]


###############################################################################
#Cleanup figure
rule cleanup_figures:
    run:
        shell("trash ./fig/*.pdf")

###############################################################################
#Plot figure 1: Oversampled magnitude/phase response of a TCR
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
        ridx = list_of_reflectors[1][0],
        azidx = list_of_reflectors[1][1]
    script:
        'scripts/figure_1.py'

###############################################################################
#Plot figure 2/3: Oversampled magnitude/phase response of all TCR
#before and after the phase correction

#this serves to select the proper channel
def select_slc(wildcards):
    proc_type = 'coreg' if int(wildcards.n) == 2 else 'corr'
    VV = new_data('slc_{name}/20160914_145059_BBBl.slc'.format(name=proc_type))
    return VV

rule fig2:
    output:
        'fig/figure_{n,(2|3)}.pdf'
    input:
        VV = select_slc,
        style = 'paper_style.rc'
    params:
        reflectors = list_of_reflectors,
    script:
        'scripts/figure_2.py'





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


