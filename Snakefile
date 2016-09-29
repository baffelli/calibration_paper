import os
import glob
import pyrat.fileutils.gpri_files as gpf
import pyrat
import pyrat.gpri_utils
from types import SimpleNamespace
import re
import json




subworkflow old_data:
    snakefile: './analyze_old_data.snake'
    configfile: './calibration_configuration_20160222.json'

subworkflow new_data:
    snakefile: './calibration_analysis.snake'
    configfile: './calibration_configuration_chutze.json'

rule all:
    input:
        old_data('outputs/img/HV_gain.pdf'),
#        old_data('outputs/img/HV_loss.pdf'),
        new_data('analysis.done'),
#        'doc/calibration_paper.pdf'










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


