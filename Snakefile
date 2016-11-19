import os
import glob
import pyrat.fileutils.gpri_files as gpf
import pyrat
import pyrat.gpri_utils
from types import SimpleNamespace
import re
import json


configfile: './calibration_configuration_chutze.json'

#Define list of reflectors from json config file
list_of_reflectors = config['list_of_reflectors']


subworkflow old_data:
    workdir: './processed'
    snakefile:  pyrat.rules['slc_to_calibrated_c']
    configfile: './calibration_configuration_20160222.json'

subworkflow new_data:
    workdir: './processed'
    snakefile:  pyrat.rules['slc_to_calibrated_c']
    configfile: './calibration_configuration_chutze.json'

#subworkflow geo:
#    workdir: './processed'
#    snakefile:  pyrat.rules['geocoding']
#    configfile: './calibration_configuration_chutze.json'



rule all:
    input:
        new_data('geo/Chutzen.mli_gc.tif'),
        'fig/figure_1.pdf',
        'fig/figure_2.pdf',
        'fig/figure_3.pdf',
        'fig/figure_4.pdf',
        'fig/figure_5.pdf',
        'fig/figure_6.pdf',
        'fig/figure_7.pdf',
        'fig/figure_8.pdf',
        'tab/table_1.csv',
        'tab/table_2.csv',
        'tab/table_3.csv',
        'tab/RMS_polcal.csv',
        'doc/calibration_paper.pdf'








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
        ridx = list_of_reflectors[1]['ridx'],
        azidx = list_of_reflectors[1]['azidx']
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
#Plot figure 4: Gain in HV response after coregistration
rule fig4:
    input:
        C_HV_new = old_data("slc_corr/20160224_130521_ABBl.mli_dec"),
        C_HV_old = old_data("slc_corr/20160224_105201_ABBl.mli_dec"),
        style = 'paper_style.rc'
    output:
        'fig/figure_4.pdf'
    params:#position of dihedral
        ridx = 720,
        azidx = 168
    script:
        'scripts/figure_4.py'

###############################################################################
#Plot figure 5/6: Polarisation signatures for two reflectors
rule fig5:
    input:
        C = new_data(expand("cov_flat/20160914_145059_l.c{i}{j}",i=range(4),j=range(4))),
        C_par = new_data("cov_flat/20160914_145059_l.par"),
        C_cal = new_data(expand("cov_cal/20160914_145059_l.c{i}{j}",i=range(4),j=range(4))),
        C_cal_par = new_data("cov_cal/20160914_145059_l.par"),
        style = 'paper_style.rc'
    output:
        'fig/figure_{n, (5|6)}.pdf'
    params:
        ref = lambda wildcards: list_of_reflectors[1] if wildcards.n == 1 else list_of_reflectors[-1]
    script:
        'scripts/figure_5.py'

###############################################################################
#Plot figure 7: Calibrated pauli with location of reflectors
rule fig7:
    input:
        C_cal = new_data(expand("cov_cal/20160914_145059_l.c{i}{j}",i=range(4),j=range(4))),
        map = 'processed/geo/pk25krel_latest_Clip.tif',
        C_cal_par = new_data("cov_cal/20160914_145059_l.par"),
        LUT = new_data('geo/Chutzen.gpri_to_dem'),
        sh_map = new_data('geo/Chutzen.sh_map_gc'),
        dem_seg_par = new_data('geo/Chutzen.dem_seg.par'),
        style = 'paper_style.rc'
    params:
        ref = list_of_reflectors
    output:
        'fig/figure_7.pdf'
    script:
        'scripts/figure_7.py'


###############################################################################
#Make table 1: Location and RCS for all reflectors
rule table1:
    input:
        slc_par = new_data("slc_coreg/20160914_145059_BBBl.slc.par"),
        slc = new_data("slc_coreg/20160914_145059_BBBl.slc"),
    output:
        'tab/table_1.csv'
    params:
        ref = list_of_reflectors
    script:
        'scripts/table_1.py'

###############################################################################
#Make table 2: Phase center estimate for all reflectors:
rule table2:
    input:
        slc_VV = new_data("slc_coreg/20160914_145059_BBBl.slc"),
        slc_VV_par = new_data("slc_coreg/20160914_145059_BBBl.slc.par"),
        slc_HH = new_data("slc_coreg/20160914_145059_AAAl.slc"),
        slc_HH_par = new_data("slc_coreg/20160914_145059_AAAl.slc.par"),
    output:
        'tab/table_2.csv'
    params:
        ref = list_of_reflectors
    script:
        'scripts/table_2.py'

###############################################################################
#Make table 2: Calibration residuals:
rule table3:
    input:
        C = new_data(expand("cov_cal/20160914_145059_l.c{i}{j}",i=range(4),j=range(4))),
        C_par = new_data("cov_cal/20160914_145059_l.par"),
    output:
        'tab/table_3.csv'
    params:
        ref = list_of_reflectors
    script:
        'scripts/table_3.py'


###############################################################################
#Make table 3: Residuals RMS:
rule RMS_residual:
    input:
        res = 'tab/table_3.csv'
    output:
        'tab/RMS_polcal.csv'
    script:
        'scripts/RMS_polcal.py'

###############################################################################
#Make figure 8: dependence of residuals with incidence angle
rule fig8:
    input:
        res = 'tab/table_3.csv',
        inc = new_data('geo/Chutzen.inc_fgc'),
        C_cal_par = new_data("cov_cal/20160914_145059_l.par"),
        style = 'paper_style.rc'
    output:
        'fig/figure_8.pdf'
    script:
        'scripts/figure_8.py'



###############################################################################
### Whenever the library is synchronized, cleanup all the {online} tags and
#the other unecessary things
rule cleanup_bibtex:
    input:
        bib = '../Texts/library.bib'
    output:
        'doc/library.bib'
    run:
        import re
        url_re = re.compile(r"url\s=\s{\S+},")
        month_re = re.compile(r"(?P<tag>month)\s=\s\{{(?P<month>\S+)\}}")
        with open(input.bib) as infile, open(output[0],'w') as outfile:
            for line in infile:
               outfile.write(re.sub(url_re, "",re.sub(month_re,r"\1 = \{\2\} ,",line)))


###############################################################################
### Remake the paper when sections or figures change
rule paper:
    input:
        sections = glob.glob('doc/sections/*.tex'),
        main_paper = 'doc/calibration_paper.tex',
        figures = expand('fig/figure_{n}.pdf', n=range(1,9)),
        library = 'doc/library.bib'
    output:
        paper_pdf = 'doc/calibration_paper.pdf'
    shell:
        """
            cd doc
            pdflatex $(basename {input.main_paper})
            bibtex $(basename {input.main_paper} .tex)
            pdflatex $(basename {input.main_paper})
        """


