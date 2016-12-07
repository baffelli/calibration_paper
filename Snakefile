import os
import glob
import pyrat.fileutils.gpri_files as gpf
import pyrat
import pyrat.gpri_utils
from types import SimpleNamespace
import re
import json
from collections import namedtuple


def load_config(config_file):
    with open(config_file) as input_file:
        return json.load(input_file)


hongg_conf = load_config('./calibration_configuration_20160222.json')
chutzen_conf = load_config('./calibration_configuration_chutze.json')
#Define list of reflectors from json config file for the scene with dihedrals
list_of_reflectors_dihedral = [ref for ref in hongg_conf['list_of_reflectors'] if ref['type'] == 'dihedral'][0]
#Define list of reflectors from json config file
list_of_reflectors = chutzen_conf['list_of_reflectors']

configfile: './calibration_configuration_chutze.json'




subworkflow old_data:
    workdir: './processed'
    snakefile:  pyrat.rules['slc_to_calibrated_c']
    configfile: './calibration_configuration_20160222.json'

subworkflow new_data:
    workdir: './processed'
    snakefile:  pyrat.rules['slc_to_calibrated_c']
    configfile: './calibration_configuration_chutze.json'


rule all:
    input:
        expand('fig/figure_{r}.{ext}',r=range(1,14), ext=['png', 'pdf']),
        'tab/table_1.csv',
        'tab/table_2.csv',
        'tab/table_3.csv',
        'tab/RMS_polcal.csv',
        'tab/table_squint.csv',
        'doc/calibration_paper.pdf'







###############################################################################
#Cleanup figure
rule cleanup_figures:
    run:
        shell("trash ./fig/*.pdf")


###############################################################################
#Select stylesheet depending on figure type
def select_style(wildcards):
    if wildcards.ext == 'pdf':
        return 'paper_style.rc'
    else:
        return 'slides_style.rc'


###############################################################################
#Plot figure 1: Oversampled magnitude/phase response of a TCR
rule fig1:
    output:
        paper_fig = 'fig/figure_1.{ext}',
        pres_fig = expand('fig/figure_1_{n}.{{ext}}',n=range(6))
    input:
        HH = new_data('slc_chan/20160914_145059_AAAl.slc'),
        VV = new_data('slc_chan/20160914_145059_BBBl.slc'),
        HH_desq = new_data('slc_desq/20160914_145059_AAAl.slc'),
        VV_desq = new_data('slc_desq/20160914_145059_BBBl.slc'),
        HH_corr = new_data('slc_corr/20160914_145059_AAAl.slc'),
        VV_corr = new_data('slc_corr/20160914_145059_BBBl.slc'),
        style = select_style
    params:
        ridx = list_of_reflectors[1]['ridx'],
        azidx = list_of_reflectors[1]['azidx']
    script:
        'scripts/figure_1.py'

###############################################################################
#Plot figure 2/3: Oversampled magnitude/phase response of all TCR
#before and after the phase correction

#this serves to select the proper channel
def select_slc_for_rule_2(wildcards):
    proc_type = 'coreg' if int(wildcards.n) == 2 else 'corr'
    name_HH = "slc_{name}/20160914_145059_AAAl.slc".format(name=proc_type)
    name_VV = "slc_{name}/20160914_145059_BBBl.slc".format(name=proc_type)
    VV = new_data(name_VV)
    HH = new_data(name_HH)
    return HH, VV

rule fig2:
    output:
        'fig/figure_{n,(2|3)}.{ext}'
    input:
        VV = select_slc_for_rule_2,
        a = new_data("slc_corr/20160914_145059_AAAl.slc"),
        b = new_data("slc_coreg/20160914_145059_AAAl.slc"),
        c = new_data("slc_corr/20160914_145059_BBBl.slc"),
        d = new_data("slc_coreg/20160914_145059_BBBl.slc"),
        style = select_style
    params:
        reflectors = list_of_reflectors,
    script:
        'scripts/figure_2.py'


###############################################################################
#Plot figure 4: Gain in HV response after coregistration
rule fig4:
    input:
        C_HV_new = old_data("slc_coreg_common/20160224_130521_ABBl.slc"),
        C_HV_old = old_data("slc_coreg_common/20160224_105201_ABBl.slc"),
        style = select_style
    output:
        'fig/figure_4.{ext}'
    params:#position of dihedral
        ref = list_of_reflectors_dihedral
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
        style = select_style
    output:
        paper_fig = 'fig/figure_{n, (5|6)}.{ext}'
    params:
        ref = lambda wildcards: list_of_reflectors[1] if int(wildcards.n) == 5 else list_of_reflectors[-1]
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
        style = select_style
    params:
        ref = list_of_reflectors
    output:
        paper_fig = 'fig/figure_7.{ext}',
        pres_fig = 'fig/figure_7_full.{ext}'
    script:
        'scripts/figure_7.py'


###############################################################################
#Make figure 8: dependence of residuals with incidence angle
rule fig8:
    input:
        res = 'tab/table_3.csv',
        inc = new_data('geo/Chutzen.inc_fgc'),
        C_cal_par = new_data("cov_cal/20160914_145059_l.par"),
        style = select_style
    output:
        'fig/figure_8.{ext}'
    script:
        'scripts/figure_8.py'

###############################################################################
#Plot figure 9/10: HH/VV phase before and after removal of topographic contribution
#this serves to select the proper channel
def select_cov_for_rule_9(wildcards):
    proc_type = 'normal' if int(wildcards.n) == 9 else 'flat'
    HHVV = "cov_{name}/20160914_145059_l.c03".format(name=proc_type)
    HH = "cov_{name}/20160914_145059_l.c00".format(name=proc_type)
    VV = "cov_{name}/20160914_145059_l.c33".format(name=proc_type)
    HHVV = new_data(HHVV)
    HH = new_data(HH)
    VV = new_data(VV)
    return HHVV, HH, VV

rule fig9:
    input:
        style = select_style,
        aui = new_data("cov_normal/20160914_145059_l.par"),#dummy
        ali = new_data("cov_flat/20160914_145059_l.par"),
        C_cal_par = new_data("cov_cal/20160914_145059_l.par"),
        HHVV_phase = select_cov_for_rule_9,
    output:
        'fig/figure_{n, (9)|(10)}.{ext}'
    params:
        ref = list_of_reflectors
    script:
        'scripts/figure_9.py'

###############################################################################
#Plot figure 11: Slc with and without processing of squint
rule fig11:
    input:
        style = select_style,
        slc = new_data("slc_chan/20160914_145059_BBBl.slc_dec"),
        slc_par = new_data("slc_chan/20160914_145059_BBBl.slc_dec.par"),
        slc_desq = new_data("slc_desq/20160914_145059_BBBl.slc_dec"),
        slc_desq_par = new_data("slc_desq/20160914_145059_BBBl.slc_dec.par"),
        slc_corr = new_data("slc_corr/20160914_145059_BBBl.slc_dec"),
        slc_corr_par = new_data("slc_corr/20160914_145059_BBBl.slc_dec.par"),
    output:
        'fig/figure_11.{ext}'
    script:
        'scripts/figure_11.py'

###############################################################################
#Plot figure 12: Squint vs azimuth for VV and HH channel
def select_raw_for_rule_12(wildcards):
    chan_name = 'AAAl' if int(wildcards.n) == 12 else 'BBBl'
    chan = "raw_chan/20160914_145059_{chan}.raw".format(chan=chan_name)
    par = "raw_chan/20160914_145059_{chan}.raw_par".format(chan=chan_name)
    chan = new_data(chan)
    par = new_data(par)
    return chan, par

chan_string = new_data("raw_{}/20160914_145059_{}.raw")
rule fig12:
    input:
        HH = chan_string.format('chan', 'AAAl'),
        VV = chan_string.format('chan', 'BBBl'),
        HH_desq = chan_string.format('desq', 'AAAl'),
        VV_desq = chan_string.format('desq', 'BBBl'),
        style = select_style,
        slc_par = new_data("slc_chan/20160914_145059_BBBl.slc.par")
    params:
        ref = list_of_reflectors[1]
    output:
        paper_fig = 'fig/figure_12.{ext}',
        pres_fig = expand('fig/figure_12_{n}.{{ext}}',n=range(4))
    script:
        'scripts/figure_12.py'


###############################################################################
#Plot figure 13: H and V patterns
rule fig13:
    input:
        H_pat = './H_mainlobe_171_GHz.txt',
        V_pat = './V_mainlobe_171_GHz.txt',
        coreg_par = old_data("diff/20160224_105201_AAAl_BBBl.off_par"),
        slc_par = old_data("slc_coreg_common/20160224_105201_BBBl.slc"),
        style = select_style,
    output:
        'fig/figure_13.{ext}'
    script:
        'scripts/figure_13.py'
###############################################################################
#Make table 1: Location and RCS for all reflectors
rule table1:
    input:
        LUT = new_data('geo/Chutzen.gpri_to_dem'),
        dem_seg_par = new_data('geo/Chutzen.dem_seg.par'),
        slc_par = new_data("slc_corr/20160914_145059_BBBl.slc_dec.par"),
        slc = new_data("slc_corr/20160914_145059_BBBl.slc_dec"),
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
#Make table 4: Squint parameters:
rule table4:
    input:
        HH = chan_string.format('chan', 'AAAl'),
        VV = chan_string.format('chan', 'BBBl'),
        HH_desq = chan_string.format('desq', 'AAAl'),
        VV_desq = chan_string.format('desq', 'BBBl'),
        slc_par = new_data("slc_chan/20160914_145059_BBBl.slc.par")
    params:
        ref = list_of_reflectors
    output:
        'tab/table_squint.csv'
    script:
        'scripts/table_squint.py'



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
### Whenever the library is synchronized, cleanup all the {online} tags and
#the other unecessary things
rule cleanup_bibtex:
    input:
        bib = '../Texts/library.bib'
    output:
        'doc/library.bib'
    run:
        import re
        url_re = re.compile(r"url = .+,")
        abstract_re = re.compile(r"abstract = .+,")
        month_re = re.compile(r"(?P<tag>month)\s=\s\{{(?P<month>\S+)\}}")
        with open(input.bib) as infile, open(output[0],'w') as outfile:
            for line in infile:
                new_line = re.sub(abstract_re, "", re.sub(url_re, "",re.sub(month_re,r"\1 = \{ \2 \},",line)))
                print(new_line)
                outfile.write(new_line)


###############################################################################
### Remake the paper when sections or figures change
rule paper:
    input:
        sections = glob.glob('doc/sections/*.tex'),
        main_paper = 'doc/calibration_paper.tex',
        figures = expand('fig/figure_{n}.pdf', n=range(1,12)),
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


