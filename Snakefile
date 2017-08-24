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


hongg_conf = load_config('calibration_configuration_20160222.json')
chutzen_conf = load_config('calibration_configuration_chutze.json')
#Define list of reflectors from json config file for the scene with dihedrals
list_of_reflectors_dihedral = [ref for ref in hongg_conf['list_of_reflectors'] if ref['type'] == 'dihedral'][0]
#Define list of reflectors from json config file
list_of_reflectors = chutzen_conf['list_of_reflectors']


workdir: '.'
configfile: './calibration_configuration_chutze.json'



#Process data of antenna shift experiment
subworkflow old_data:
    workdir: './data'
    snakefile:  pyrat.rules['slc_to_calibrated_c']
    configfile: 'calibration_configuration_20160222.json'

#Process data of Chutzen calibration campaign
subworkflow new_data:
    workdir: './data'
    snakefile:  pyrat.rules['slc_to_calibrated_c']
    configfile: 'calibration_configuration_chutze.json'

#Recreate pdf illustrations
subworkflow illustrations:
    workdir: './drawings'
    snakefile: './drawings/Snakefile'




rule all:
    input:
        #figures from data
        expand('fig/{ext}/figure_{r}.{ext}',r=range(1,13), ext=[ 'pdf']),
        #expand('fig/{ext}/figure_17.{ext}', ext='pdf'),
        #tables
        'tab/table_1.csv',
        'tab/table_2.csv',
        'tab/table_3.csv',
        'tab/RMS_polcal.csv',
        'tab/table_squint.csv',
        'doc/calibration_paper.pdf',
        #figures
        illustrations('pdf/squint_correction.pdf'),
        illustrations('pdf/squint_correction_interpolation.pdf'),
        illustrations('pdf/antenna_squint.pdf'),
        illustrations('pdf/real_aperture_signal_model_geometry.pdf'),
        illustrations('pdf/antenna_offset.pdf'),
        illustrations('pdf/kapri_antenna_arrangement.pdf')

#Figures for presentation
rule presentation:
    input:
        expand('fig/{ext}/figure_{r}.{ext}',r=range(1,17), ext=['png']),
        expand('fig/{ext}/figure_signature_{r}.{ext}',r=range(0,6), ext=['png']),



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
#Utility rule to select location of TCR
rule select_TCR_location:
    input:
        slc = new_data("slc_desq/20160914_145059_AAAl.slc_dec"),
        slc_par = new_data("slc_desq/20160914_145059_AAAl.slc_dec.par"),
        map = 'data/geo/pk25krel_latest_Clip.tif',
        C_cal_par = new_data("cov_cal/20160914_145059_l.par"),
        LUT = new_data('geo/Chutzen.gpri_to_dem'),
        sh_map = new_data('geo/Chutzen.sh_map_gc'),
        dem_seg_par = new_data('geo/Chutzen.dem_seg.par'),
        LUT1 = new_data('geo/Chutzen.dem_to_gpri'),
        ref_mli_par = new_data('geo/Chutzen.mli.par'),
        dem_seg = new_data('geo/Chutzen.dem_seg.tif'),
    params:
        ref = list_of_reflectors
#    output:
#        paper_fig = 'fig/{ext}/figure_7.{ext}',
#        pres_fig = 'fig/{ext}/figure_7_full.{ext}'
    script:
        'scripts/select_reflectors.py'



###############################################################################
#Plot figure 1: Oversampled magnitude/phase response of a TCR
rule fig1:
    output:
        paper_fig = 'fig/{ext}/figure_1.{ext}',
        pres_fig = expand('fig/{{ext}}/figure_1_{n}.{{ext}}',n=range(6))
    input:
        HH = new_data('slc_chan/20160914_145059_AAAl.slc'),
        VV = new_data('slc_chan/20160914_145059_BBBl.slc'),
        HH_desq = new_data('slc_desq/20160914_145059_AAAl.slc'),
        VV_desq = new_data('slc_desq/20160914_145059_BBBl.slc'),
        HH_corr = new_data('slc_corr/20160914_145059_AAAl.slc'),
        VV_corr = new_data('slc_corr/20160914_145059_BBBl.slc'),
        style = select_style,
        script = 'scripts/figure_1.py'
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
        'fig/{ext}/figure_{n,(2|3)}.{ext}'
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
        'fig/{ext}/figure_4.{ext}'
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
        paper_fig = 'fig/{ext}/figure_{n, (5|6)}.{ext}'
    params:
        ref = lambda wildcards: list_of_reflectors[1] if int(wildcards.n) == 5 else list_of_reflectors[-1],
        sw = config['calibration']['search_window'],
        aw = config['calibration']['averaging_window']
    script:
        'scripts/figure_5.py'

###############################################################################
#Plot figure 7: Calibrated pauli with location of reflectors
rule fig7:
    input:
        C_cal = new_data(expand("cov_cal/20160914_145059_l.c{i}{j}",i=range(4),j=range(4))),
        map = 'data/geo/pk25krel_latest_Clip.tif',
        C_cal_par = new_data("cov_cal/20160914_145059_l.par"),
        LUT = new_data('geo/Chutzen.gpri_to_dem'),
        LUT_inv = new_data('geo/Chutzen.dem_to_gpri'),
        sh_map = new_data('geo/Chutzen.sh_map_gc'),
        dem_seg_par = new_data('geo/Chutzen.dem_seg.par'),
        style = select_style,
        dem_seg = new_data('geo/Chutzen.dem_seg.tif'),
    params:
        ref = list_of_reflectors
    output:
        paper_fig = 'fig/{ext}/figure_7.{ext}',
        pres_fig = 'fig/{ext}/figure_7_full.{ext}',
        tif = 'fig/{ext}/gc_pauli.tif'
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
        'fig/{ext}/figure_8.{ext}'
    script:
        'scripts/figure_8.py'

###############################################################################
#Plot figure 9/10: HH/VV phase before and after removal of topographic contribution
#this serves to select the proper channel
def select_cov_for_rule_9(wildcards):
    proc_type = 'normal' if int(wildcards.n) == 9 else 'normal'
    HHVV = "cov_{name}/20160914_145059_l.c03".format(name=proc_type)
    HH = "cov_{name}/20160914_145059_l.c00".format(name=proc_type)
    VV = "cov_{name}/20160914_145059_l.c33".format(name=proc_type)
    data = {'HHVV':HHVV, 'HH':HH, 'VV':VV}
    return data

HH = "cov_{name}/20160914_145059_l.c00"
VV = "cov_{name}/20160914_145059_l.c33"
HHVV = "cov_{name}/20160914_145059_l.c03"


rule fig9:
    input:
        HH = new_data(HH.format(name='normal')),
        VV = new_data(VV.format(name='normal')),
        HHVV = new_data(HHVV.format(name='normal')),
        style = select_style,
        C_cal_par = new_data("cov_cal/20160914_145059_l.par"),
        theta = new_data('geo/Chutzen.lv_theta_fgc'),
        script = 'scripts/figure_9.py'
    output:
        fig_a ='fig/{ext}/figure_9.{ext}',
        fig_b = 'fig/{ext}/figure_9_b.{ext}'
    params:
        ref = list_of_reflectors
    script:
        'scripts/figure_9.py'

rule fig10:
    input:
        HH = new_data(HH.format(name='flat')),
        VV = new_data(VV.format(name='flat')),
        HHVV = new_data(HHVV.format(name='flat')),
        style = select_style,
        C_cal_par = new_data("cov_cal/20160914_145059_l.par"),
        theta = new_data('geo/Chutzen.lv_theta_fgc'),
        script = 'scripts/figure_9.py'
    output:
        fig_a ='fig/{ext}/figure_10.{ext}',
        fig_b = 'fig/{ext}/figure_10_b.{ext}'
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
        'fig/{ext}/figure_11.{ext}'
    script:
        'scripts/figure_11.py'

###############################################################################
#Plot figure 12: Squint vs azimuth for VV and HH channel
def select_raw_for_rule_12(wildcards):
    chan_name = 'AAAl' if int(wildcards.n) == 12 else 'BBBl'
    chan = "raw_chan/20160914_145059_{chan}.raw".format(chan=chan_name)
    par = "raw_chan/20160914_145059_{chan}.raw_par".format(chan=chan_name)
    data = {'raw':chan,'raw_par':par}
    chan = new_data(chan)
    par = new_data(par)
    return chan, par

chan_string = new_data("raw_{}/20160914_145059_{}.raw")
rule fig12:
    input:
        HH = new_data("raw_chan/20160914_145059_AAAl.raw"),
        VV = new_data("raw_chan/20160914_145059_BBBl.raw"),
        HH_desq = new_data("raw_desq/20160914_145059_AAAl.raw"),
        VV_desq = new_data("raw_desq/20160914_145059_BBBl.raw"),
        style = select_style,
        slc_par = new_data("slc_chan/20160914_145059_BBBl.slc.par"),
        script = 'scripts/figure_12.py'
    params:
        ref = list_of_reflectors[1]
    output:
        paper_fig = 'fig/{ext}/figure_12.{ext}',
        pres_fig = expand('fig/{{ext}}/figure_12_{n}.{{ext}}',n=range(4))
    script:
        'scripts/figure_12.py'


###############################################################################
#Plot figure 13: H and V patterns
rule fig13:
    input:
        H_pat = 'H_mainlobe_171_GHz.txt',
        V_pat = 'V_mainlobe_171_GHz.txt',
        coreg_par = old_data("diff/20160224_105201_AAAl_BBBl.off_par"),
        slc_par = old_data("slc_coreg_common/20160224_105201_BBBl.slc"),
        style = select_style,
    output:
        'fig/{ext}/figure_13.{ext}'
    script:
        'scripts/figure_13.py'

###############################################################################
#Plot figure 14: H and alpha before and after calibration
#this serves to select the proper channel
def select_c_for_rule_14(wildcards):
    proc_type = 'flat' if int(wildcards.n) == 14 else 'cal'
    C_cal = new_data(expand("cov_{type}/20160914_145059_l.c{{i}}{{j}}".format(type=proc_type),i=range(4),j=range(4))),
    return C_cal[0]
rule fig14:
    input:
        C = select_c_for_rule_14,
        map = 'data/geo/pk25krel_latest_Clip.tif',
        C_cal_par = new_data("cov_cal/20160914_145059_l.par"),
        LUT = new_data('geo/Chutzen.gpri_to_dem'),
        sh_map = new_data('geo/Chutzen.sh_map_gc'),
        dem_seg_par = new_data('geo/Chutzen.dem_seg.par'),
        style = select_style
    params:
        ref = list_of_reflectors
    wildcard_constraints:
        n = '(14)|(15)'
    output:
        rgb_fig = 'fig/{ext}/figure_{n}.{ext}',
        H_fig = 'fig/{ext}/figure_{n}_H.{ext}',
        alpha_fig = 'fig/{ext}/figure_{n}_alpha.{ext}',
    script:
        'scripts/figure_14.py'


###############################################################################
#Plot figure: Polarisation signatures for all reflectors
rule fig_ref_n:
    input:
        C = new_data(expand("cov_flat/20160914_145059_l.c{i}{j}",i=range(4),j=range(4))),
        C_par = new_data("cov_flat/20160914_145059_l.par"),
        C_cal = new_data(expand("cov_cal/20160914_145059_l.c{i}{j}",i=range(4),j=range(4))),
        C_cal_par = new_data("cov_cal/20160914_145059_l.par"),
        style = select_style
    output:
        paper_fig = 'fig/{ext}/figure_signature_{n}.{ext}'
    params:
        ref = lambda wildcards: list_of_reflectors[int(wildcards.n)]
    script:
        'scripts/figure_5.py'



###############################################################################
#Plot figure: polarisation isolation at the calibration reflector
rule fig16:
    input:
        C = new_data(expand("cov_flat/20160914_145059_l.c{i}{j}",i=range(4),j=range(4))),
        C_par = new_data("cov_flat/20160914_145059_l.par"),
        style = select_style
    output:
        paper_fig = 'fig/{ext}/figure_16.{ext}'
    params:
        ref = lambda wildcards: list_of_reflectors[int(config['calibration']['reflector_index'])]
    script:
        'scripts/figure_16.py'

###############################################################################
#Plot figure: phase center location as a function of frequency
rule fig17:
    input:
        slc_VV = new_data("slc_coreg/20160914_145059_BBBl.slc"),
        slc_VV_par = new_data("slc_coreg/20160914_145059_BBBl.slc.par"),
        slc_HH = new_data("slc_coreg/20160914_145059_AAAl.slc"),
        slc_HH_par = new_data("slc_coreg/20160914_145059_AAAl.slc.par"),
        style = select_style
    output:
        paper_fig = 'fig/{ext}/figure_17.{ext}'
    params:
        ref = lambda wildcards: list_of_reflectors
    script:
        'scripts/figure_17.py'


###############################################################################
#Make table 1: Location and RCS for all reflectors
rule table1:
    input:
        LUT = new_data('geo/Chutzen.gpri_to_dem'),
        LUT_inv = new_data('geo/Chutzen.dem_to_gpri'),
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
        ref = list_of_reflectors,
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
        ref = list_of_reflectors,
        sw = config['calibration']['search_window'],
        aw = config['calibration']['averaging_window']
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
#Cleanup biblography
rule clean_bib:
    shell:
        """
        cd doc
        rm  *.bbl
        rm *.aux
        rm *.blg
        """

rule pull_bib:
    shell:
        """
        cd doc/biblography
        git pull
        """

################################################################################
#### Whenever the library is synchronized, cleanup all the {online} tags and
##the other unecessary things
rule cleanup_bibtex:
    input:
        bib = 'doc/library.bib'
    output:
        'doc/library_clean.bib'
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
        #figures = expand('fig/pdf/figure_{n}.pdf', n=range(1,12)),
        tables = glob.glob('tab/*.csv'),
        figures = glob.glob('fig/*/*.pdf'),
        drawings = glob.glob('drawings/pdf/*.pdf'),
        library = 'doc/library_clean.bib'
    output:
        paper_pdf = 'doc/calibration_paper.pdf'
    shell:
        """
            cd doc
            pdflatex $(basename {input.main_paper})
            bibtex $(basename {input.main_paper} .tex)
            pdflatex $(basename {input.main_paper})
        """


