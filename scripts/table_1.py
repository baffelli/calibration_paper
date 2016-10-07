import csv

import pyrat.fileutils.gpri_files as gpf
import pyrat.gpri_utils.calibration as cal
import pyrat.core.corefun as cf

def table_1(inputs, outputs, threads, config, params, wildcards):
    # list of reflectors
    refl_list = [ref for ref in params['ref'] if ref['type'] == "cubic" or ref['type'] == 'triangular']
    refl_list = sorted(refl_list, key=lambda x: x['ridx'])
    # Load matrix and convert to coherency
    slc = gpf.gammaDataset(inputs['slc_par'], inputs['slc'])
    # list of results
    res_list = []
    with open(outputs[0], 'w+') as of:
        tabwrite = csv.writer(of, delimiter=',')
        tabwrite.writerow(
            ["rsl", 'type', 'RCS', "ridx", "azidx"])
        for current_reflector in refl_list:
            RCS = cf.dB(cal.cr_rcs(current_reflector['side'], slc.radar_frequency[0], type=current_reflector['type']))
            row = [slc.r_vec[current_reflector['ridx']], current_reflector['type'], RCS, current_reflector['ridx'],
                   current_reflector['azidx'] // slc.azimuth_looks]
            res_list.append(row)
            tabwrite.writerow(row)
            # tabwrite.writerow(['average r_ph', 'weighted average r_ph'])
            # res_list = np.array(res_list)
            # r_ph_opt = np.mean(res_list[:,0])
            # r_ph_opt_w = np.average(res_list[:, 0], weights = 1/res_list[:,1])
            # tabwrite.writerow([r_ph_opt, r_ph_opt_w])


table_1(snakemake.input, snakemake.output, snakemake.threads, snakemake.config, snakemake.params, snakemake.wildcards)
