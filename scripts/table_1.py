import csv

import pyrat.fileutils.gpri_files as gpf
import pyrat.gpri_utils.calibration as cal
import pyrat.core.corefun as cf

import pyrat.geo.geofun as geo

def table_1(inputs, outputs, threads, config, params, wildcards):
    #Geocoding table to compute location of reflectors
    #Geocode
    LUT = geo.GeocodingTable(inputs.dem_seg_par, inputs.LUT)

    # list of reflectors
    refl_list = [ref for ref in params['ref'] if ref['type'] == "cubic" or ref['type'] == 'triangular']
    refl_list = sorted(refl_list, key=lambda x: x['ridx'])
    # Load matrix and convert to coherency
    slc = gpf.gammaDataset(inputs['slc_par'], inputs['slc'])
    print(slc._params)
    # list of results
    res_list = []
    with open(outputs[0], 'w+') as of:
        tabwrite = csv.writer(of, delimiter=',')
        tabwrite.writerow(
            ["name","rsl", 'type', 'RCS', "ridx", "azidx", 'easting', 'northing'])
        for current_reflector in refl_list:
            #compute rcs
            RCS = cf.dB(cal.cr_rcs(current_reflector['side'], slc.radar_frequency, type=current_reflector['type']))
            #compute indices
            azidx_dec = current_reflector['azidx'] // slc.GPRI_decimation_factor
            ridx =  current_reflector['ridx']
            #Compute geographical location
            geo_coord = LUT.dem_coord_to_geo_coord(LUT.radar_coord_to_dem_coord([ridx, azidx_dec]))
            row = [current_reflector.get('name', "N/A"), slc.r_vec[current_reflector['ridx']], current_reflector['type'], RCS, current_reflector['ridx'],
                   azidx_dec, geo_coord[0], geo_coord[1]]
            res_list.append(row)
            tabwrite.writerow(row)
            # tabwrite.writerow(['average r_ph', 'weighted average r_ph'])
            # res_list = np.array(res_list)
            # r_ph_opt = np.mean(res_list[:,0])
            # r_ph_opt_w = np.average(res_list[:, 0], weights = 1/res_list[:,1])
            # tabwrite.writerow([r_ph_opt, r_ph_opt_w])


table_1(snakemake.input, snakemake.output, snakemake.threads, snakemake.config, snakemake.params, snakemake.wildcards)
