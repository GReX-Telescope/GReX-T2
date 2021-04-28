#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# dsahead python 3.7
from astropy.io import ascii
from astropy.io.ascii.core import InconsistentTableError
import numpy as np
import hdbscan
import json
import os.path
from dsautils import dsa_store
ds = dsa_store.DsaStore()
import dsautils.dsa_syslog as dsl
logger = dsl.DsaSyslogger()
logger.subsystem('software')
logger.app('T2')

# half second at heimdall time resolution (after march 18)
offset = 1907
downsample = 4


def parse_candsfile(candsfile):
    """ Takes standard MBHeimdall giants output and returns full table, classifier inputs and snr tables.
    (Can add cleaning here, eventually)
    """

    if os.path.exists(candsfile):
        logger.debug(f'Candsfile {candsfile} is path, so opening it')
        candsfile = open(candsfile, 'r').read()
    else:
        ncands = len(candsfile.split('\n'))-1
        logger.debug(f'Received {ncands} candidates')
#    candsfile = '\n'.join([line for line in candsfile.split('\n') if line.count(' ') == 7])
#    print(f'Received {ncands0} candidates, removed {ncands0-ncands} lines.')
    col_heimdall = ['snr', 'if', 'itime', 'mjds', 'ibox', 'idm', 'dm', 'ibeam']
    col_T2old = ['snr', 'if', 'itime', 'mjds', 'ibox', 'idm', 'dm', 'ibeam', 'cl', 'cntc', 'cntb']
    col_T2 = ['snr', 'if', 'itime', 'mjds', 'ibox', 'idm', 'dm', 'ibeam', 'cl', 'cntc', 'cntb', 'trigger']

    try:
        tab = ascii.read(candsfile, names=col_heimdall, guess=True, fast_reader=False, format='no_header')
        logger.debug('Read with heimdall columns')
    except InconsistentTableError:
        try:
            tab = ascii.read(candsfile, names=col_T2, guess=True, fast_reader=False, format='no_header')
            logger.debug('Read with T2 columns')
        except InconsistentTableError:
            try:
                tab = ascii.read(candsfile, names=col_T2old, guess=True, fast_reader=False, format='no_header')
                logger.debug('Read with old style T2 columns')
            except InconsistentTableError:
                logger.warning('Inconsistent table. Skipping...')              
                return ([], [], [])

    tab['ibeam'] = tab['ibeam'].astype(int)
    try:
        ret_time = ds.get_dict('/mon/snap/1/armed_mjd')['armed_mjd']+float(ds.get_dict('/mon/snap/1/utc_start')['utc_start'])*4.*8.192e-6/86400.
    except:
        ret_time = 55000.0
    tab['mjds'] = tab['mjds']/86400.+ret_time
    
#
#    snrs = tab['snr']
    # how to use ibeam?
   
#    return tab, data, snrs
    return tab

def cluster_data(tab, selectcols=['itime', 'idm', 'ibox', 'ibeam'], min_cluster_size=2, min_samples=5, metric='hamming', return_clusterer=False, allow_single_cluster=True):
    """ Take data from parse_candsfile and identify clusters via hamming metric.
    selectcols will take a subset of the standard MBHeimdall output
    """

    data = np.lib.recfunctions.structured_to_unstructured(tab[selectcols].as_array())  # ok for single dtype (int)
    try:
        clusterer = hdbscan.HDBSCAN(metric=metric, min_cluster_size=min_cluster_size,
                                    min_samples=min_samples, cluster_selection_method='eom',
                                    allow_single_cluster=allow_single_cluster).fit(data) 

        nclustered = np.max(clusterer.labels_ + 1) 
        nunclustered = len(np.where(clusterer.labels_ == -1)[0]) 
        cl = clusterer.labels_
    except ValueError:
        logger.info("Clustering did not run. Each point assigned to unique cluster.")
        cl = np.arange(len(data))
        nclustered = 0
        nunclustered = len(cl)

    logger.info(f'Found {nclustered} clustered and {nunclustered} unclustered rows')

    # hack assumes fixed columns
    bl = data[:, 3]
    cntb, cntc = np.zeros((len(data), 1), dtype=int), np.zeros((len(data), 1), dtype=int)
    ucl = np.unique(cl)

    for i in ucl:
        ww = np.where(i == cl)
        cntc[ww] = len(ww[0]) 
        ubl = np.unique(bl[ww])
        cntb[ww] = len(ubl) 

    # append useful metastats to original data
#    data_labeled = np.hstack((data, cl[:,None], cntb, cntc))
    # modifies tab in place
    tab['cl'] = cl.tolist()
    tab['cntc'] = cntc.flatten().tolist()
    tab['cntb'] = cntb.flatten().tolist()

    if return_clusterer:
#        return clusterer, data_labeled
        return clusterer
#    else:
#        return data_labeled


def get_peak(tab):
    """ Given labeled data, find max snr row per cluster
    Adds in count of candidates in same beam and same cluster.
    Note that unclustered candidates are ignored.
    """

#    clsnr = []
#    cl = datal[:, 4].astype(int)   # hack. should really use table.
#    cnt_beam = datal[:, 5].astype(int)
#    cnt_cl = datal[:, 6].astype(int)
    cl = tab['cl'].astype(int)
#    cnt_beam = tab['cntb'].astype(int)
#    cnt_cl = tab['cntc'].astype(int)
    snrs = tab['snr']
    ipeak = []
    for i in np.unique(cl):
        if i == -1:
            continue
        clusterinds = np.where(i == cl)[0]
        maxsnr = snrs[clusterinds].max()
        imaxsnr = np.where(snrs == maxsnr)[0][0]
        ipeak.append(imaxsnr)
#        clsnr.append((imaxsnr, maxsnr, cnt_beam[imaxsnr], cnt_cl[imaxsnr]))
    ipeak += [i for i in range(len(tab)) if cl[i] == -1]  # append unclustered
    logger.info(f"Cluster peaks at {ipeak}:\n{tab[ipeak]}")
    print(f"Cluster peaks at {ipeak}:\n{tab[ipeak]}")

    return tab[ipeak]


def filter_clustered(tab, min_snr=None, min_dm=None, max_ibox=None, min_cntb=None, max_cntb=None, min_cntc=None,
                     max_cntc=None, target_params=None):
    """ Function to select a subset of clustered output.
    Can set minimum SNR, min/max number of beams in cluster, min/max total count in cluster.
    target_params is a tuple (min_dmt, max_dmt, min_snrt) for custom snr threshold for target.
    """

    if target_params is not None:
        min_dmt, max_dmt, min_snrt = target_params
    else:
        min_dmt, max_dmt, min_snrt = None, None, None

    good = [True] * len(tab)

    if min_snr is not None:
        if min_snrt is None:
            good *= tab['snr'] > min_snr
        else:
            #print(f'min_snr={min_snr}, min_snrt={min_snrt}, min_dmt={min_dmt}, max_dmt={max_dmt}, tab={tab[["snr", "dm"]]}')
            good0 = (tab['snr'] > min_snr)*(tab['dm'] > max_dmt)
            good1 = (tab['snr'] > min_snr)*(tab['dm'] < min_dmt)
            good2 = (tab['snr'] > min_snrt)*(tab['dm'] > min_dmt)*(tab['dm'] < max_dmt)
            good *= good0 + good1 + good2
            print('good0, good1, good2, good:')
            print(good0, good1, good2, good)
    if min_dm is not None:
        good *= tab['dm'] > min_dm
    if max_ibox is not None:
        good *= tab['ibox'] < max_ibox
    if min_cntb is not None:
        good *= tab['cntb'] > min_cntb
    if max_cntb is not None:
        good *= tab['cntb'] < max_cntb
    if min_cntc is not None:
        good *= tab['cntc'] > min_cntc
    if max_cntc is not None:
        good *= tab['cntc'] < max_cntc

    #    clsnr_out.append((imaxsnr, snr, cntb, cntc))
    tab_out = tab[good]

    logger.info(f'Filtering from {len(tab)} to {len(tab_out)} candidates.')
    print(f'Filtering from {len(tab)} to {len(tab_out)} candidates.')

    return tab_out


def dump_cluster_results_json(tab, outputfile, output_cols=['mjds', 'snr', 'ibox', 'dm', 'ibeam', 'cntb', 'cntc'], trigger=False, max_ncl=10):
    """   
    Takes tab from parse_candsfile and clsnr from get_peak, 
    output columns output_cols into a jason file outputfile. 
    trigger is bool to update DsaStore to trigger data dump.
    """

    itimes = tab['itime']
    maxsnr = tab['snr'].max()
    imaxsnr = np.where(tab['snr'] == maxsnr)[0][0]
    itime = str(itimes[imaxsnr])
    row = tab[imaxsnr]
    output_dict = {itime: {}}

#    # find time
#    try:
#        ret_time = ds.get_dict('/mon/snap/1/armed_mjd')['armed_mjd']+float(ds.get_dict('/mon/snap/1/utc_start')['utc_start'])*4.*8.192e-6/86400.
#    except:
#        ret_time = 55000.0
        
    for col in output_cols:
        if type(row[col]) == np.int64:
            output_dict[itime][col] = int(row[col])
        else:
            output_dict[itime][col] = row[col]

#    output_dict[itime]['mjds'] = ret_time
            
    if trigger and len(tab) and (len(tab) < max_ncl):
        print(f'Triggering on candidate {imaxsnr} with SNR={maxsnr}')
        logger.info(f'Triggering on candidate {imaxsnr} with SNR={maxsnr}')
        specnum = (int(itimes[imaxsnr])-offset)*downsample
        ds.put_dict('/cmd/corr/0', {'cmd': 'trigger', 'val': f'{specnum}'})
        # add output_dict to etcd
        ds.put_dict('/mon/corr/1/trigger',output_dict)
        
        with open(outputfile, 'w') as f: #encoding='utf-8'
            logger.info(f'Writing trigger info to json {outputfile}')
            json.dump(output_dict, f, ensure_ascii=False, indent=4) 
        return row
    elif len(tab) >= max_ncl:
        logger.info(f'Not triggering on block with {len(tab)} > {max_ncl} candidates')
        return None


def dump_cluster_results_heimdall(tab, outputfile): 
    """   
    Takes tab from parse_candsfile and clsnr from get_peak, 
    output T2-clustered results with the same columns as heimdall.cand into a file outputfile.
    The output is in pandas format with column names in the 1st row.
    """

# comment out to simplify parsing to two types
# 1) heimdall output
# 2) post-clustered output (adds cl, cntc, cntb)
#    output = tab['snr', 'if', 'itime', 'mjds', 'ibox', 'idm', 'dm']
#    output['members'] = tab['cntc']
#    output['ibeam'] = tab['ibeam']

    tab['itime'] = (tab['itime']-offset)*downsample  # transform to specnum

    tab.write(outputfile, format='ascii.no_header')
