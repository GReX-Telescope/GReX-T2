import json
import os.path
import socket
import hdbscan
import numpy as np
from astropy.io import ascii
from astropy.io.ascii.core import InconsistentTableError
from astropy.time import Time
from numpy.lib.recfunctions import structured_to_unstructured
from grex_t2 import triggering, names
import logging as logger

logger.basicConfig(filename="logs/output.log", encoding="utf-8", level=logger.DEBUG)

# half second at heimdall time resolution (after march 18)
OFFSET = 1907
DOWNSAMPLE = 4


def parse_candsfile(candsfile):
    """Takes standard MBHeimdall giants output and returns full table,
    classifier inputs and snr tables.
    (Can add cleaning here, eventually)
    """

    if os.path.exists(candsfile):
        logger.debug(f"Candsfile {candsfile} is path, so opening it")
        candsfile = open(candsfile, "r").read()
    else:
        ncands = len(candsfile.split("\n")) - 1
        logger.debug(f"Received {ncands} candidates")
    col_heimdall = ["snr", "if", "itime", "mjds", "ibox", "idm", "dm", "ibeam"]
    col_T2old = [
        "snr",
        "if",
        "itime",
        "mjds",
        "ibox",
        "idm",
        "dm",
        "ibeam",
        "cl",
        "cntc",
        "cntb",
    ]
    col_T2 = [
        "snr",
        "if",
        "itime",
        "mjds",
        "ibox",
        "idm",
        "dm",
        "ibeam",
        "cl",
        "cntc",
        "cntb",
        "trigger",
    ]

    # flag for heimdall file
    hdfile = False

    try:
        tab = ascii.read(
            candsfile,
            names=col_heimdall,
            guess=True,
            fast_reader=False,
            format="no_header",
        )
        hdfile = True
        logger.debug("Read with heimdall columns")
    except InconsistentTableError:
        try:
            tab = ascii.read(
                candsfile,
                names=col_T2,
                guess=True,
                fast_reader=False,
                format="no_header",
            )
            hdfile = False
            logger.debug("Read with T2 columns")
        except InconsistentTableError:
            try:
                tab = ascii.read(
                    candsfile,
                    names=col_T2old,
                    guess=True,
                    fast_reader=False,
                    format="no_header",
                )
                hdfile = False
                logger.debug("Read with old style T2 columns")
            except InconsistentTableError:
                logger.warning("Inconsistent table. Skipping...")
                return ([], [], [])

    tab["ibeam"] = tab["ibeam"].astype(int)

    start_time_mjd = requests.get("http://localhost:8083/start_time").json()
    tab["mjds"] = tab["mjds"] / 86400.0 + start_time_mjd

    return tab


def dm_range(dm_max, dm_min=5.0, frac=0.2):
    """Generate list of DM-windows in which
    to search for single pulse groups.

    Parameters
    ----------
    dm_max : float
        max DM
    dm_min : float
        min DM
    frac : float
        fractional size of each window

    Returns
    -------
    dm_list : list
        list of tuples containing (min, max) of each
        DM window
    """

    dm_list = []
    prefac = (1 - frac) / (1 + frac)

    while dm_max > dm_min:
        if dm_max < 100.0:
            prefac = (1 - 2 * frac) / (1 + 2 * frac)
        if dm_max < 50.0:
            prefac = 0.0

        dm_list.append((int(prefac * dm_max), int(dm_max)))
        dm_max = int(prefac * dm_max)

    return dm_list


def cluster_data(
    tab,
    selectcols=["itime", "idm", "ibox", "ibeam"],
    min_cluster_size=2,
    min_samples=5,
    metric="hamming",
    return_clusterer=False,
    allow_single_cluster=True,
):
    """Take data from parse_candsfile and identify clusters
    via hamming metric.
    selectcols will take a subset of the standard MBHeimdall output
    """

    data = structured_to_unstructured(
        tab[selectcols].as_array()
    )  # ok for single dtype (int)
    clusterer = None
    try:
        clusterer = hdbscan.HDBSCAN(
            metric=metric,
            min_cluster_size=min_cluster_size,
            min_samples=min_samples,
            cluster_selection_method="eom",
            allow_single_cluster=allow_single_cluster,
        ).fit(data)

        cl = clusterer.labels_
    except ValueError:
        print(
            "Clustering did not run. Each point \
               assigned to unique cluster."
        )
        logger.info(
            "Clustering did not run. Each point \
               assigned to unique cluster."
        )
        cl = np.arange(len(data))

    # hack assumes fixed columns
    bl = data[:, 3]
    cntb, cntc = np.zeros((len(data), 1), dtype=int), np.zeros(
        (len(data), 1), dtype=int
    )
    ucl = np.unique(cl)

    for i in ucl:
        ww = np.where(i == cl)
        cntc[ww] = len(ww[0])
        ubl = np.unique(bl[ww])
        cntb[ww] = len(ubl)

    # modifies tab in place
    tab["cl"] = cl.tolist()
    tab["cntc"] = cntc.flatten().tolist()
    tab["cntb"] = cntb.flatten().tolist()

    if return_clusterer:
        return clusterer


def get_peak(tab):
    """Given labeled data, find max snr row per cluster
    Adds in count of candidates in same beam and same cluster.
    Puts unclustered candidates in as individual events.
    """

    cl = tab["cl"].astype(int)
    snrs = tab["snr"]
    ipeak = []
    for i in np.unique(cl):
        if i == -1:
            continue
        clusterinds = np.where(i == cl)[0]
        maxsnr = snrs[clusterinds].max()
        imaxsnr = np.where(snrs == maxsnr)[0][0]
        ipeak.append(imaxsnr)
    # Append unclustered
    ipeak += [i for i in range(len(tab)) if cl[i] == -1]
    logger.info(f"Found {len(ipeak)} cluster peaks")
    print(f"Found {len(ipeak)} cluster peaks")

    return tab[ipeak]


def filter_clustered(
    tab,
    min_snr=None,
    min_dm=None,
    max_ibox=None,
    min_cntb=None,
    max_cntb=None,
    min_cntc=None,
    max_cntc=None,
    max_ncl=None,
    target_params=None,
):
    """Function to select a subset of clustered output.
    Can set minimum SNR, min/max number of beams in cluster,
    min/max total count in cluster.
    target_params is a tuple (min_dmt, max_dmt, min_snrt)
    for custom snr threshold for target.
    max_ncl is maximum number of clusters returned (sorted by SNR).
    """

    if target_params is not None:
        min_dmt, max_dmt, min_snrt = target_params
    else:
        min_dmt, max_dmt, min_snrt = None, None, None

    good = [True] * len(tab)

    if min_snr is not None:
        if min_snrt is None:
            good *= tab["snr"] > min_snr
        else:
            good0 = (tab["snr"] > min_snr) * (tab["dm"] > max_dmt)
            good1 = (tab["snr"] > min_snr) * (tab["dm"] < min_dmt)
            good2 = (
                (tab["snr"] > min_snrt) * (tab["dm"] > min_dmt) * (tab["dm"] < max_dmt)
            )
            good *= good0 + good1 + good2

    if min_dm is not None:
        good *= tab["dm"] > min_dm
    if max_ibox is not None:
        good *= tab["ibox"] < max_ibox
    if min_cntb is not None:
        good *= tab["cntb"] > min_cntb
    if max_cntb is not None:
        good *= tab["cntb"] < max_cntb
    if min_cntc is not None:
        good *= tab["cntc"] > min_cntc
    if max_cntc is not None:
        good *= tab["cntc"] < max_cntc

    tab_out = tab[good]

    if max_ncl is not None:
        if len(tab_out) > max_ncl:
            min_snr_cl = sorted(tab_out["snr"])[-max_ncl]
            good = tab_out["snr"] >= min_snr_cl
            tab_out = tab_out[good]
            print(
                f"Limiting output to {max_ncl} \
                    clusters with snr>{min_snr_cl}."
            )

    logger.info(f"Filtering clusters from {len(tab)} to {len(tab_out)} candidates.")
    print(f"Filtering clusters from {len(tab)} to {len(tab_out)} candidates.")

    return tab_out


def dump_cluster_results_json(
    tab,
    outputfile=None,
    output_cols=["mjds", "snr", "ibox", "dm", "ibeam", "cntb", "cntc"],
    trigger=False,
    lastname=None,
    cat=None,
    coords=None,
    snrs=None,
    outroot="./",
    injectionfile=None,
    last_trigger_time=0.0,
):
    """
    Takes tab from parse_candsfile and clsnr from get_peak,
    json file will be named with generated name, unless outputfile is set
    candidate name and specnum is calculated. name is unique.
    trigger is bool to update DsaStore to trigger data dump.
    cat is path to source catalog (default None)
    beam_model is pre-calculated beam model (default None)
    coords and snrs are parsed source file input
    injectionfile is path to info on injects and controls whether trigger is compared to that
    returns row of table that triggered, along with name generated for candidate.
    """

    if coords is None or snrs is None:
        coords, snrs = triggering.parse_catalog(cat)

    itimes = tab["itime"]
    maxsnr = tab["snr"].max()
    imaxsnr = np.where(tab["snr"] == maxsnr)[0][0]
    # itime = str(itimes[imaxsnr])
    specnum = (int(itimes[imaxsnr]) - OFFSET) * DOWNSAMPLE
    mjd = tab["mjds"][imaxsnr]
    # _snr = tab["snr"][imaxsnr]
    dm = tab["dm"][imaxsnr]
    ibeam = tab["ibeam"][imaxsnr]

    # if no injection file or no coincident injection
    candname = names.increment_name(mjd, lastname=lastname)

    isinjection = False
    if injectionfile is not None:
        # check candidate against injectionfile
        tab_inj = ascii.read(injectionfile)
        assert all(
            [col in tab_inj.columns for col in ["MJD", "Beam", "DM", "SNR", "FRBno"]]
        )

        # is candidate proximal to any in tab_inj?
        t_close = 15  # seconds  TODO: why not 1 sec?
        dm_close = 10  # pc/cm3
        beam_close = 2  # number
        sel_t = np.abs(tab_inj["MJD"] - mjd) < t_close / (3600 * 24)
        sel_dm = np.abs(tab_inj["DM"] - dm) < dm_close
        sel_beam = np.abs(tab_inj["Beam"] - ibeam) < beam_close
        sel = sel_t * sel_dm * sel_beam
        if len(np.where(sel)[0]):
            isinjection = True

        if isinjection:
            basename = names.increment_name(mjd, lastname=lastname)
            candname = f"{basename}_inj{tab_inj[sel]['FRBno'][0]}"
            print(f"Candidate identified as injection. Naming it {candname}")
            if len(sel) > 1:
                print(
                    f"Found {len(sel)} injections coincident with this event. Using first."
                )

    output_dict = {candname: {}}
    if outputfile is None:
        outputfile = f"{outroot}{candname}.json"

    row = tab[imaxsnr]
    red_tab = tab[imaxsnr : imaxsnr + 1]
    for col in output_cols:
        if type(row[col]) == np.int64:
            output_dict[candname][col] = int(row[col])
        else:
            output_dict[candname][col] = row[col]

    output_dict[candname]["specnum"] = specnum

    if len(tab) > 0:
        print("\n", red_tab, "\n")
        if cat is not None and red_tab is not None:
            tab_checked = triggering.check_clustered_sources(
                red_tab, coords, snrs, do_check=False
            )
            if len(tab_checked):
                with open(outputfile, "w") as f:  # encoding='utf-8'
                    print(f"Writing trigger file for index {imaxsnr} with SNR={maxsnr}")
                    logger.info(
                        f"Writing trigger file for index {imaxsnr} with SNR={maxsnr}"
                    )
                    json.dump(output_dict, f, ensure_ascii=False, indent=4)

                trigger = True

                if trigger:
                    if (Time.now().mjd - last_trigger_time) < 90.0/86400.:
                        print(
                            "Not triggering on source in beam as last trigger was < 90s ago"
                        )
                        logger.info(
                            "Not triggering on source in beam as last trigger was < 90s ago"
                        )
                        return None, lastname, last_trigger_time
                    
                    print(output_dict)
                    send_trigger(output_dict=output_dict)
                    last_trigger_time = Time.now().mjd

                return row, candname, last_trigger_time

            else:
                print("Not triggering on source in beam")
                logger.info("Not triggering on source in beam")
                return None, lastname, last_trigger_time

        else:
            with open(outputfile, "w") as f:  # encoding='utf-8'
                print(f"Writing trigger file for index {imaxsnr} with SNR={maxsnr}")
                logger.info(
                    f"Writing trigger file for index {imaxsnr} with SNR={maxsnr}"
                )
                json.dump(output_dict, f, ensure_ascii=False, indent=4)

            if trigger:  # and not isinjection ?
                if (Time.now().mjd - last_trigger_time) < 90.0/86400.:
                    print(
                        "Not triggering on source in beam as last trigger was < 90s ago"
                    )
                    logger.info(
                        "Not triggering on source in beam as last trigger was < 90s ago"
                    )
                    return None, lastname, last_trigger_time
                print(output_dict)
                send_trigger(output_dict=output_dict)
                last_trigger_time = Time.now().mjd

            return row, candname, last_trigger_time

    else:
        print(f"Not triggering on block with {len(tab)} candidates")
        logger.info(f"Not triggering on block with {len(tab)} candidates")
        return None, lastname, last_trigger_time


def send_trigger(output_dict=None, outputfile=None):
    """Use either json file or dict to send trigger for voltage dumps via udp."""

    if outputfile is not None:
        print("Overloading output_dict trigger info with that from outputfile")
        logger.info("Overloading output_dict trigger info with that from outputfile")
        with open(outputfile, "w") as f:
            output_dict = json.load(f)

    if output_dict is not None:
        candname = list(output_dict)[0]
        val = output_dict.get(candname)
        nowmjd = Time.now().mjd
        print(f"Sending trigger for candidate {candname} with specnum {val['specnum']} at Time", nowmjd)
        logger.info(
            f"Sending trigger for candidate {candname} with specnum {val['specnum']}"
        )
        trigger_info = {"key": candname, "spectrum": val['spectrum']}
        trigger_message = json.dumps(trigger_info).encode('utf-8')
        
        UDP_PORT = 65432
        UDP_IP = "127.0.0.1"
        sock = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)  # Internet  # UDP
        sock.sendto(trigger_message, (UDP_IP, UDP_PORT))


def dump_cluster_results_heimdall(tab, outputfile, 
                                  min_snr_t2out=None, 
                                  max_ncl=None):
    """
    Takes tab from parse_candsfile and clsnr from get_peak,
    output T2-clustered results with the same columns as
    heimdall.cand into a file outputfile.
    The output is in pandas format with column names
    in the 1st row.
    min_snr_t2out is a min snr on candidates to write.
    max_ncl is number of rows to write.
    """

    # transform to specnum
    tab["itime"] = (tab["itime"] - OFFSET) * DOWNSAMPLE

    if min_snr_t2out is not None:
        good = [True] * len(tab)
        good *= tab["snr"] > min_snr_t2out
        tab = tab[good]
        if not all(good) and len(tab):
            print(f"Limiting output to SNR>{min_snr_t2out} with {len(tab)} clusters.")

    if max_ncl is not None:
        if len(tab) > max_ncl:
            min_snr_cl = sorted(tab["snr"])[-max_ncl]
            good = (tab["snr"] >= min_snr_cl) + [
                str(tt) != "0" for tt in tab["trigger"]
            ]  # keep trigger
            tab = tab[good]
            print(f"Limiting output to {max_ncl} clusters with snr>{min_snr_cl}.")
    else:
        print("max_ncl not set. Not filtering heimdall output file.")

    if len(tab) > 0:
        tab.write(outputfile, format="ascii.no_header", overwrite=True)
        return True

    return False
