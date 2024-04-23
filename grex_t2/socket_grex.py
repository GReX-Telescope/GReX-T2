import os
import numpy as np
import time
from astropy.time import Time
from astropy.io import ascii
from grex_t2 import cluster_heimdall, names
from collections import deque
import logging as logger
import requests

logger.basicConfig(filename="output.log", encoding="utf-8", level=logger.DEBUG)

nbeams_queue = deque(maxlen=10)

def filter_candidates(candsfile, output=True, trigger=True, last_trigger_time=0.0, min_dm=50, max_ibox=64, min_snr=12.5, min_snr_t2out=15.0):
    """Take a single gulp of candidates,
    parse, cluster, and then filter to
    produce highest S/N candidate and save
    to a json file
    """
    outroot = "/hdd/data/candidates/T2/"

    col_heimdall = ["snr", "if", "itime", "mjds", "ibox", "idm", "dm", "ibeam"]
    max_ncl = np.inf
    max_cntb = np.inf
    max_cntb0 = np.inf
    target_params = (50.0, 100.0, 20.0)  # Galactic bursts

    tab = ascii.read(
        candsfile, names=col_heimdall, guess=True, fast_reader=False, format="no_header"
    )

    # Ensure that the candidate table is not empty
    if not len(tab):
        return

    cluster_heimdall.cluster_data(
        tab, metric="euclidean", allow_single_cluster=True, return_clusterer=False
    )

    tab2 = cluster_heimdall.get_peak(tab)
    col_trigger = np.zeros(len(tab2), dtype=int)

    # Ensure that the candidate table is not empty
    if not len(tab2):
        return

    tab3 = cluster_heimdall.filter_clustered(
        tab2,
        min_snr=min_snr,
        min_dm=min_dm,
        max_ibox=max_ibox,
        max_cntb=max_cntb,
        max_ncl=max_ncl,
        target_params=target_params,
    )

    # Ensure that the candidate table is not empty
    if not len(tab3):
        return

    # itimes = tab3["itime"]
    # maxsnr = tab3["snr"].max()
    # imaxsnr = np.where(tab3["snr"] == maxsnr)[0][0]
    # itime_imax = str(itimes[imaxsnr])
    #   mjd = tab3["mjds"][imaxsnr]
    lastname = names.get_lastname_grex(outroot)
    cat = None
    coords = None
    snrs = None
    # prev_trig_time = None
    # min_timedelt = 60.0

    # Query for the start-time in MJD
    start_time = requests.get("http://localhost:8083/start_time").json()
    tab3["mjds"] = tab3["mjds"] / 86400. + start_time
    tab2["mjds"] = tab2["mjds"] / 86400. + start_time    

    tab4, lastname, last_trigger_time = cluster_heimdall.dump_cluster_results_json(
        tab3,
        trigger=trigger,
        lastname=lastname,
        cat=cat,
        coords=coords,
        snrs=snrs,
        outroot=outroot,
        last_trigger_time=last_trigger_time)

    if tab4 is not None and trigger:
        col_trigger = np.where(tab4 == tab2, lastname, 0)  # if trigger, then overload

    # write T2 clustered/filtered results
    if outroot is not None and len(tab2):
        tab2["trigger"] = col_trigger
        output_file = (
            outroot
            + "cluster_output"
            + str(np.floor(time.time()).astype("int"))
            + ".cand"
        )
        outputted = cluster_heimdall.dump_cluster_results_heimdall(
            tab2, output_file, min_snr_t2out=min_snr_t2out, max_ncl=max_ncl
        )

        # aggregate files
        if outputted:
            a = Time.now().mjd
            output_mjd = str(int(a))
            old_mjd = str(int(a) - 1)

            os.system("cat " + output_file + " >> " + outroot + output_mjd + ".csv")
            os.system(
                "if ! grep -Fxq 'snr,if,specnum,mjds,ibox,idm,dm,ibeam,cl,cntc,cntb,trigger' "
                + outroot
                + output_mjd
                + ".csv; then sed -i '1s/^/snr,if,specnum,mjds,ibox,idm,dm,ibeam,cl,cntc,cntb,trigger\\n/' "
                + outroot
                + output_mjd
                + ".csv; fi"
            )

            os.system(
                "echo 'snr,if,specnum,mjds,ibox,idm,dm,ibeam,cl,cntc,cntb,trigger' > "
                + outroot
                + "cluster_output.csv"
            )
            os.system(
                "test -f "
                + outroot
                + old_mjd
                + ".csv && tail -n +2 "
                + outroot
                + old_mjd
                + ".csv | tr ' ' ',' >> "
                + outroot
                + "cluster_output.csv"
            )
            os.system(
                "tail -n +2 "
                + outroot
                + output_mjd
                + ".csv | tr ' ' ',' >> "
                + outroot
                + "cluster_output.csv"
            )


# def cluster_and_plot(
#     tab,
#     globct,
#     selectcols=["itime", "idm", "ibox", "ibeam"],
#     outroot=None,
#     plot_dir=None,
#     trigger=False,
#     lastname=None,
#     max_ncl=None,
#     cat=None,
#     coords=None,
#     snrs=None,
# ):
#     """
#     Run clustering and plotting on read data.
#     Can optionally save clusters as heimdall candidate table before filtering and json version of buffer trigger.
#     lastname is name of previously triggered/named candidate.
#     cat: path to source catalog (default None)
#     beam_model: pre-calculated beam model (default None)
#     coords and snrs: from source catalog (default None)
#     """

#     # TODO: put these in json config file
#     min_dm = t2_cnf["min_dm"]  # smallest dm in filtering
#     max_ibox = t2_cnf["max_ibox"]  # largest ibox in filtering
#     min_snr = t2_cnf["min_snr"]  # smallest snr in filtering
#     min_snr_t2out = t2_cnf["min_snr_t2out"]  # smallest snr to write T2 output cand file
#     if max_ncl is None:
#         max_ncl = t2_cnf["max_ncl"]  # largest number of clusters allowed in triggering
#     max_cntb0 = t2_cnf["max_ctb0"]
#     max_cntb = t2_cnf["max_ctb"]
#     target_params = (50.0, 100.0, 20.0)  # Galactic bursts

#     # cluster
#     cluster_heimdall.cluster_data(
#         tab,
#         metric="euclidean",
#         allow_single_cluster=True,
#         return_clusterer=False,
#     )
#     tab2 = cluster_heimdall.get_peak(tab)
#     nbeams_gulp = cluster_heimdall.get_nbeams(tab2)
#     nbeams_queue.append(nbeams_gulp)

#     # Liam edit to preserve real FRBs during RFI storm:
#     # if nbeam > 100 and frac_wide < 0.8: do not discard
#     maxsnr = tab["snr"].max()
#     imaxsnr = np.where(tab["snr"] == maxsnr)[0][0]
#     cl_max = tab["cl"][imaxsnr]
#     frac_wide = np.sum(tab["ibox"][tab["cl"] == cl_max] >= 32) / float(
#         len(tab["ibox"][tab["cl"] == cl_max])
#     )

#     if len(tab["ibox"][tab["cl"] == cl_max]) == 1:
#         frac_wide = 0.0

#     # Width filter for false positives
#     ibox64_filter = False
#     if len(tab2):
#         ibox64_cnt = np.sum(tab2["ibox"] == 64) / float(len(tab2["ibox"]))
#         if ibox64_cnt > 0.85 and len(tab2["ibox"]) > 15:
#             ibox64_filter = True
#             print("ibox64 filter")

#     # Done

#     tab3 = cluster_heimdall.filter_clustered(
#         tab2,
#         min_snr=min_snr,
#         min_dm=min_dm,
#         max_ibox=max_ibox,
#         max_cntb=max_cntb,
#         max_cntb0=max_cntb0,
#         max_ncl=max_ncl,
#         target_params=target_params,
#     )  # max_ncl rows returned

#     col_trigger = np.zeros(len(tab2), dtype=int)
#     if outroot is not None and len(tab3) and not ibox64_filter:
#         tab4, lastname = cluster_heimdall.dump_cluster_results_json(
#             tab3,
#             trigger=trigger,
#             lastname=lastname,
#             cat=cat,
#             coords=coords,
#             snrs=snrs,
#             outroot=outroot,
#             nbeams=sum(nbeams_queue),
#             frac_wide=frac_wide,
#         )
#         if tab4 is not None and trigger:
#             col_trigger = np.where(
#                 tab4 == tab2, lastname, 0
#             )  # if trigger, then overload

#     # write T2 clustered/filtered results
#     if outroot is not None and len(tab2):
#         tab2["trigger"] = col_trigger
#         cluster_heimdall.dump_cluster_results_heimdall(
#             tab2,
#             outroot + str(np.floor(time.time()).astype("int")) + ".cand",
#             min_snr_t2out=min_snr_t2out,
#             max_ncl=max_ncl,
#         )

#     return lastname


def recvall(sock, n):
    """
    helper function to receive all bytes from a socket
    sock: open socket
    n: maximum number of bytes to expect. you can make this ginormous!
    """

    data = bytearray()
    while len(data) < n:
        packet = sock.recv(n - len(data))
        if not packet:
            return data
        data.extend(packet)

    return data
