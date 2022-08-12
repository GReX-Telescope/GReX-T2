import sys

import matplotlib.pyplot as plt
import numpy as np
import time

from T2 import cluster_heimdall, plotting

# my_cnf = cnf.Conf(use_etcd=True)
# t2_cnf = my_cnf.get('t2')
# min_dm = t2_cnf['min_dm']  # smallest dm in filtering
# max_ibox = t2_cnf['max_ibox']  # largest ibox in filtering
# min_snr = t2_cnf['min_snr']  # smallest snr in filtering
# min_snr_t2out = t2_cnf['min_snr_t2out']  # smallest snr to write T2 output cand file
# #if max_ncl is None:
# max_ncl = t2_cnf['max_ncl']  # largest number of clusters allowed in triggering
# max_cntb0 = t2_cnf['max_ctb0']
# max_cntb = t2_cnf['max_ctb']

candsfile = sys.argv[1]

max_cntb = np.inf
min_dm = 25.0
max_ibox = 64
nbeams_queue = []
tab = cluster_heimdall.parse_candsfile(candsfile)
tabdumb = cluster_heimdall.cluster_dumb(tab, t_window=0.5)

cluster_heimdall.cluster_data(
    tab, metric="euclidean", allow_single_cluster=True, return_clusterer=False
)

tab2 = cluster_heimdall.get_peak(tab)
nbeams_gulp = cluster_heimdall.get_nbeams(tab2)
nbeams_queue.append(nbeams_gulp)
tab3 = cluster_heimdall.filter_clustered(
    tab2,
    min_snr=9.3,
    min_dm=min_dm,
    max_ibox=max_ibox,
    max_cntb=max_cntb,
    target_params=target_params,
)


tab4, lastname = cluster_heimdall.dump_cluster_results_json(
    tab3,
    trigger=False,
    lastname="171010aaaa",
    cat=None,
    beam_model=None,
    coords=None,
    snrs=None,
    outroot=None,
    nbeams=sum(nbeams_queue[-10:]),
)
max_point_size = 5.0
plt.scatter(
    (tab["mjds"] - tab["mjds"].min()) * 86400.0,
    tab["dm"],
    s=(max_point_size * tab["snr"] / tab["snr"].mean()),
    alpha=0.5,
)
plt.scatter(
    (tab2["mjds"] - tab2["mjds"].min()) * 86400.0,
    tab2["dm"],
    s=3.5,
)
plt.scatter(
    (tabdumb["mjds"] - tabdumb["mjds"].min()) * 86400.0,
    tabdumb["dm"],
    s=3.5,
)
plt.scatter((tab3["mjds"] - tab3["mjds"].min()) * 86400.0, tab3["dm"], s=3.5, c="k")

plt.xlabel("Time (s)")
plt.ylabel("DM (pc cm**-3)")
plt.legend(["Heimdall candidates", "Clustered (HDBSCAN)", "Clustered (Dumb)"])
plt.savefig("here.png")
