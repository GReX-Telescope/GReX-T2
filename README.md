[![build status](https://img.shields.io/github/workflow/status/GReX-Telescope/GReX-T2/CI/main?style=flat-square&logo=github)](https://github.com/GReX-Telescope/GReX-T2/actions)
[![codecov](https://codecov.io/gh/GReX-Telescope/GReX-T2/branch/master/graph/badge.svg)](https://codecov.io/gh/GReX-Telescope/GReX-T2)

# GReX-T2

Module for parsing Heimdall output and clustering candidates to trigger voltage recording.

## Install

`python setup.py install`

## Dependencies

- python 3.6+
- astropy
- numpy
- hdbscan

## Running tests

1. install development version
2. `pytest` (all should pass)
3. `python scripts/run_T2.py` (waits for socket connection on port 12345)
4. run heimdall
   `./heimdall -f ~gcchen/fake_data/test_32beams_10p_5w.fil -output_dir ~/ -dm 10 2000 -dm_tol 1.25 -nsamps_gulp 4096 -detect_thresh 8 -boxcar_max 12 -nbeams 32 -v -coincidencer 127.0.0.1:12345`

## Test

`pytest`
