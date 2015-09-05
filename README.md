[![Build Status](https://travis-ci.org/guruucsd/PING.svg?branch=master)](https://travis-ci.org/guruucsd/PING)
[![Coverage Status](https://coveralls.io/repos/guruucsd/PING/badge.svg?branch=master&service=github)](https://coveralls.io/github/guruucsd/PING?branch=master)

PING automation scripts
=======================

This is a set of scripts to use with the Pediatric Imaging, Neurocognition, and Genetics (PING) study. Please see http://pingstudy.ucsd.edu/ for more information about this study.

# Installation

All installation dependencies are defined in the `requirements.txt` file.

Installation steps:

1. `git clone git@github.com:guruucsd/PING.git`
2. `cd PING`
3. `sudo pip install -r requirements.txt`  If this fails:
    * Make sure you have pip installed: https://pip.pypa.io/en/latest/installing.html
    * If it fails installing numpy or scipy, I suggest installing [miniconda](http://conda.pydata.org/miniconda.html) and installing them via `conda install numpy scipy pip`
4. Visit https://ping-dataportal.ucsd.edu/applications/User/requestLogin.php and request access to the PING data portal.


# Scripts

The following functions of the data portal are available through this scripting interface:

* `export.py` - Compute some derived measures from the PING data, and export to a local CSV
* `grouping.py` - Plot regressions, separating data by the chosen groupings.
* `gwas.py` - **PRELIMINARY--you must request access before using!** download GWAS results, or launch a new GWAS study.
* `scatter.py` - Show scatter plots, filtering columns by prefix and choosing computations for x, y, and point size. 
* `similarity.py` - Script interface for viewing similarity matrices between different PING measures. Select measures by prefix, and use correlation or partial correlation.
* `snps.py` - Script interface for the SNP browser; find genes associated with SNPs, SNPs associated with genes, and download user data for specific genes. **NOTE: PING limits SNP downloads to 5000, so use them wisely!**
* `upload.py` - Upload your local CSV with derived measures to the PING data portal, so you can access in the Data browser.

### Authenticating in scripts

Scripts must be run with a valid username and password to the PING portal--if you don't have one, you'll need to [request one](https://ping-dataportal.ucsd.edu/applications/User/requestLogin.php).

There are two ways to pass your username and password to the scripts:

* Via the command-line: use parameters `--username=[USERNAME] --password=[PASSWORD]`
* Via the shell: set the following environment variables:
    * `PING_USERNAME` - your PING data portal username (*not* your email address).
    * `PING_PASSWORD` - your PING data portal password


### Example script usage

These examples assume that you've set your username and password via shell environment variables. If you haven't, you'll need to add the `--username=[USERNAME]` and `--password=[PASSWORD]` options to the commands below.

* `python export.py` - exports data sheet, including computed measures, to a local CSV file. 
* `python grouping.py MRI_cort_area.ctx Gender` - show linear regression for each cortical area measure, with regression by gender overlaid on the same plot.
* `python scatter.py MRI_cort_area AI:mean AI:std LH_PLUS_RH:mean` - show a scatter plot over cortical area measures, of asymmetry index mean vs. standard deviation, with dot size given by the total area (LH+RH) 
* `python similarity.py MRI_cort_area partial-correlation "Left Hemisphere"` - show a similarity matrix between all measures with the MRI_cort_area prefix, using partial correlation, for the left hemisphere only.
* `python snps.py view STK31` - show all SNPs associated with gene STK31, according to the PING genetics DB.


