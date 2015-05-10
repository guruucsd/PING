PING automation scripts
=======================

This is a set of scripts to use with the Pediatric Imaging, Neurocognition, and Genetics (PING) study. Please see http://pingstudy.ucsd.edu/ for more information about this study.

# Installation

All installation dependencies are defined in the `requirements.txt` file.

Installation steps:

1. `git clone git@github.com/bcipolli/PING.git`
2. `cd PING`
3. `sudo pip install -r requirements.txt`
4. Visit https://ping-dataportal.ucsd.edu/applications/User/requestLogin.php and request access to the PING data portal.


# Usage

Scripts can be modified to accept your username and password, or can be run as-is if you set the following environment variables:

* `PING_USERNAME` - your PING data portal username (*not* your email address).
* `PING_PASSWORD` - your PING data portal password

The following functions of the data portal are available through this scripting interface:

* `gwas.py` - **PRELIMINARY--you must request access before using!** download GWAS results, or launch a new GWAS study.
* `export.py` - Compute some derived measures from the PING data, and export to a local CSV
* `upload.py` - Upload your local CSV with derived measures to the PING data portal, so you can access in the Data browser.
* `regress.py` - Script interface for the Data Visualization tool (computing regressions, downloading CSV results)
* `snps.py` - Script interface for the SNP browser; find genes associated with SNPs, SNPs associated with genes, and download user data for specific genes. **NOTE: PING limits SNP downloads to 5000, so use them wisely!**

