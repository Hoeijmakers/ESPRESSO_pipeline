# ESPRESSO pipeline
A script that executes the ESPRESSO-pipeline using Esorex on time-series observations.




This is a single python script that passes a raw time-series dataset obtained with the [ESPRESSO spectrograph](https://www.eso.org/sci/facilities/paranal/instruments/espresso.html) through the various Data Reduction Software (DRS) recipes. The input data can be a single science exposure, or a sequence of exposures obtained using the same instrument setup (i.e. using the same calibration files), which is customarily the case for exoplanet transit observations.

It uses [Esorex](https://www.eso.org/sci/software/cpl/esorex.html) to execute the cascade of DRS recipes written for ESPRESSO, [which can be downloaded and installed here](https://www.eso.org/sci/software/pipelines/index.html#source_kit). This is perhaps not as user friendly as ESO's GUI tool named [Reflex](https://www.eso.org/sci/software/esoreflex/), but is meant as a solution for users who are unwilling or unable to use it. The rest of this ReadMe assumes that you have the pipeline and Esorex installed, and can execute pipeline recipes in the command line. The ReadMe assumes that you have not downloaded ESO data very often, so it proceeds to explain how to download a dataset in the correct format and then how to run the script. Command-line commands provided below assume that you are on a UNIX environment. On Windows, things will probably be similar. *Warning:* Running the ESPRESSO pipeline is resource-intensive. The pipeline assumes that your system has access to 16GB of RAM. My laptop doesn't, but there is a work-around using swap memory (see point 9 below).

## To get started
1. Go to [ESO's raw data archive](http://archive.eso.org/eso/eso_archive_main.html) and search for the target you wish to observe (select ESPRESSO as the instrument and increase the maximum number of returned rows to a large number to make sure you retrieve all science exposures).
2. If the target has been observed multiple times, there may be many e.g. transit datasets in the archive. Select the observations that were taken on a single date, while making sure that they indeed use the same instrument mode.
3. Request these while ticking the '+ raw calibrations' button, to use ESO's Calselector tool to automatically fish out the right calibration files that the pipeline will need (flats, darks, etc.). Download the bash-script provided by ESO that will allow you to download the data from the command line.
4. For some reason, the Calselector does not include CONTAM_FP calibration exposures, which we need to download manually. Go back to the [raw data archive](http://archive.eso.org/eso/eso_archive_main.html), and use the Start & End dates, query all ESPRESSO observations taken within a day of the observations you just requested (i.e. without specifying a target). With CTRL-F, search for CONTAM_FP. There should be multiple such exposures on the day of your observations, each corresponding to one of the two detector binning and resolution-modes (1x1, 2x1, UHR, etc.). Request these without the calselector option, giving you a second bash script.
5. Execute both bash scripts in the command line (i.e. code::`bash downloadRequest1234568.sh`). It will require authorisation using the ESO user account with which you requested the data. This creates a folder named `data_with_raw_calibs`, as well as the two CONTAM_FP files that you requested separately. Move these into the `data_with_raw_calibs` folder. The script will figure out which of the two binning factors to use.
6. The fits files are compressed (they have the extension .fits.Z). Uncompress them by running code::`uncompress *.Z`
or code::`gunzip *.Z` in the command line.
7. Determine whether the observations were taken with the secondary fiber B on sky or with the Fabry-Perot, and in what
binning mode the data were taken. To check which read-out mode your exposures were obtained in, open a random science frame and search for the `HIERARCH ESO DET BINX` and `HIERARCH ESO DET BINY` keywords in the fits header. These tell you the binning factor.
8. Run the script as code::`python3 espresso_pipeline.py inpath outpath binning fiber_B`. In this example, the inpath variable would point to the `data_with_raw_calibs` folder; outpath to a *local* folder on your machine (remote folders or folders on external drives are sometimes problematic because they may use different file systems, giving you an OS-error in python). The `binning` variable should be set to the correct binning factor of your observations, e.g. `1x1` or `2x1`.
Set the fiber_B variable to FP for observations taken with the Fabry-Perot, or sky if fiber B is on sky.
By default, observations are assumed to be in `2x1` mode and with Fiber B on sky. If that is correct, the last two
input variables can be omitted, e.g. as code::`python3 espresso_pipeline.py inpath outpath`.
9. Reducing a single dataset on my laptop takes hours, many GBs of disk space and close to 16 GB of RAM. The latter could be a problem (my laptop has 4GB of RAM only): If your computer does not have sufficient RAM available, the code will crash halfway through. To alleviate this, you can assign swap memory to increase your RAM capacity, as follows.  (adopted from <https://linuxize.com/post/create-a-linux-swap-file/>). Although this is much slower than using RAM, at least it will allow you to run the pipeline even if you don't have enough RAM.<br>
   `sudo dd if=/dev/zero of=/swapfile bs=1024 count=10000000`<br>
   `sudo chmod 600 /swapfile`<br>
   `sudo mkswap /swapfile`<br>
   `sudo swapon /swapfile`<br>
   To see that it has worked, hit `sudo swapon -show`.

<br>
<br>
After a few hours, this should have provided you with pipeline-reduced ESPRESSO data!
