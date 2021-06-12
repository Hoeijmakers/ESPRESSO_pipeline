#This script is a wrapper for the execution of the ESPRESSO pipeline using ESOREX.
#ESOREX is a command-line utility that allows the user to call data reduction recipes
#one by one, allowing the step-wise execution of the entire reduction cascade.
#That would be very neat if the process of getting all the files in the right order wasn't
#so convoluted. The purpose of this script is to make the latter as easy as possible.

#IMPORTANT: the execution of this code generates a lot of temporary files that are later
#removed (by-product of using ESOREX). You are adviced to run this python script in an
#empty work folder (not the output folder).

#First, you need to make sure that both ESOREX and the ESPRESSO pipeline are installed on
#your system. To check whether ESOREX is installed correctly, you can open a terminal and type:

#>>>  esorex

#If this command is recognised, ESOREX has been installed.
#This script was tested to work on ESOREX version 3.13.1.

#Next, test wether the ESPRESSO data reduction recipes have been installed. In the terminal, type:

#>>>   esorex -recipes

#This will return a list of all the recipes installed on your system. It should include such recipes
#as espdr_mdark, espdr_mbias, and a bunch of others starting with espdr.
#This script was built around version 1.2.2 of the ESPRESSO pipeline and was upgraded to run with version 2.2.1 ***

#A recipe needs to be called with a single argument: The path to a text-file that lists the locations
#of the files necessary for the execution of the script. For example, the espdr_mbias recipe needs to
#be fed by a list that details the locations of the BIAS frames. These files are called "SOF files".

#The script assumes that the user has created *one* folder on their system that contains *all*
#of the relevant frames *of a single observing mode*. ESPRESSO has multiple observing modes,
#varying the spectral resolution, how many UTs were used and the detector binning.

#This script assumes that observations are taken in SINGLEHR mode. The user can choose the binning
#factor; but all the relevant data of a single mode needs to be in the same folder, *without any files
#obtained in different modes, if any*! I did not do the effort of checking all the relevant keywords
#of every single file... so make sure you do this right.

#Version 2 of the ESPRESSO pipeline no longer contains all static calibration files. Instead, these are provided
#by the ESO calselector tool when downloading the data. This means that you NEED to use the ESO calselector to
# download the data.

#These two constraints imply the following work flow:
#1) In the RAW data archive, find a science observing sequence (i.e. during a single night, obtained in the same mode/setup) that you would like to reduce.
#2) Request all these Science frames, and enable the calselector to provide them + associated raw calibrations.
#3) Download these to a single folder (for an exoplanet transit observation, this can easily amount to >10GB).
#4) Below (at the end of the script), provide the path to this folder as inpath, and the path to a folder as outpath.
#5) For some reason, the calselector does not select the CONTAM,OFF,FP calibration files. These are obtained every night in 1x1 and 2x1 binning
# modes. To download these, you need to query the entire night of observations with the Raw data selection form, and find it (ctrl-F) to
#request and download it manually. Notice that there might be multiple CONTAM,OFF,FP files, but only one can
#be used. Likely either will work (as long as they are in the right binning).
#5) Make sure that your calib folder is set correctly, and run this script.


#The script first creates the necessary SOF files with which the recipes are called, before executing
#the entire cascade one recipe after another, all the while moving the necessary intermediate products
#to the output folder.


#That's all nice and well, but be aware that the pipeline requires a large amount of both *computer
#resources and time*. Each script can take one or even several hours to run. Espdr_mflat took over
#8000 seconds to run on my laptop; which becomes nearly unusable during this time.
#So set aside a decent number of hours on your machine to do this.

#More importantly however, some recipes are very heavy on RAM and diskspace. The pipeline is designed
#to run on systems with 16 GB of RAM. My run-off-the-mill laptop has 4. So I am condemned to using a lot
#of swap memory, which is slow; causing the recipes to take longer. However, adding swap space does do the
#trick, and given a couple of days, this script will successfully reduce your data.

#Adding ~10GB additional swap memory to a file called swapfile located at / is done as follows,
#adopted from https://linuxize.com/post/create-a-linux-swap-file/:

#sudo dd if=/dev/zero of=/swapfile bs=1024 count=10000000
#sudo chmod 600 /swapfile
#sudo mkswap /swapfile
#sudo swapon /swapfile

#To verify that it has worked, check the swap space using:
#sudo swapon -show

#To unload the swap file from memory:
#sudo swapoff /swapfile
#Remove the entry 'swapfile swap swap defaults 0 0' from /etc/fstab
#sudo rm /swapfile

#The swap file is dismounted during a restart, so when restarting you need to mount again with swapon.




# ***
#As the pipeline is still in flux, certain filenames may be named differently in the future. I foresee
#this to be the main potential challenge of the usage of this script in the future. If this happens,
#the code below needs to be modified to take into account the newer filenames. That's will be a bit of
#a tedious job, but should be doable with the pipeline manual in hand. The main things that would probably be affected in
#such a scenario would be the static calibration files, which are located somewhere in your installation
#of esorex. In my case, for example, this would be in /home/jens/ESOREX/calib/espdr-2.1.1/cal/.













#==============================================================================================#
#                             Start of function definitions
#==============================================================================================#


def create_sof(inpath,outpath,path_cdb,binning,sky=True):
    """This script creates the file association lists (sof files) that are the main inputs
    to the pipeline recipes when called with esorex. The user provides the path of the raw data files
    (inpath) as downloaded from the ESO archive. These must be sorted by instrument mode
    (i.e. resolution and binning). The user provides the output path, to which the sof files will be
    written (outpath).  This is also the location in which the pipeline reduction products will be
    stored. The user provies the path to the static calibration files (path_cdb) which are
    installed with the pipeline. By default this could be something like:
    /home/jens/esorex_install_folder/calib/espdr-1.2.2/cal/. Finally, the user manually provides the
    binning factor as a string (binning), as a string. Valid inputs are '1x1' or '2x1'. Other entries
    will result in a crash when executing the recipes.

    This recipe assumes that the object is taken with fiber B on sky. If fiber B is FB, change the sky keyword to False.

    """
    import os
    import numpy as np
    import astropy.io.fits as fitsio
    import pdb
    import sys
    from pathlib import Path
    import glob
    #Start by defining the input and output paths:

    path_out=outpath

    try:
        os.mkdir(outpath)
    except:
        pass

    #The rest works automatically.
    file_list = os.listdir(inpath)
    file_list = [str(i) for i in Path(inpath).glob('ESPRE*.fits')]
    static_list = [str(i) for i in Path(inpath).glob('M.ESPRESSO*.fits')]
    mask_list =[]
    fits_list=[]
    type_list=[]
    dits_list=[]
    binx_list=[]
    biny_list=[]
    static_type_list=[]
    binstring = binning.split('x')
    binx = int(binstring[0])
    biny = int(binstring[1])


    for file in static_list:
        with fitsio.open(file) as fu:
            static_type_list.append(fu[0].header['ESO PRO CATG'])


    static_dict = dict()#We save the statics in a dictionary so that they can be parsed easily later.
    for i,s in enumerate(static_type_list):
        if s == 'MASK_TABLE':
            mask_list.append(static_list[i])
        else:
            static_dict[s]=static_list[i]

    #for i in static_dict:
    #    print(i)
    #for i in mask_list:
    #    print(i)
    #sys.exit()
    for file in file_list:
        fits_list=np.append(fits_list,file)
        with fitsio.open(file) as fu:
            type_list=np.append(type_list,fu[0].header['HIERARCH ESO DPR TYPE'])
            dits_list=np.append(dits_list,fu[0].header['EXPTIME'])
            binx_list=np.append(binx_list,fu[0].header['HIERARCH ESO DET BINX'])
            biny_list=np.append(biny_list,fu[0].header['HIERARCH ESO DET BINY'])

    for i in range(len(type_list)):
        print(type_list[i]+'  %s x %s' % (int(binx_list[i]),int(biny_list[i])))
    binning2 = 'VOID' #This will cause it to crash when trying to find calib 11 or 21 files down the road, unless binning is set properly to 1x1 or 2x1.


    #The following is to switch between different sky modes.
    if sky == True:
        object_keyword='OBJECT,SKY'
        object_tag='OBJ_SKY'
    else:
        object_keyword='OBJECT,FP'
        object_tag='OBJ_FP'


    #Define the lists in which the frame types will be sorted.
    bias_list=[]
    dark_list=[]
    LED_list=[]
    orderdefA_list=[]
    orderdefB_list=[]
    flatA_list=[]
    flatB_list=[]
    contam_list=[]
    eff_list=[]
    std_list=[]
    sci_list=[]
    FP_FP_list=[]
    FP_TH_list=[]
    TH_FP_list=[]



    #The following goes through the list of user-supplied files and sorts them along type.
    for i in range(len(fits_list)):
        if type_list[i] == 'BIAS' and binx_list[i] == binx and biny_list[i] == biny:
            bias_list = np.append(bias_list,fits_list[i]+'   '+type_list[i])
        if type_list[i] == 'DARK' and binx_list[i] == binx and biny_list[i] == biny:
            dark_list = np.append(dark_list,fits_list[i]+'   '+type_list[i])#+'   %s' % dits_list[i])
            #For the DARKS we dont care about the exptime. It is 3600 for all of them....
        if type_list[i] == 'LED' and binx_list[i] == binx and biny_list[i] == biny:
            LED_list = np.append(LED_list,fits_list[i]+'   '+type_list[i]+'_FF')
        if type_list[i] == 'ORDERDEF,LAMP,OFF' and binx_list[i] == binx and biny_list[i] == biny:
            orderdefA_list = np.append(orderdefA_list,fits_list[i]+'   '+'ORDERDEF_A')
        if type_list[i] == 'ORDERDEF,OFF,LAMP' and binx_list[i] == binx and biny_list[i] == biny:
            orderdefB_list = np.append(orderdefB_list,fits_list[i]+'   '+'ORDERDEF_B')
        if type_list[i] == 'FLAT,LAMP,OFF' and binx_list[i] == binx and biny_list[i] == biny:
            flatA_list = np.append(flatA_list,fits_list[i]+'   '+'FLAT_A')
        if type_list[i] == 'FLAT,OFF,LAMP' and binx_list[i] == binx and biny_list[i] == biny:
            flatB_list = np.append(flatB_list,fits_list[i]+'   '+'FLAT_B')
        if type_list[i] == 'WAVE,FP,FP' and binx_list[i] == binx and biny_list[i] == biny:
            FP_FP_list = np.append(FP_FP_list,fits_list[i]+'   '+'FP_FP')
        if type_list[i] == 'WAVE,FP,THAR' and binx_list[i] == binx and biny_list[i] == biny:
            FP_TH_list = np.append(FP_TH_list,fits_list[i]+'   '+'FP_THAR')
        if type_list[i] == 'WAVE,THAR,FP' and binx_list[i] == binx and biny_list[i] == biny:
            TH_FP_list = np.append(TH_FP_list,fits_list[i]+'   '+'THAR_FP')
        if type_list[i] == 'CONTAM,OFF,FP' and binx_list[i] == binx and biny_list[i] == biny:
            contam_list = np.append(contam_list,fits_list[i]+'   '+'RAW_CONTAM_FP')
        if type_list[i] == 'EFF,SKY,SKY' and binx_list[i] == binx and biny_list[i] == biny:
            eff_list = np.append(eff_list,fits_list[i]+'   '+'EFF_AB')
        if type_list[i] == 'FLUX,STD,SKY' and binx_list[i] == binx and biny_list[i] == biny:
            std_list = np.append(std_list,fits_list[i]+'   '+'FLUX')
        if type_list[i] == object_keyword  and binx_list[i] == binx and biny_list[i] == biny:
            sci_list = np.append(sci_list,fits_list[i]+'   '+object_tag)


    #The following checks that all these types are populated.
    if len(bias_list) == 0:
        print('ERROR: No BIAS frames detected. Check that you downloaded them properly.')
        sys.exit()
    if len(dark_list) == 0:
        print('ERROR: No DARK frames detected. Check that you downloaded them properly.')
        sys.exit()
    if len(LED_list) == 0:
        print('ERROR: No LED frames detected. Check that you downloaded them properly.')
        sys.exit()
    if len(orderdefA_list) == 0:
        print('ERROR: No ORDERDEF,LAMP,OFF frames detected. Check that you downloaded them properly.')
        sys.exit()
    if len(orderdefB_list) == 0:
        print('ERROR: No ORDERDEF,OFF,LAMP frames detected. Check that you downloaded them properly.')
        sys.exit()
    if len(flatA_list) == 0:
        print('ERROR: No FLAT,LAMP,OFF frames detected. Check that you downloaded them properly.')
        sys.exit()
    if len(flatB_list) == 0:
        print('ERROR: No FLAT,OFF,LAMP frames detected. Check that you downloaded them properly.')
        sys.exit()
    if len(FP_FP_list) == 0:
        print('ERROR: No WAVE,FP,FP frames detected. Check that you downloaded them properly.')
        sys.exit()
    if len(FP_TH_list) == 0:
        print('ERROR: No WAVE,FP,THAR frames detected. Check that you downloaded them properly.')
        sys.exit()
    if len(TH_FP_list) == 0:
        print('ERROR: No WAVE,THAR,FP frames detected. Check that you downloaded them properly.')
        sys.exit()
    if len(contam_list) == 0:
        print('ERROR: No CONTAM,OFF,FP frames detected. Check that you downloaded them properly.')
        sys.exit()
    if len(contam_list) > 1:
        print('ERROR: More than one CONTAM,OFF,FP frames detected. Please remove so you have only one left. The files dected are:')
        print(contam_list)
        sys.exit()
    if len(eff_list) == 0:
        print('ERROR: No EFF,SKY,SKY frames detected. Check that you downloaded them properly.')
        sys.exit()
    if len(std_list) == 0:
        print('ERROR: No FLUX,STD,SKY frames detected. Check that you downloaded them properly.')
        sys.exit()
    if len(sci_list) == 0:
        print('ERROR: No '+object_keyword+' frames detected. Check that you downloaded them properly.')
        sys.exit()
    if len(std_list) >= 2:
        print("WARNING: There is more than 1 FLUX,STD,SKY frame. Only using the first one.")
    if 'CCD_GEOM' not in static_dict.keys():
        print("ERROR: CCD_GEOM is missing. Was it downloaded correctly by the calselector?")
        sys.exit()
    if 'INST_CONFIG' not in static_dict.keys():
        print("ERROR: CCD_GEOM is missing. Was it downloaded correctly by the calselector?")
        sys.exit()




    #==============================================================================================#
    #The following blocks write out the SOF files for the different frame types.
    #Big wall of text with a lot of repetition but the functionality should be clear.
    #All of this comes from the ESPRESSO Pipeline manual v1.2.2, modulo the typos
    #that exist in there (some files are named differently in reality versus what's
    #written in the manual. Sorting this out once and for all is the purpose of this script.
    #Before the recipes are called, there is another script that checks whether the required files
    #as written in the SOF files are actually present.
    #==============================================================================================#

    outF = open(outpath+"BIAS.txt", "w")
    for line in bias_list:
        outF.write(line)
        outF.write("\n")
    import pdb
    pdb.set_trace()
    outF.write(static_dict['CCD_GEOM']+'   CCD_GEOM')
    outF.write("\n")
    outF.write(static_dict['INST_CONFIG']+'   INST_CONFIG')
    outF.close()

    outF = open(outpath+"DARK.txt", "w")
    for line in dark_list:
        outF.write(line)
        outF.write("\n")
    outF.write(static_dict['CCD_GEOM']+'   CCD_GEOM')
    outF.write("\n")
    outF.write(static_dict['INST_CONFIG']+'   INST_CONFIG')
    outF.write("\n")
    outF.write(path_out+'ESPRESSO_master_bias_res.fits MASTER_BIAS_RES')
    outF.close()

    outF = open(outpath+"LED.txt", "w")
    for line in LED_list:
        outF.write(line)
        outF.write("\n")
    outF.write(static_dict['CCD_GEOM']+'   CCD_GEOM')
    outF.write("\n")
    outF.write(static_dict['INST_CONFIG']+'   INST_CONFIG')
    outF.write("\n")
    outF.write(static_dict['LED_FF_GAIN_WINDOWS']+'   LED_FF_GAIN_WINDOWS')
    outF.write("\n")
    outF.write(path_out+'ESPRESSO_master_bias_res.fits MASTER_BIAS_RES')
    outF.write("\n")
    outF.write(path_out+'ESPRESSO_hot_pixels.fits HOT_PIXEL_MASK')
    outF.close()

    outF = open(outpath+"ORDERDEF.txt", "w")
    for line in orderdefA_list:
        outF.write(line)
        outF.write("\n")
    for line in orderdefB_list:
        outF.write(line)
        outF.write("\n")
    outF.write(static_dict['CCD_GEOM']+'   CCD_GEOM')
    outF.write("\n")
    outF.write(static_dict['INST_CONFIG']+'   INST_CONFIG')
    outF.write("\n")
    outF.write(path_out+'ESPRESSO_master_bias_res.fits MASTER_BIAS_RES')
    outF.write("\n")
    outF.write(path_out+'ESPRESSO_hot_pixels.fits HOT_PIXEL_MASK')
    outF.write("\n")
    outF.write(path_out+'ESPRESSO_bad_pixels.fits BAD_PIXEL_MASK')
    outF.close()


    outF = open(outpath+"FLAT.txt", "w")
    for line in flatA_list:
        outF.write(line)
        outF.write("\n")
    for line in flatB_list:
        outF.write(line)
        outF.write("\n")
    outF.write(static_dict['CCD_GEOM']+'   CCD_GEOM')
    outF.write("\n")
    outF.write(static_dict['INST_CONFIG']+'   INST_CONFIG')
    outF.write("\n")
    outF.write(path_out+'ESPRESSO_hot_pixels.fits HOT_PIXEL_MASK')
    outF.write("\n")
    outF.write(path_out+'ESPRESSO_bad_pixels.fits BAD_PIXEL_MASK')
    outF.write("\n")
    outF.write(path_out+'ESPRESSO_master_bias_res.fits MASTER_BIAS_RES')
    outF.write("\n")
    outF.write(path_out+'ESPRESSO_ORDER_TABLE_A.fits ORDER_TABLE_A')
    outF.write("\n")
    outF.write(path_out+'ESPRESSO_ORDER_TABLE_B.fits ORDER_TABLE_B')
    outF.write("\n")
    outF.write(static_dict['STATIC_WAVE_MATRIX_A']+'   STATIC_WAVE_MATRIX_A')
    outF.write("\n")
    outF.write(static_dict['STATIC_WAVE_MATRIX_B']+'   STATIC_WAVE_MATRIX_B')
    outF.close()

    outF = open(outpath+"WAVE_FP_FP.txt", "w")
    for line in FP_FP_list:
        outF.write(line)
        outF.write("\n")
    outF.write(static_dict['CCD_GEOM']+'   CCD_GEOM')
    outF.write("\n")
    outF.write(static_dict['INST_CONFIG']+'   INST_CONFIG')
    outF.write("\n")
    outF.write(path_out+'ESPRESSO_master_bias_res.fits MASTER_BIAS_RES')
    outF.write("\n")
    outF.write(path_out+'ESPRESSO_hot_pixels.fits HOT_PIXEL_MASK')
    outF.write("\n")
    outF.write(path_out+'ESPRESSO_bad_pixels.fits BAD_PIXEL_MASK')
    outF.write("\n")
    outF.write(path_out+'ESPRESSO_ORDER_TABLE_A.fits ORDER_TABLE_A')
    outF.write("\n")
    outF.write(path_out+'ESPRESSO_ORDER_TABLE_B.fits ORDER_TABLE_B')
    outF.write("\n")
    outF.write(path_out+'ESPRESSO_BLAZE_A.fits BLAZE_A')
    outF.write("\n")
    outF.write(path_out+'ESPRESSO_BLAZE_B.fits BLAZE_B')
    outF.write("\n")
    outF.write(path_out+'ESPRESSO_FLAT_A.fits FSPECTRUM_A')
    outF.write("\n")
    outF.write(path_out+'ESPRESSO_FLAT_B.fits FSPECTRUM_B')
    outF.write("\n")
    outF.write(path_out+'ESPRESSO_ORDER_PROFILE_A.fits ORDER_PROFILE_A')
    outF.write("\n")
    outF.write(path_out+'ESPRESSO_ORDER_PROFILE_B.fits ORDER_PROFILE_B')
    outF.close()


    outF = open(outpath+"WAVE_FP_TH.txt", "w")
    for line in FP_TH_list:
        outF.write(line)
        outF.write("\n")
    outF.write(static_dict['CCD_GEOM']+'   CCD_GEOM')
    outF.write("\n")
    outF.write(static_dict['INST_CONFIG']+'   INST_CONFIG')
    outF.write("\n")
    outF.write(path_out+'ESPRESSO_master_bias_res.fits MASTER_BIAS_RES')
    outF.write("\n")
    outF.write(path_out+'ESPRESSO_hot_pixels.fits HOT_PIXEL_MASK')
    outF.write("\n")
    outF.write(path_out+'ESPRESSO_bad_pixels.fits BAD_PIXEL_MASK')
    outF.write("\n")
    outF.write(path_out+'ESPRESSO_ORDER_TABLE_A.fits ORDER_TABLE_A')
    outF.write("\n")
    outF.write(path_out+'ESPRESSO_ORDER_TABLE_B.fits ORDER_TABLE_B')
    outF.write("\n")
    outF.write(path_out+'ESPRESSO_BLAZE_A.fits BLAZE_A')
    outF.write("\n")
    outF.write(path_out+'ESPRESSO_BLAZE_B.fits BLAZE_B')
    outF.write("\n")
    outF.write(path_out+'ESPRESSO_FLAT_A.fits FSPECTRUM_A')
    outF.write("\n")
    outF.write(path_out+'ESPRESSO_FLAT_B.fits FSPECTRUM_B')
    outF.write("\n")
    outF.write(path_out+'ESPRESSO_ORDER_PROFILE_A.fits ORDER_PROFILE_A')
    outF.write("\n")
    outF.write(path_out+'ESPRESSO_ORDER_PROFILE_B.fits ORDER_PROFILE_B')
    outF.write("\n")
    outF.write(static_dict['REF_LINE_TABLE_A']+'   REF_LINE_TABLE_A')
    outF.write("\n")
    outF.write(static_dict['REF_LINE_TABLE_B']+'   REF_LINE_TABLE_B')
    outF.write("\n")
    outF.write(path_out+'ESPRESSO_FP_SEARCHED_LINE_TABLE_A.fits FP_SEARCHED_LINE_TABLE_FP_FP_A')
    outF.write("\n")
    outF.write(path_out+'ESPRESSO_FP_SEARCHED_LINE_TABLE_B.fits FP_SEARCHED_LINE_TABLE_FP_FP_B')
    outF.write("\n")
    outF.write(path_out+'ESPRESSO_S2D_BLAZE_FP_FP_A.fits S2D_BLAZE_FP_FP_A')
    outF.write("\n")
    outF.write(path_out+'ESPRESSO_S2D_BLAZE_FP_FP_B.fits S2D_BLAZE_FP_FP_B')
    outF.write("\n")
    outF.write(static_dict['STATIC_DLL_MATRIX_A']+'   STATIC_DLL_MATRIX_A')
    outF.write("\n")
    outF.write(static_dict['STATIC_DLL_MATRIX_B']+'   STATIC_DLL_MATRIX_B')
    outF.write("\n")
    outF.write(static_dict['STATIC_WAVE_MATRIX_A']+'   STATIC_WAVE_MATRIX_A')
    outF.write("\n")
    outF.write(static_dict['STATIC_WAVE_MATRIX_B']+'   STATIC_WAVE_MATRIX_B')
    outF.close()


    outF = open(outpath+"WAVE_TH_FP.txt", "w")
    for line in TH_FP_list:
        outF.write(line)
        outF.write("\n")
    outF.write(static_dict['CCD_GEOM']+'   CCD_GEOM')
    outF.write("\n")
    outF.write(static_dict['INST_CONFIG']+'   INST_CONFIG')
    outF.write("\n")
    outF.write(path_out+'ESPRESSO_master_bias_res.fits MASTER_BIAS_RES')
    outF.write("\n")
    outF.write(path_out+'ESPRESSO_hot_pixels.fits HOT_PIXEL_MASK')
    outF.write("\n")
    outF.write(path_out+'ESPRESSO_bad_pixels.fits BAD_PIXEL_MASK')
    outF.write("\n")
    outF.write(path_out+'ESPRESSO_ORDER_TABLE_A.fits ORDER_TABLE_A')
    outF.write("\n")
    outF.write(path_out+'ESPRESSO_ORDER_TABLE_B.fits ORDER_TABLE_B')
    outF.write("\n")
    outF.write(path_out+'ESPRESSO_BLAZE_A.fits BLAZE_A')
    outF.write("\n")
    outF.write(path_out+'ESPRESSO_BLAZE_B.fits BLAZE_B')
    outF.write("\n")
    outF.write(path_out+'ESPRESSO_FLAT_A.fits FSPECTRUM_A')
    outF.write("\n")
    outF.write(path_out+'ESPRESSO_FLAT_B.fits FSPECTRUM_B')
    outF.write("\n")
    outF.write(path_out+'ESPRESSO_ORDER_PROFILE_A.fits ORDER_PROFILE_A')
    outF.write("\n")
    outF.write(path_out+'ESPRESSO_ORDER_PROFILE_B.fits ORDER_PROFILE_B')
    outF.write("\n")
    outF.write(static_dict['REF_LINE_TABLE_A']+'   REF_LINE_TABLE_A')
    outF.write("\n")
    outF.write(static_dict['REF_LINE_TABLE_B']+'   REF_LINE_TABLE_B')
    outF.write("\n")
    outF.write(path_out+'ESPRESSO_FP_SEARCHED_LINE_TABLE_A.fits FP_SEARCHED_LINE_TABLE_FP_FP_A')
    outF.write("\n")
    outF.write(path_out+'ESPRESSO_FP_SEARCHED_LINE_TABLE_B.fits FP_SEARCHED_LINE_TABLE_FP_FP_B')
    outF.write("\n")
    outF.write(path_out+'ESPRESSO_S2D_BLAZE_FP_FP_A.fits S2D_BLAZE_FP_FP_A')
    outF.write("\n")
    outF.write(path_out+'ESPRESSO_S2D_BLAZE_FP_FP_B.fits S2D_BLAZE_FP_FP_B')
    outF.write("\n")
    outF.write(static_dict['STATIC_DLL_MATRIX_A']+'   STATIC_DLL_MATRIX_A')
    outF.write("\n")
    outF.write(static_dict['STATIC_DLL_MATRIX_B']+'   STATIC_DLL_MATRIX_B')
    outF.write("\n")
    outF.write(static_dict['STATIC_WAVE_MATRIX_A']+'   STATIC_WAVE_MATRIX_A')
    outF.write("\n")
    outF.write(static_dict['STATIC_WAVE_MATRIX_B']+'   STATIC_WAVE_MATRIX_B')
    outF.close()


    outF = open(outpath+"CONTAM.txt", "w")
    for line in contam_list:
        outF.write(line)
        outF.write("\n")
    outF.write(static_dict['CCD_GEOM']+'   CCD_GEOM')
    outF.write("\n")
    outF.write(static_dict['INST_CONFIG']+'   INST_CONFIG')
    outF.write("\n")
    outF.write(path_out+'ESPRESSO_master_bias_res.fits MASTER_BIAS_RES')
    outF.write("\n")
    outF.write(path_out+'ESPRESSO_hot_pixels.fits HOT_PIXEL_MASK')
    outF.write("\n")
    outF.write(path_out+'ESPRESSO_bad_pixels.fits BAD_PIXEL_MASK')
    outF.write("\n")
    outF.write(path_out+'ESPRESSO_ORDER_TABLE_A.fits ORDER_TABLE_A')
    outF.write("\n")
    outF.write(path_out+'ESPRESSO_ORDER_TABLE_B.fits ORDER_TABLE_B')
    outF.write("\n")
    outF.write(path_out+'ESPRESSO_ORDER_PROFILE_A.fits ORDER_PROFILE_A')
    outF.write("\n")
    outF.write(path_out+'ESPRESSO_ORDER_PROFILE_B.fits ORDER_PROFILE_B')
    outF.write("\n")
    outF.write(path_out+'ESPRESSO_FLAT_A.fits FSPECTRUM_A')
    outF.write("\n")
    outF.write(path_out+'ESPRESSO_FLAT_B.fits FSPECTRUM_B')
    outF.close()

    outF = open(outpath+"EFF_SKY.txt", "w")
    #for line in eff_list:
        #outF.write(line)
        #outF.write("\n")
    outF.write(eff_list[0])
    outF.write("\n")
    outF.write(static_dict['CCD_GEOM']+'   CCD_GEOM')
    outF.write("\n")
    outF.write(static_dict['INST_CONFIG']+'   INST_CONFIG')
    outF.write("\n")
    outF.write(path_out+'ESPRESSO_hot_pixels.fits HOT_PIXEL_MASK')
    outF.write("\n")
    outF.write(path_out+'ESPRESSO_bad_pixels.fits BAD_PIXEL_MASK')
    outF.write("\n")
    outF.write(path_out+'ESPRESSO_ORDER_TABLE_A.fits ORDER_TABLE_A')
    outF.write("\n")
    outF.write(path_out+'ESPRESSO_ORDER_TABLE_B.fits ORDER_TABLE_B')
    outF.write("\n")
    outF.write(path_out+'ESPRESSO_ORDER_PROFILE_A.fits ORDER_PROFILE_A')
    outF.write("\n")
    outF.write(path_out+'ESPRESSO_ORDER_PROFILE_B.fits ORDER_PROFILE_B')
    outF.write("\n")
    outF.write(path_out+'ESPRESSO_FLAT_A.fits FSPECTRUM_A')
    outF.write("\n")
    outF.write(path_out+'ESPRESSO_FLAT_B.fits FSPECTRUM_B')
    outF.close()


    outF = open(outpath+"FLUX_STD.txt", "w")
    if len(std_list) >= 1:
        outF.write(std_list[0])
        outF.write("\n")
    outF.write(static_dict['CCD_GEOM']+'   CCD_GEOM')
    outF.write("\n")
    outF.write(static_dict['INST_CONFIG']+'   INST_CONFIG')
    outF.write("\n")
    outF.write(static_dict['STD_TABLE']+'   STD_TABLE')
    outF.write("\n")
    outF.write(static_dict['EXT_TABLE']+'   EXT_TABLE')
    outF.write("\n")
    outF.write(path_out+'ESPRESSO_master_bias_res.fits MASTER_BIAS_RES')
    outF.write("\n")
    outF.write(path_out+'ESPRESSO_hot_pixels.fits HOT_PIXEL_MASK')
    outF.write("\n")
    outF.write(path_out+'ESPRESSO_bad_pixels.fits BAD_PIXEL_MASK')
    outF.write("\n")
    outF.write(path_out+'ESPRESSO_ORDER_TABLE_A.fits ORDER_TABLE_A')
    outF.write("\n")
    outF.write(path_out+'ESPRESSO_ORDER_TABLE_B.fits ORDER_TABLE_B')
    outF.write("\n")
    outF.write(path_out+'ESPRESSO_ORDER_PROFILE_A.fits ORDER_PROFILE_A')
    outF.write("\n")
    outF.write(path_out+'ESPRESSO_ORDER_PROFILE_B.fits ORDER_PROFILE_B')
    outF.write("\n")
    outF.write(path_out+'ESPRESSO_FLAT_A.fits FSPECTRUM_A')
    outF.write("\n")
    outF.write(path_out+'ESPRESSO_FLAT_B.fits FSPECTRUM_B')
    outF.write("\n")
    outF.write(path_out+'ESPRESSO_BLAZE_A.fits BLAZE_A')
    outF.write("\n")
    outF.write(path_out+'ESPRESSO_BLAZE_B.fits BLAZE_B')
    outF.write("\n")
    outF.write(path_out+'ESPRESSO_WAVE_MATRIX_A.fits WAVE_MATRIX_FP_THAR_B')
    outF.write("\n")
    outF.write(path_out+'ESPRESSO_WAVE_MATRIX_B.fits WAVE_MATRIX_THAR_FP_A')
    outF.close()




    #Almost there...
    outF = open(outpath+"SCI_OBJ_part1.txt", "w")
    outF.write('EMPTY LINE')
    for line in sci_list:
        outF.write(line)
        outF.write("\n")
    outF.close()

    outF = open(outpath+"SCI_OBJ_part2.txt", "w")
    outF.write(static_dict['CCD_GEOM']+'   CCD_GEOM')
    outF.write("\n")
    outF.write(static_dict['INST_CONFIG']+'   INST_CONFIG')
    outF.write("\n")
    outF.write(static_dict['EXT_TABLE']+'   EXT_TABLE')
    outF.write("\n")
    outF.write(static_dict['MASK_LUT']+'   MASK_LUT')
    outF.write("\n")
    for line in mask_list:
        outF.write(line+'   MASK_TABLE')
        outF.write("\n")
    outF.write(static_dict['STD_TABLE']+'   STD_TABLE')
    outF.write("\n")
    outF.write(path_out+'ESPRESSO_master_bias_res.fits MASTER_BIAS_RES')
    outF.write("\n")
    outF.write(path_out+'ESPRESSO_hot_pixels.fits HOT_PIXEL_MASK')
    outF.write("\n")
    outF.write(path_out+'ESPRESSO_bad_pixels.fits BAD_PIXEL_MASK')
    outF.write("\n")
    outF.write(path_out+'ESPRESSO_ORDER_TABLE_A.fits ORDER_TABLE_A')
    outF.write("\n")
    outF.write(path_out+'ESPRESSO_ORDER_TABLE_B.fits ORDER_TABLE_B')
    outF.write("\n")
    outF.write(path_out+'ESPRESSO_ORDER_PROFILE_A.fits ORDER_PROFILE_A')
    outF.write("\n")
    outF.write(path_out+'ESPRESSO_ORDER_PROFILE_B.fits ORDER_PROFILE_B')
    outF.write("\n")
    outF.write(path_out+'ESPRESSO_FLAT_A.fits FSPECTRUM_A')
    outF.write("\n")
    outF.write(path_out+'ESPRESSO_FLAT_B.fits FSPECTRUM_B')
    outF.write("\n")
    outF.write(path_out+'ESPRESSO_BLAZE_A.fits BLAZE_A')
    outF.write("\n")
    outF.write(path_out+'ESPRESSO_BLAZE_B.fits BLAZE_B')
    outF.write("\n")
    outF.write(path_out+'ESPRESSO_S2D_BLAZE_THAR_FP_A.fits S2D_BLAZE_THAR_FP_A')
    outF.write("\n")
    outF.write(path_out+'ESPRESSO_S2D_BLAZE_THAR_FP_B.fits S2D_BLAZE_THAR_FP_B')
    outF.write("\n")
    outF.write(path_out+'ESPRESSO_WAVE_MATRIX_A.fits WAVE_MATRIX_FP_THAR_B')
    outF.write("\n")
    outF.write(path_out+'ESPRESSO_WAVE_MATRIX_B.fits WAVE_MATRIX_THAR_FP_A')
    outF.write("\n")
    outF.write(path_out+'ESPRESSO_DLL_MATRIX_B.fits DLL_MATRIX_FP_THAR_B')
    outF.write("\n")
    outF.write(path_out+'ESPRESSO_DLL_MATRIX_A.fits DLL_MATRIX_THAR_FP_A')
    outF.write("\n")
    outF.write(static_dict['FLUX_TEMPLATE']+'   FLUX_TEMPLATE')
    outF.write("\n")
    outF.write(path_out+'ESPRESSO_CONTAM_FP_B.fits CONTAM_FP')
    outF.write("\n")
    outF.write(path_out+'ESPRESSO_REL_EFF_B.fits REL_EFF_B')
    outF.write("\n")
    outF.write(path_out+'ESPRESSO_ABS_EFF_A.fits ABS_EFF_A')
    outF.close()


    #==============================================================================================#
    #==============================================================================================#
    #This is the end of the create_sof script. What follows are three small utility functions, after
    #which there is a set of wrappers for the esorex recipes.
    #==============================================================================================#
    #==============================================================================================#






def move_to(filename,outpath,newname=None):
    import pdb
    """This short script moves a file at location filename to the folder outpath.
    If the newname keyword is set to a string, the file will also be renamed in the process.
    This moving overwrites existing files."""
    #I was lazy to type shutil.move all the time...
    import shutil
    if newname == None:
        shutil.move(filename,outpath+filename)
    else:
        shutil.move(filename,outpath+newname)


def clean_trash():
    """This program deletes left-over files at the end of a recipe."""
    import os
    os.system('rm -rf *fits *.log *.dump')


def check_files_exist(sof_file):
    "This program reads the sof file prior to execution of the recipe, to make sure that all the dependent files actually exist. This is to prevent the recipe running for 3 hours and then crashing due to a missing file or a wrongly spelled filename somewhere. If the tag is spelled wrongly, well then hopefully the recipe itself will crash at the start."""
    import csv
    import os
    import sys


    f=open(sof_file,'r').read().splitlines()
    for line in f:
        if len(line.split()) > 2:#If there are spaces in the main path (DONT DO THIS) then we split on the .fits extension instead.
            filename = line.split('.fits')[0]+'.fits'
        else:
            filename=line.split()[0]
        exists=os.path.isfile(filename)
        if exists != True:
            print("ERROR IN RECIPE PATH FILE: "+filename+" is required for "+sof_file+" but it doesn't exist. Check:")
            print("  1. Whether this is a static calibration file in your esorex installation that is named wrongly. (Could happen if those building the pipeline have changed the name of their static calibration files in a new version).")
            print("  2. That all required previous recipes were executed.")
            print("  3. That previous recipes produced the right output files, and that these were moved to the right (outpath) folder.")
            print("  4. That the recipe path file file was created correctly by create_sof (i.e. without typos).")
            sys.exit()





    #==============================================================================================#
    #==============================================================================================#
    #What follows are the wrappers for the esorex recipes. These are executed one by one when
    #calling this script, all the way at the end of this file.
    #==============================================================================================#
    #==============================================================================================#






def master_bias(outpath):
    """This is a wrapper for the mbias recipe."""
    import os
    print('==========>>>>> CREATING MASTER BIAS<<<<<==========')
    check_files_exist(outpath+"BIAS.txt")
    os.system("esorex espdr_mbias "+outpath+"BIAS.txt")
    move_to('ESPRESSO_master_bias.fits',outpath)
    move_to('esorex.log',outpath,newname='esorex_masterbias.log')
    move_to('ESPRESSO_master_bias_res.fits',outpath)


def master_dark(outpath):
    """This is a wrapper for the mdark recipe."""
    import os
    print('==========>>>>> CREATING MASTER DARK AND HOT PIXEL MAP<<<<<==========')
    check_files_exist(outpath+"DARK.txt")
    os.system("esorex espdr_mdark "+outpath+"DARK.txt")
    move_to('ESPRESSO_master_dark.fits',outpath)
    move_to('ESPRESSO_hot_pixels.fits',outpath)
    move_to('esorex.log',outpath,newname='esorex_masterdark.log')
    clean_trash()


def bad_pixels(outpath):
    """This is a wrapper for the led_ff recipe."""
    import os
    print('==========>>>>> CREATING BAD PIXEL MAP<<<<<==========')
    check_files_exist(outpath+"LED.txt")
    os.system("esorex espdr_led_ff "+outpath+"LED.txt")
    move_to('ESPRESSO_bad_pixels.fits',outpath)
    move_to('esorex.log',outpath,newname='esorex_badpixels.log')
    clean_trash()


def orderdef(outpath):
    """This is a wrapper for the orderdef recipe."""
    import os
    print('==========>>>>> FIND ORDER TRACES<<<<<==========')
    check_files_exist(outpath+"ORDERDEF.txt")
    os.system("esorex espdr_orderdef "+outpath+"ORDERDEF.txt")
    move_to('ESPRESSO_ORDER_TABLE_A.fits',outpath)
    move_to('ESPRESSO_ORDER_TABLE_B.fits',outpath)
    move_to('esorex.log',outpath,newname='esorex_orderdef.log')
    clean_trash()

def master_flat(outpath):
    """This is a wrapper for the mflat recipe."""
    import os
    print('==========>>>>> CREATE MASTER FLAT<<<<<==========')
    check_files_exist(outpath+"FLAT.txt")
    os.system("esorex espdr_mflat "+outpath+"FLAT.txt")
    move_to('ESPRESSO_ORDER_PROFILE_A.fits',outpath)
    move_to('ESPRESSO_ORDER_PROFILE_B.fits',outpath)
    move_to('ESPRESSO_BLAZE_A.fits',outpath)
    move_to('ESPRESSO_BLAZE_B.fits',outpath)
    move_to('ESPRESSO_FLAT_A.fits',outpath)
    move_to('ESPRESSO_FLAT_B.fits',outpath)
    move_to('ESPRESSO_background_map_A.fits',outpath)
    move_to('ESPRESSO_background_map_B.fits',outpath)
    move_to('ESPRESSO_spectrum_extracted_A.fits',outpath)
    move_to('ESPRESSO_spectrum_extracted_B.fits',outpath)
    move_to('esorex.log',outpath,newname='esorex_mflat.log')
    clean_trash()


def wave_FP_FP(outpath):
    """This is a wrapper for the wave_FP_FP recipe."""
    import os
    print('==========>>>>> CREATE WAVE FP_FP <<<<<==========')
    check_files_exist(outpath+'WAVE_FP_FP.txt')
    os.system("esorex espdr_wave_FP "+outpath+"WAVE_FP_FP.txt")
    move_to('ESPRESSO_S2D_FP_FP_A.fits',outpath)
    move_to('ESPRESSO_S2D_FP_FP_B.fits',outpath)
    move_to('ESPRESSO_S2D_BLAZE_FP_FP_A.fits',outpath)
    move_to('ESPRESSO_S2D_BLAZE_FP_FP_B.fits',outpath)
    move_to('ESPRESSO_FP_SEARCHED_LINE_TABLE_A.fits',outpath)
    move_to('ESPRESSO_FP_SEARCHED_LINE_TABLE_B.fits',outpath)
    move_to('esorex.log',outpath,newname='esorex_wave_fp_fp.log')
    clean_trash


def wave_FP_TH(outpath):
    """This is a wrapper for the wave_FP_THAR recipe."""
    import os
    print('==========>>>>> CREATE WAVE FP_THAR<<<<<==========')
    check_files_exist(outpath+'WAVE_FP_TH.txt')
    os.system("esorex espdr_wave_THAR "+outpath+"WAVE_FP_TH.txt")
    move_to('ESPRESSO_AIR_DLL_MATRIX_B.fits',outpath)
    move_to('ESPRESSO_AIR_WAVE_MATRIX_B.fits',outpath)
    move_to('ESPRESSO_DLL_MATRIX_B.fits',outpath)
    move_to('ESPRESSO_FP_FITTED_LINE_TABLE_B.fits',outpath)
    move_to('ESPRESSO_LINE_TABLE_RAW_B.fits',outpath)
    move_to('ESPRESSO_S2D_BLAZE_FP_THAR_A.fits',outpath)
    move_to('ESPRESSO_S2D_BLAZE_FP_THAR_B.fits',outpath)
    move_to('ESPRESSO_S2D_FP_THAR_A.fits',outpath)
    move_to('ESPRESSO_S2D_FP_THAR_B.fits',outpath)
    move_to('ESPRESSO_WAVE_MATRIX_B.fits',outpath)
    move_to('ESPRESSO_WAVE_TABLE_B.fits',outpath)
    move_to('ESPRESSO_THAR_LINE_TABLE_B.fits',outpath)
    move_to('esorex.log',outpath,newname='esorex_wave_fp_thar.log')
    clean_trash

def wave_TH_FP(outpath):
    """This is a wrapper for the wave_THAR_FP recipe."""
    import os
    print('==========>>>>> CREATE WAVE THAR_FP<<<<<==========')
    check_files_exist(outpath+'WAVE_TH_FP.txt')
    os.system("esorex espdr_wave_THAR "+outpath+"WAVE_TH_FP.txt")
    move_to('ESPRESSO_AIR_DLL_MATRIX_A.fits',outpath)
    move_to('ESPRESSO_AIR_WAVE_MATRIX_A.fits',outpath)
    move_to('ESPRESSO_DLL_MATRIX_A.fits',outpath)
    move_to('ESPRESSO_FP_FITTED_LINE_TABLE_A.fits',outpath)
    move_to('ESPRESSO_LINE_TABLE_RAW_A.fits',outpath)
    move_to('ESPRESSO_S2D_BLAZE_THAR_FP_A.fits',outpath)
    move_to('ESPRESSO_S2D_BLAZE_THAR_FP_B.fits',outpath)
    move_to('ESPRESSO_S2D_THAR_FP_A.fits',outpath)
    move_to('ESPRESSO_S2D_THAR_FP_B.fits',outpath)
    move_to('ESPRESSO_WAVE_MATRIX_A.fits',outpath)
    move_to('ESPRESSO_WAVE_TABLE_A.fits',outpath)
    move_to('ESPRESSO_THAR_LINE_TABLE_A.fits',outpath)
    move_to('esorex.log',outpath,newname='esorex_wave_fp_thar.log')
    clean_trash()

def contamination(outpath):
    """This is a wrapper for the contam recipe."""
    import os
    print('==========>>>>> CREATE CROSS-FIBER CONTAMINATION FRAMES <<<<<==========')
    check_files_exist(outpath+'CONTAM.txt')
    os.system('esorex espdr_cal_contam '+outpath+'CONTAM.txt')
    move_to('ESPRESSO_CONTAM_FP_B.fits',outpath)
    move_to('ESPRESSO_CONTAM_S2D_A.fits',outpath)
    move_to('ESPRESSO_CONTAM_S2D_B.fits',outpath)
    move_to('esorex.log',outpath,newname='esorex_cal_contam.log')
    clean_trash()

def relative_efficiency(outpath):
    """This is a wrapper for the eff_ab recipe."""
    import os
    print('==========>>>>> CREATE RELATIVE FIBER EFFICIENCY FRAMES <<<<<==========')
    check_files_exist(outpath+"EFF_SKY.txt")
    os.system('esorex espdr_cal_eff_ab '+outpath+"EFF_SKY.txt")
    move_to('ESPRESSO_S2D_BLAZE_EFF_A.fits',outpath)
    move_to('ESPRESSO_S2D_BLAZE_EFF_B.fits',outpath)
    move_to('ESPRESSO_REL_EFF_B.fits',outpath)
    move_to('esorex.log',outpath,newname='esorex_cal_eff_ab.log')
    clean_trash()

def flux_calibration(outpath):
    """This is a wrapper for the  recipe."""
    import os
    print('==========>>>>> CREATE FLUX CALIBRATION FRAMES <<<<<==========')
    check_files_exist(outpath+"FLUX_STD.txt")
    os.system('esorex espdr_cal_flux '+outpath+"FLUX_STD.txt")
    move_to('ESPRESSO_S2D_STD_A.fits',outpath)
    move_to('ESPRESSO_S1D_STD_A.fits',outpath)
    move_to('ESPRESSO_S1D_ENERGY_STD_A.fits',outpath)
    move_to('ESPRESSO_S2D_BLAZE_STD_A.fits',outpath)
    move_to('ESPRESSO_AVG_FLUX_STD_A.fits',outpath)
    move_to('ESPRESSO_ABS_EFF_RAW_A.fits',outpath)
    move_to('ESPRESSO_ABS_EFF_A.fits',outpath)
    move_to('esorex.log',outpath,newname='esorex_cal_flux.log')
    clean_trash()

def reduce_science(outpath):
    import os
    import astropy.io.ascii as ascii
    import pdb
    import shutil

    print('==========>>>>> PRODUCE REDUCED SCIENCE SPECTRA <<<<<==========')
#    check_files_exist(outpath+'SCI_OBJ_part1.txt')
    check_files_exist(outpath+'SCI_OBJ_part2.txt')

    F=ascii.read(outpath+'SCI_OBJ_part1.txt',names=['paths','tags'])
    N=len(F['paths'])
    if not os.path.exists(outpath+'SCIENCE_PRODUCTS'):
        os.mkdir(outpath+'SCIENCE_PRODUCTS')

    for i in range(N):
        shutil.copy(outpath+'SCI_OBJ_part2.txt',outpath+'SCI_OBJ_combined.txt')
        filename=os.path.splitext(os.path.basename(F['paths'][i]))[0]
        with open(outpath+'SCI_OBJ_combined.txt','a') as SOF:
            SOF.write('\n')
            SOF.write(F['paths'][i]+' '+F['tags'][i])
        #pdb.set_trace()
        print('>>>> RUNNING FILE '+F['paths'][i])
        os.system('esorex espdr_sci_red --background_sw=off '+outpath+'SCI_OBJ_combined.txt')
        move_to('ESPRESSO_CCF_A.fits',outpath+'SCIENCE_PRODUCTS/',newname=filename+'_CCF_A.fits')
        move_to('ESPRESSO_CCF_B.fits',outpath+'SCIENCE_PRODUCTS/',newname=filename+'_CCF_B.fits')
        move_to('ESPRESSO_CCF_RESIDUALS_A.fits',outpath+'SCIENCE_PRODUCTS/',newname=filename+'_CCF_RESIDUALS_A.fits')
        move_to('ESPRESSO_CCF_SKYSUB_A.fits',outpath+'SCIENCE_PRODUCTS/',newname=filename+'_CCF_SKYSUB_A.fits')
        move_to('ESPRESSO_S1D_A.fits',outpath+'SCIENCE_PRODUCTS/',newname=filename+'_S1D_A.fits')
        move_to('ESPRESSO_S1D_B.fits',outpath+'SCIENCE_PRODUCTS/',newname=filename+'_S1D_B.fits')
        move_to('ESPRESSO_S1D_FLUXCAL_A.fits',outpath+'SCIENCE_PRODUCTS/',newname=filename+'_S1D_FLUXCAL_A.fits')
        move_to('ESPRESSO_S1D_SKYSUB_A.fits',outpath+'SCIENCE_PRODUCTS/',newname=filename+'_S1D_SKYSUB_A.fits')
        move_to('ESPRESSO_S1D_SKYSUB_FLUXCAL_A.fits',outpath+'SCIENCE_PRODUCTS/',newname=filename+'_S1D_SKYSUB_FLUXCAL_A.fits')
        move_to('ESPRESSO_S2D_A.fits',outpath+'SCIENCE_PRODUCTS/',newname=filename+'_S2D_A.fits')
        move_to('ESPRESSO_S2D_B.fits',outpath+'SCIENCE_PRODUCTS/',newname=filename+'_S2D_B.fits')
        move_to('ESPRESSO_S2D_BLAZE_A.fits',outpath+'SCIENCE_PRODUCTS/',newname=filename+'_S2D_BLAZE_A.fits')
        move_to('ESPRESSO_S2D_BLAZE_B.fits',outpath+'SCIENCE_PRODUCTS/',newname=filename+'_S2D_BLAZE_B.fits')
        move_to('ESPRESSO_S2D_SKYSUB_A.fits',outpath+'SCIENCE_PRODUCTS/',newname=filename+'_S2D_SKYSUB_A.fits')

        if sky:#The following files dont exist if spectra were taken with the FP on fiber B:
            move_to('ESPRESSO_CCF_B.fits',outpath+'SCIENCE_PRODUCTS/',newname=filename+'_CCF_B.fits')
            move_to('ESPRESSO_CCF_SKYSUB_A.fits',outpath+'SCIENCE_PRODUCTS/',newname=filename+'_CCF_SKYSUB_A.fits')
            move_to('ESPRESSO_S2D_SKYSUB_A.fits',outpath+'SCIENCE_PRODUCTS/',newname=filename+'_S2D_SKYSUB_A.fits')
            move_to('ESPRESSO_S1D_SKYSUB_A.fits',outpath+'SCIENCE_PRODUCTS/',newname=filename+'_S1D_SKYSUB_A.fits')
            move_to('ESPRESSO_S1D_SKYSUB_FLUXCAL_A.fits',outpath+'SCIENCE_PRODUCTS/',newname=filename+'_S1D_SKYSUB_FLUXCAL_A.fits')














inpath='/run/media/jens/SAMSUNG/Raw_data/data_with_raw_calibs/'
outpath='/run/media/jens/10AC118E10AC118E/Ilm/Echelle_reduction/ESPRESSO/reduced/'
path_cdb='/home/jens/ESO_PIPELINES/calib/espdr-2.2.1/'
binning='2x1'

create_sof(inpath,outpath,path_cdb,binning)
master_bias(outpath)
master_dark(outpath)
bad_pixels(outpath)
orderdef(outpath)
master_flat(outpath)
wave_FP_FP(outpath)
wave_FP_TH(outpath)
wave_TH_FP(outpath)
contamination(outpath)
relative_efficiency(outpath)
flux_calibration(outpath)
reduce_science(outpath)
