#This code relies on asto_script to be INSIDE the "QUESTdata" folder, next to all of
#the folders (darks,flats,science).

#The aim of the code is to extract data from fits files found inside their
#corresponding folders and do image correction until the user is left with
#sky subtracted images that are also x and y adjusted

import numpy as np
from astropy.io import fits
from scipy.stats import mode
from scipy.ndimage import interpolation
import matplotlib.pyplot as plt
import os

###Some very useful functions to cut down on coding and computational time###

def file_name_nohidden(folder):
    """
    Inputs: folder - str; The folder where the modified listdir will operate

    Returns: lst   - list of str; All valid file names in a list.

    This function is a modified version of the os.listdir method which excludes
    fits files starting with '.' (These exist due to .tar file extraction)
    """
    lst = []
    for f in os.listdir(folder):
        if not f.startswith('.') and f.endswith('.fits'):
            lst.append(f)
    return lst

def exposure_check(dark_folder, flat_folder, science_folder):
    """
    Inputs: dark_folder - str; The folder where dark files are located

            flat_folder - str; The folder where flat files are located

            science_folder - str; The folder where science files are located

    Returns: Prints out the exposure times of the files found in the given
             folders.
    """
    dark_files = file_name_nohidden(dark_folder)
    flat_files = file_name_nohidden(flat_folder)
    science_files = file_name_nohidden(science_folder)

    for file in dark_files:
        dark_read = fits.getheader(dark_folder+'/'+file)
        print('Dark file',file,'has an exposure of',dark_read['EXPTIME'],'s')

    for file in flat_files:
        flat_read = fits.getheader(flat_folder+'/'+file)
        print('Flat file',file,'has exposure of',flat_read['EXPTIME'],'s')

    for file in science_files:
        science_read = fits.getheader(science_folder+'/'+file)
        print('science file',file,'has an exposure of',science_read['EXPTIME'],'s')

def bias_overscan(bias_folder,bias_file,overscan_start=None,overscan_end=None,\
                  orientation='vertical'):
    """
    Inputs: bias_folder - str; The folder of the file

            bias_file - str; The file name

            overscan_start - int(default value is None); The beginning of the overscan
                             region in pixels
            overscan_end   - int(default value is None); The end of the overscan
                             region in pixels

            orientation    - str(default value is 'vertical'); This instructs
                             whether the code should be looking at a vertical or
                             horizontal overscan strip.


    Returns: A bias reduced file in a seperate folder

    IMPORTANT: Assumes that the overscan region is at the right end of the image
    for 'vertical' and at the bottom for 'horizontal' for orientation varible

    Used when no bias files available, likely when a higher quality CCD was used
    """
    if overscan_start is None:
        return 'Overscan value error: Provide an overscan area to proceed'

    fits_data = fits.getdata(bias_folder+'/'+bias_file)
    fits_head = fits.getheader(bias_folder+'/'+bias_file)
    fits_final = fits_data
    overscan_med = np.array([]) #Container for median values

    #Looking at the columns in this case
    if orientation == 'vertical':
        if overscan_end is None:
            overscan_reg = fits_data[:,overscan_start:]
        else:
            overscan_reg = fits_data[:,overscan_start:overscan_end]
        for os_row in range(overscan_reg.shape[0]):
            bias_med_row = np.median(overscan_reg[os_row])
            overscan_med = np.append(overscan_med, bias_med_row)
            fits_final[os_row] = fits_final[os_row] - bias_med_row
        overscan_med = overscan_med[:, None] #Rotates the array to a 1D vertical array

    #Looking at the rows in this case
    elif orientation == 'horizontal':
        if overscan_end is None:
            overscan_reg = fits_data[overscan_start:,:]
        else:
            overscan_reg = fits_data[overscan_start:overscan_end,:]
        for os_col in range(overscan_reg.shape[1]):
            bias_med_col = np.median(overscan_reg[:,os_col])
            overscan_med = np.append(overscan_med, bias_med_col)
            fits_final[os_col] = fits_final[os_col] - bias_med_col
    else:
        return "Invalid orientation. Check spelling and use only 'vertical' and\
                'horizontal'."
    #New filtered fits image
    if os.path.isdir('bias_filtered/') is False:
        os.mkdir('bias_filtered/')
    fits.writeto(filename='bias_filtered/'+'bias_filt'+bias_file+'.fits',\
                 data=fits_final,header=fits_head,output_verify='silentfix',\
                 overwrite=True)


def sky_subtract(redux_folder):
    """
    Inputs: redux_folder - str; The folder where the reduced images are located

    Returns: Sky-subtracted fits files of the ones found in the folder selected.
    """
    redux_files = file_name_nohidden(redux_folder)
    for i,file in enumerate(redux_files):
        redux_raw = fits.getdata(redux_folder+'/'+file)
        redux_header = fits.getheader(redux_folder+'/'+file)
        skymode, _ = mode(redux_raw[(redux_raw>0) & (redux_raw<1e10)].flatten())
        redux_final  = redux_raw - skymode
        if os.path.isdir('skysub/') is False:
            os.mkdir('skysub/')
        fits.writeto(filename=('skysub/skysub_'+file),data=redux_final,\
                     header=redux_header,overwrite=True)


#==========================Script part starts here==============================

#Setting folder names and bias file naming list
dark_folder = 'darks'
flat_folder = 'flats'
science_folder = 'science'
bias_files_flats = file_name_nohidden(flat_folder)
bias_files_sci   = file_name_nohidden(science_folder)
redux_folder   = 'redux'

#Overscan region for both flats and science images
os_start = 599
os_end = 635

#Important for a overscan calculation function
orientation = 'vertical'

#Checks exposure
#exposure_check(dark_folder,flat_folder,science_folder)
#^^^^ flats have exposure of 10s, science have 180s ^^^^

#Overscan for all of flat files
for file in bias_files_flats:
    bias_overscan(bias_folder=flat_folder,bias_file=file,\
                     overscan_start=os_start,overscan_end=os_end,\
                     orientation=orientation)


#Remove the overscan region from all corresponding flats ^^^
dark_flats = fits.getdata('darks/dark_10.C22.fits')
flat_mod   = file_name_nohidden('bias_filtered')
size_findr = fits.getdata('bias_filtered/'+flat_mod[0])
big_flat_arr = np.empty((size_findr.shape[0],size_findr.shape[1],0))


for file in flat_mod:
    #Find import only necessary files (to avoid mixup with bias corrected science corrected)
    if file.endswith('e.C22.fits.fits') or file.endswith('m.C22.fits.fits'):
        flat_ext     = fits.getdata('bias_filtered/'+file)
        flat_ext     = flat_ext - dark_flats
        big_flat_arr = np.dstack((big_flat_arr,flat_ext))


#Median along the 3D axis to get a master flat file
flat_med = np.median(big_flat_arr,axis=2)


#Normalisation of master flat
flat_master = flat_med / np.median(flat_med.flatten())


#Removing redundant os region
flat_master = flat_master[:,:597]


#Flat pixel value histo creation
plt.hist(flat_master.flatten(),bins=100,facecolor='#2ab0ff',edgecolor='#e0e0e0')
plt.xlabel('Pixel value')
plt.ylabel('Frequency')
plt.xlim(0,2)
plt.show()

#Saving master flat
fits.writeto('flat_master.fits',data=flat_master,overwrite=True)

#Bias remove from science images (overscan region)
for file in bias_files_sci:
    bias_overscan(bias_folder=science_folder,bias_file=file,\
                     overscan_start=os_start,overscan_end=os_end,\
                     orientation=orientation)

#Flat field them and remove dark field
flat_master = fits.getdata('flat_master.fits')
dark_science = fits.getdata('darks/dark_180.C22.fits')
dark_science = dark_science[:,:597] #Removing the overscan region to keep the resolution consistent

if os.path.isdir('redux/') is False:
    os.mkdir('redux/')
for file in file_name_nohidden('bias_filtered'):
    #Filter out the bias corrected flat files
    if file.endswith('s.C22.fits.fits'):
        sci_ext     = fits.getdata('bias_filtered/'+file)
        sci_ext_hdu = fits.getheader('bias_filtered/'+file)
        sci_ext = sci_ext[:,:597] #Removing os region
        sci_mod = sci_ext - dark_science
        sci_mod = sci_mod / flat_master
        fits.writeto(filename='redux/'+'reduced'+file,data=sci_mod,\
                     header=sci_ext_hdu,overwrite=True)

#Sky reduce
sky_subtract(redux_folder)

#Correct x and y
ref_img = fits.getdata('skysub/skysub_reducedbias_filt20130910234901s.C22.fits.fits')
shift_img_1 = fits.getdata('skysub/skysub_reducedbias_filt20130911020246s.C22.fits.fits')
shift_head_1 = fits.getheader('skysub/skysub_reducedbias_filt20130911020246s.C22.fits.fits')
shift_adj = interpolation.shift(shift_img_1,(16,-4))
fits.writeto(filename='shift_adj246s.fits',data=shift_adj,header=shift_head_1,overwrite=True)

master_img = ref_img + shift_adj #Combined img creation
shift_img_2 = fits.getdata('skysub/skysub_reducedbias_filt20130911040543s.C22.fits.fits')
shift_head_2 = fits.getheader('skysub/skysub_reducedbias_filt20130911040543s.C22.fits.fits')
shift_adj = interpolation.shift(shift_img_2,(40,-6))
fits.writeto(filename='shift_adj543s.fits',data=shift_adj,header=shift_head_2,overwrite=True)

#Combined img finalisation
master_img += shift_adj
fits.writeto(filename='combined.fits',data=master_img,overwrite=True)
