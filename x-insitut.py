from __future__ import division
from lyse import *
from pylab import *
from time import time
from scipy.ndimage import *
from mpl_toolkits.axes_grid1 import make_axes_locatable
from common.OD_handler_fl import ODShot
import os
import pandas as pd
import numpy as np
import numexpr as ne
import matplotlib.gridspec as gridspec



# Time stamp
print '\nRunning %s' % os.path.basename(__file__)
t = time()

# Load dataframe
run = Run(path)

# Methods
def print_time(text):
    print 't = %6.3f : %s' % ((time()-t), text)
 
def raw_to_OD(fpath):
    with h5py.File(fpath) as h5_file:
        image_group = {}
        if '/data/imagesYZ_2_Flea3' in h5_file:
            image_group['imagesYZ_2_Flea3'] = 'x-MOT_Fluorescence'
            Isat = 100          
            alpha = 1 #1.645
        if  '/data/imagesXY_1_Flea3' in h5_file:
            image_group['imagesXY_1_Flea3'] = 'z-TOF'
            Isat_TOF = 297          # In counts // Paco:04/21/2016
            alpha = 1 #1.645
        if  '/data/imagesXY_2_Flea3' in h5_file:
            image_group['imagesXY_2_Flea3'] = 'z-insitu'
            Isat_insitu = 395      # In counts // Paco:04/20/2016
            alpha= 1 #4.3
        if '/data/imagesYZ_1_Flea3' in h5_file:
            image_group['imagesYZ_1_Flea3'] = 'x-insitu'
            Isat = 2248     # In counts // Paco:05/06/2016
            alpha = 1 #0.441
        atoms, probe, bckg, div, calculated_OD = {}, {}, {}, {}, {}
        for appended, id in image_group.iteritems():
            # Check if data is present
            if len(h5_file['data'][appended]['Raw'][:]) != 0:
                if id == 'x-MOT_Fluorescence':
                    MOT_Fluorescence_atoms1 = (np.array(h5_file['data'][appended]['Raw'])[0])
                    MOT_Fluorescence_dark = (np.array(h5_file['data'][appended]['Raw'])[1])
                    MOT_Fluorescence_atoms2 = (np.array(h5_file['data'][appended]['Raw'])[2])
                    
                    
                
                    
                    
        
        atoms1_mask = np.ma.masked_less_equal(MOT_Fluorescence_atoms1/100., 0.)
        
        
        dark_mask = np.ma.masked_less_equal(MOT_Fluorescence_dark/100., 0.)
        
        atoms2_mask = np.ma.masked_less_equal(MOT_Fluorescence_atoms2/100., 0.)
        
        atoms1 = np.matrix(atoms1_mask,dtype = np.float16)
        dark = np.matrix(dark_mask,dtype = np.float16)
        atoms2 = np.matrix(atoms2_mask,dtype = np.float16)
                
                
        return atoms1, dark, atoms2
        
                        
# Main
try:
    #plt.xkcd()
    with h5py.File(path) as h5_file:
        if '/data' in h5_file:
            print_time('Calculating count...')
        # Get OD
            
            atoms1_, dark_, atoms2_  = raw_to_OD(path)
            print_time('Get count...')
            atoms1 = ODShot(atoms1_)
            dark = ODShot(dark_)
            atoms2 = ODShot(atoms2_)
            #ODrot = OD.rotate_OD(angle=4.5)
            F_red, mF_red, _ROI_atoms1, BCK_atoms1 =  atoms1.get_ROI(sniff=False, get_background=False) 
            
            ROI_atoms1 = ODShot(_ROI_atoms1)
            
            BCK_atoms1 = np.mean(BCK_atoms1)*np.ones(_ROI_atoms1.shape)
           
            run.save_result( 'pk_count_atoms1', (np.amax(atoms1_.astype(float16))))
            
			#[100:280,60:300] magnetic trap position
            #[120:300,130:500] Mot position
            
        # Compute number
            if True: #stored == "z-TOF":
                count_atoms1 = np.sum(atoms1_[100:280,60:300])
                
                run.save_result('count_atoms1' , count_atoms1)
               
            else:
                count = 0
                
            
            
        # Display figure with shot, slices and fits and N
            fig = figure(1, figsize=(16, 16), frameon=False)
            gs = gridspec.GridSpec(2, 2, width_ratios=[1,1], height_ratios=[1,1])
            subplot(gs[3])
            
            
            gca().axison = False
            tight_layout()
        # OD and ROI display
            subplot(gs[0])
            im0= imshow(atoms1_[100:280,60:300]/np.max(atoms1_), vmin= -0.0, vmax = 1.0, cmap='RdYlBu_r', aspect='auto', interpolation='none')
            divider = make_axes_locatable(gca())
            cax = divider.append_axes("right", "5%", pad="3%")
            colorbar(im0, cax=cax)
            title('atoms1')
            tight_layout()
            
            subplot(gs[1])
            im0= imshow(dark_/np.max(dark_), vmin= -0.0, vmax = 1.0, cmap='RdYlBu_r', aspect='auto', interpolation='none')
            divider = make_axes_locatable(gca())
            cax = divider.append_axes("right", "5%", pad="3%")
            colorbar(im0, cax=cax)
            title('dark')
            tight_layout()
            
            subplot(gs[2])
            im0= imshow(atoms2_/np.max(atoms2_), vmin= -0.0, vmax = 1.0, cmap='RdYlBu_r', aspect='auto', interpolation='none')
            divider = make_axes_locatable(gca())
            cax = divider.append_axes("right", "5%", pad="3%")
            colorbar(im0, cax=cax)
            title('atoms2')
            tight_layout()
            
            
		
			
			
except Exception as e:
    print '%s' %e +  os.path.basename(path)
    print '\n ********** Not Successful **********\n\n'        