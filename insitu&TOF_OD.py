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
from common.wave_process import single_wave, multi_wave
from common.traces import bimodal1D


# Parameters
pixel_size = 5.6e-6/1.19# Divided by Magnification Factor as measured 12/20/2017
            # 5.6e-6/5.33 for z in situ                        # Yuchen and Paco: 08/19/2016
            #5.6e-6/3.44 for z TOF                           # Yuchen and Paco: 08/19/2016
            #5.6e-6/2.72 for x-in situ                        # Paco: 05/06/2016
sigma0 = 3*(780e-9)**2/(2*3.14)
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
            image_group['imagesYZ_2_Flea3'] = 'MOT_abs'
            Isat = 100          
            alpha = 1 #1.645
        if  '/data/imagesXY_1_Flea3' in h5_file:
            image_group['imagesXY_1_Flea3'] = 'z-TOF'
            Isat_TOF = 3000          # In counts // Emine:01/17/2018
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
                # if id == 'z-insitu':
                    # probe_insitu = (np.array(h5_file['data'][appended]['Raw'])[0])
                    # atoms_red = (np.array(h5_file['data'][appended]['Raw'])[1])
                    # atoms_blue = (np.array(h5_file['data'][appended]['Raw'])[2])
                    # bckg_insitu =  (np.array(h5_file['data'][appended]['Raw'])[3])
                    
                if id == 'z-TOF':
                    atoms_TOF = (np.array(h5_file['data'][appended]['Raw'])[0])
                    probe_TOF = (np.array(h5_file['data'][appended]['Raw'])[1])
                    bckg_TOF =  (np.array(h5_file['data'][appended]['Raw'])[2])
                if id == 'x-insitu':
                    atoms_xinsitu = (np.array(h5_file['data'][appended]['Raw'])[0])
                    probe_xinsitu = (np.array(h5_file['data'][appended]['Raw'])[1])
                    bckg_xinsitu =  (np.array(h5_file['data'][appended]['Raw'])[2])
                if id == 'MOT_abs':
                    atoms_MOT = (np.array(h5_file['data'][appended]['Raw'])[0])
                    probe_MOT = (np.array(h5_file['data'][appended]['Raw'])[1])
                    bckg_MOT =  (np.array(h5_file['data'][appended]['Raw'])[2])
                    
        # div_red = np.ma.masked_invalid((atoms_red - bckg_insitu)/(probe_insitu-bckg_insitu))
        # div_red = np.ma.masked_less_equal(div_red, 0.)
        # another_term_red = (probe_insitu-atoms_red)/(Isat_insitu)
        
        # div_blue = np.ma.masked_invalid((atoms_blue - bckg_insitu)/(probe_insitu-bckg_insitu))
        # div_blue = np.ma.masked_less_equal(div_blue, 0.)
        # another_term_blue = (probe_insitu-atoms_blue)/(Isat_insitu)
        
        numerator = atoms_TOF - bckg_TOF
        numerator[numerator <= 0.] = 1.
        denominator = probe_TOF-bckg_TOF
        denominator[denominator <= 0.] = 1.
        div_TOF = numerator/denominator
        #div_TOF = np.ma.masked_invalid(div_TOF)
        #div_TOF = np.ma.masked_less_equal(div_TOF, 1e-8)
        
        #div_TOF = np.ma.masked_invalid((atoms_TOF - bckg_TOF)/(probe_TOF-bckg_TOF))
        #div_TOF = np.ma.masked_less_equal(div_TOF, 0.)
        another_term_TOF = (probe_TOF-atoms_TOF)/(Isat_TOF)
        
        #div_TOF = np.ma.masked_invalid((atoms_TOF - bckg_TOF)/100)
        
        # div_MOT = np.ma.masked_invalid((atoms_MOT - bckg_MOT)/(probe_MOT-bckg_MOT))
        # div_MOT = np.ma.masked_less_equal(div_MOT, 0.)
        # another_term_MOT = (probe_MOT-atoms_MOT)/(Isat)
        
        # div_xinsitu = np.ma.masked_invalid((atoms_xinsitu - bckg_xinsitu)/(probe_xinsitu-bckg_xinsitu))
        # div_xinsitu = np.ma.masked_less_equal(div_xinsitu, 0.)
        # another_term_xinsitu = (probe_xinsitu-atoms_xinsitu)/(Isat)
        
        # calculated_OD_red = np.matrix(-alpha*np.log(div_red)+0*another_term_red)
        # calculated_OD_blue = np.matrix(-alpha*np.log(div_blue)+0*another_term_blue)
        #calculated_OD_MOT = np.matrix(-alpha*np.log(div_MOT)+another_term_MOT)
        calculated_OD_TOF = np.matrix(-alpha*np.log(np.ma.masked_invalid(div_TOF))+another_term_TOF)
                
                
       # return calculated_OD_red, calculated_OD_blue, calculated_OD_TOF
        return calculated_OD_TOF, atoms_TOF, probe_TOF, bckg_TOF
        
                        
# Main
try:
    #plt.xkcd()
    with h5py.File(path) as h5_file:
        if '/data' in h5_file:
            print_time('Calculating OD...')
        # Get OD
            
            _OD_TOF, atoms_TOF, probe_TOF, bckg_TOF  = raw_to_OD(path)
            print_time('Get OD...')
            # OD_red = ODShot(_OD_red)
            # OD_blue = ODShot(_OD_blue)
            OD_TOF= ODShot(_OD_TOF)
            
            #ODrot = OD.rotate_OD(angle=4.5)
            # F_red, mF_red, _ROI_red, BCK_a_red =  OD_red.get_ROI(sniff=False, get_background=False) 
            # F_blue, mF_blue, _ROI_blue, BCK_a_blue =  OD_blue.get_ROI(sniff=False, get_background=False) 
            # ROI_red = ODShot(_ROI_red)
            # ROI_blue = ODShot(_ROI_blue)
            # BCK_red = np.mean(BCK_a_red)*np.ones(_ROI_red.shape)
            # BCK_blue = np.mean(BCK_a_blue)*np.ones(_ROI_blue.shape)
            # run.save_result( 'pkOD_red', (np.amax(_ROI_red.astype(float16))))
            # run.save_result( 'pkOD_blue', (np.amax(_ROI_blue.astype(float16))))
            
			
            N = np.sum((_OD_TOF.T[::-1,][20:550,100:400]/sigma0)*pixel_size**2)
            run.save_result('N', N)
            max_OD = np.amax(_OD_TOF.T[::-1,][100:350,150:400],axis=None)
            run.save_result('max_OD', max_OD)
            #ROI = [0:400,300:650]
            
        # Compute number
            # if False: #stored == "z-TOF":
                # N_red = (np.sum((_ROI_red-BCK_red)/sigma0)*pixel_size**2)
                # N_blue = (np.sum((_ROI_blue-BCK_blue)/sigma0)*pixel_size**2)
                # run.save_result(('N_red(' + str(F_red) +',' +str(mF_red)+')'), N_red)
                # run.save_result(('N_blue(' + str(F_blue) +',' +str(mF_blue)+')'), N_blue)
            # else:
                # N = 0
                
            #imb =np.float( (N_red-N_blue)/(N_blue+N_red))
            
            #flip image due to camera orientation
            
            _OD_fit = _OD_TOF[::-1] # Reverse to match gravity
            
            
            
           
        
        # Display figure with shot, slices and fits and N
            fig  = figure(1, figsize=(16,16), frameon=False)
            gs = gridspec.GridSpec(2, 2, width_ratios=[1,1], height_ratios=[1,1])
            subplot(gs[0])
            im0= imshow(_OD_TOF.T[::-1,], vmin= -0.0, vmax =2.8,  cmap='RdYlBu_r', aspect=0.75, interpolation='none')
            colorbar(im0)
            title('OD_TOF atom#=' + str(N))
            
            subplot(gs[1])
            im1= imshow(atoms_TOF.T[::-1,], vmin= -0.0, cmap='RdYlBu_r', aspect=0.75, interpolation='none')
            colorbar(im1)
            title('atoms_TOF')
            
            subplot(gs[2])
            im2= imshow(probe_TOF.T[::-1,], vmin= -0.0, cmap='RdYlBu_r', aspect=0.75, interpolation='none')
            colorbar(im2)
            title('probe_TOF')
            
            subplot(gs[3])
            im3= imshow(bckg_TOF.T[::-1,], vmin= -0.0, cmap='RdYlBu_r', aspect=0.75, interpolation='none')
            colorbar(im3)
            title('bckg_TOF')
            
            
            
            # fig = figure(1, figsize=(16, 16), frameon=False)
            # gs = gridspec.GridSpec(2, 2, width_ratios=[1,1], height_ratios=[1,1])
            # subplot(gs[3])
            # str_imb = r'imbalance = %s' % imb		
            # str_red = r'atoms_red = %s' % N_red
            # str_blue = r'atoms_blue = %s' % N_blue
            # text(0.2, 0.5, str_imb, ha='left', va='bottom',fontsize=18)
            # text(0.2, 0.7 , str_red, ha='left', va='top',fontsize=18)
            # text(0.2, 0.6, str_blue, ha='left', va='center',fontsize=18)
            # gca().axison = False
            # tight_layout()
        # # OD and ROI display
            # subplot(gs[0])
            # im0= imshow(_OD_red, vmin= -0.0, vmax = 0.8, cmap='RdYlBu_r', aspect='auto', interpolation='none')
            # divider = make_axes_locatable(gca())
            # cax = divider.append_axes("right", "5%", pad="3%")
            # colorbar(im0, cax=cax)
            # title('OD_red')
            # tight_layout()
            
            # subplot(gs[1])
            # im0= imshow(_OD_blue, vmin= -0.0, vmax = 0.8, cmap='RdYlBu_r', aspect='auto', interpolation='none')
            # divider = make_axes_locatable(gca())
            # cax = divider.append_axes("right", "5%", pad="3%")
            # colorbar(im0, cax=cax)
            # title('OD_blue')
            # tight_layout()
            
            # subplot(gs[2])
            # im0= imshow(_OD_TOF, vmin= -0.0, vmax = 0.8, cmap='RdYlBu_r', aspect='auto', interpolation='none')
            # divider = make_axes_locatable(gca())
            # cax = divider.append_axes("right", "5%", pad="3%")
            # colorbar(im0, cax=cax)
            # title('OD_TOF')
            # tight_layout()
            
            
		
			
			
except Exception as e:
    print '%s' %e +  os.path.basename(path)
    print '\n ********** Not Successful **********\n\n'        