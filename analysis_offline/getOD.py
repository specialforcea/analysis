from __future__ import division
from os import listdir
from os.path import isfile, join
from scipy import optimize, signal
from common.traces import gaussian, fit_gaussian, lorentzian, fit_lorentzian, bimodal, fit_bimodal
#from yuchen_analysis.common.OD_handler_fl import ODShot
from common.fluctuation_methods import fluctuations_from_mean, covariance
import h5py
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import image
import matplotlib.gridspec as gridspec
from pylab import *
from time import time
from scipy.ndimage import *
from common.OD_handler import ODslice
import fit_table
from common.traces import bimodal1D

# Parameters
pixel_size = 5.6e-6/1.19  # Divided by Magnification Factor
sigma0 = 3*(780.24e-9)**2/(2*np.pi)  # Atomic cross section (resonant absorption imaging)

def get_all_files(sequence_path, ongoing=True):
    root = r'L:\internal\Spielman-group\RbChip\Data\Experiments\rb_chip'
    filepath = root + sequence_path
    all_files = [f for f in listdir(filepath) if isfile(join(filepath, f)) if f.endswith(".h5")]
    if ongoing:
        del all_files[-1]
    return filepath, all_files
 
# List of shots
path_1, shots_list_1 = get_all_files(r'\2018_02_09_RecoverBEC_MoveUp_Evap_Decompress1_2_DipoleEvapExpo\2018\02\09\0038', False)

def raw_to_OD(filep, fpath, global_name=None):
    with h5py.File(join(filep, fpath)) as h5_file:
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
        return calculated_OD_TOF, atoms_TOF, probe_TOF, bckg_TOF, h5_file['globals'].attrs[global_name] 
        
        
def get_data_from_shots(path, shot_list, return_global=False):
    OD, atoms, probe, bckg, glob = [], [], [], [], []
    for shot in range(len(shot_list)):
        OD.append(raw_to_OD(path, shot_list[shot], 'short_TOF')[0])
        atoms.append(raw_to_OD(path, shot_list[shot], 'short_TOF')[1])
        probe.append(raw_to_OD(path, shot_list[shot], 'short_TOF')[2])
        bckg.append(raw_to_OD(path, shot_list[shot], 'short_TOF')[3])
        glob.append(raw_to_OD(path, shot_list[shot], 'short_TOF')[4])
    if return_global:
        return np.array(OD), np.array(atoms), np.array(probe), np.array(bckg), np.array(glob)
    else: 
        return np.array(OD), np.array(atoms), np.array(probe), np.array(bckg)
#operations:
plot_OD = True
plot_N = False
vs_variable = " "
get_results = False


OD, atoms, probe, bckg = get_data_from_shots(path_1, shots_list_1, return_global=False)
if get_results:
    l = len(OD)
    N = np.zeros(l, dtype=int)
    max_OD = np.zeros(l, dtype=int)
    fit_x = np.zeros((2,l))#0 is center, 1 is width
    fit_y = np.zeros((2,l))
    for sh in range(l):
        
        N[sh] = np.sum((OD[sh].T[::-1,][300:540,100:400]/sigma0)*pixel_size**2)
        max_OD[sh] = np.amax(OD[sh].T[::-1,][100:350,150:400],axis=None)
            
            
        ycolOD, y_ax = np.mean(np.array(OD[sh][:, 260:280]), axis=1), np.linspace(0, 648, 648)
       
        _yslice_ = ODslice(slice_OD=ycolOD, slice_axis=y_ax)
        try:
            y_gauss_pars, y_dense_gauss, y_gaussian_fit = _yslice_.fit_gauss()
            #y_tf_pars, x_dense_tf, y_TF_fit = _yslice_.fit_pure_thomas_fermi()
            #y_bimodal_pars, y_dense_bimodal, y_bimodal_fit =_yslice_.fit_bimodal()
            fit_y_success = True
            print 'y slice fit'
            fit_table.get_params(y_gauss_pars)
            fit_y[0,sh] = np.abs(y_gauss_pars[1]*pixel_size/1e-6)
            fit_y[1,sh] = np.abs(y_gauss_pars[2]*pixel_size/1e-6)
            
                
        except Exception as e:
            fit_y_success = False
            print 'Fit of y slice unsuccessful, %s' %e
            
        xcolOD, x_ax = np.mean(np.array(OD[sh][int(y_gauss_pars[1])-5:int(y_gauss_pars[1])+5, :]), axis=0), np.linspace(0,488, 488)
        _xslice_ = ODslice(slice_OD=xcolOD, slice_axis=x_ax)
        try:
            x_gauss_pars, x_dense_gauss, x_gaussian_fit =_xslice_.fit_gauss()
            #x_tf_pars, x_dense_tf, x_TF_fit = _xslice_.fit_pure_thomas_fermi()
            #x_bimodal_pars, x_dense_bimodal, x_bimodal_fit =_xslice_.fit_bimodal()
            fit_x_success = True
            print 'x slice fit'
            fit_table.get_params(x_gauss_pars)
            fit_x[0,sh] = np.abs(x_gauss_pars[1]*pixel_size/1e-6)
            fit_x[1,sh] = np.abs(x_gauss_pars[2]*pixel_size/1e-6)
            
            
        except Exception as e:
            fit_x_success = False
            print 'Fit of x slice unsuccessful, %s' %e
        
        
        

if plot_OD:
        
    
    
    
    
    fig  = figure(1, figsize=(16,16), frameon=False)
    gs = gridspec.GridSpec(2, 2, width_ratios=[1,1], height_ratios=[1,1])
    subplot(gs[0])
    im0= plt.imshow(OD[0].T[::-1,], vmin= -0.0, vmax =2.8,  cmap='RdYlBu_r', aspect=0.75, interpolation='none')
    colorbar(im0)
    title('OD_TOF atom#= + str(N)')
    
    subplot(gs[1])
    im1= imshow(atoms[0].T[::-1,], vmin= -0.0, cmap='RdYlBu_r', aspect=0.75, interpolation='none')
    colorbar(im1)
    title('atoms_TOF')
    
    subplot(gs[2])
    im2= imshow(probe[0].T[::-1,], vmin= -0.0, cmap='RdYlBu_r', aspect=0.75, interpolation='none')
    colorbar(im2)
    title('probe_TOF')
    
    subplot(gs[3])
    im3= imshow(bckg[0].T[::-1,], vmin= -0.0, cmap='RdYlBu_r', aspect=0.75, interpolation='none')
    colorbar(im3)
    title('bckg_TOF')

if plot_N:
    



    
