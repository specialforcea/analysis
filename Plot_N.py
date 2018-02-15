
from __future__ import division
from lyse import *
from time import time
from matplotlib import pyplot as plt
from common.traces import *
from common.plot_pretty import plot_pretty
from common.wave_process import single_wave, multi_wave
import os
import pandas as pd
import numpy as np

print '\nRunning %s' % os.path.basename(__file__)
t = time()

@plot_pretty(marker='o', 
             markeredgecolor='royalblue', 
             markerfacecolor='lightsteelblue', 
             markeredgewidth=1.5, 
             markersize=6, 
             linestyle='-', 
             color='midnightblue', 
             linewidth=1.5)
def plot_waves(wave_x, wave_y, show_errorbar=False, is_data_fit=False):
    if show_errorbar:
        plotline, caplines, barlinecols = plt.errorbar(wave_x['wave_data'], 
                                                       wave_y['wave_data'],
                                                       xerr = 0*wave_x['wave_uncer'], 
                                                       yerr = wave_y['wave_uncer'], 
                                                       capsize=0, 
                                                       ls='None')
        line_objects = plotline, caplines, barlinecols
    elif is_data_fit:
        line_objects = plt.plot(wave_x['wave_data'], 
                                wave_y['wave_data'], 
                                label= wave_y['wave_name'] + ' ' + str(wave_y['wave_fit_params']))
    else:
        line_objects = plt.plot(wave_x['wave_data'], 
                                wave_y['wave_data'], 
                                label= wave_y['wave_name'])    
    plt.xlabel(wave_x['wave_name'], fontsize=12)
    plt.ylabel(wave_y['wave_name'], fontsize=12)
    plt.title('Plot: ' + wave_y['wave_name'] + ' vs ' + wave_x['wave_name'], fontsize=12)
    return line_objects, show_errorbar, is_data_fit
    
try:
    # Call waves from main dataframe
    if True:
        # _RunNo_ = single_wave(wave_name='RunNo').get_single()
        # number_1 = single_wave(wave_name='N_(1,-1)', single_shot_routine='get_OD').normalize(norm='val', val=1.0)
        # run_no = multi_wave(_RunNo_, number_1).build_RunNo_wave()
        # time = single_wave(wave_name='short_TOF').normalize(norm='val', val=1.0)
        # pkOD =  single_wave(wave_name='pkOD',  single_shot_routine='get_OD').get_single()
        # temp = single_wave(wave_name='temperature', single_shot_routine='get_OD').normalize(norm='val', val=1e-9)
        # width = single_wave(wave_name='x_gauss_width', single_shot_routine='get_OD').get_single()
        # com = single_wave(wave_name='x_gauss_center', single_shot_routine='get_OD').get_single()
        #counts = single_wave(wave_name='count_atoms1', single_shot_routine='MOT fluor count').get_single()
        #magtime = single_wave(wave_name='MagTrapTimeAfterMoveUp').get_single()
        run_number = single_wave(wave_name='run number').get_single()
        f_max = single_wave(wave_name='f_max').get_single()
        EndRF = single_wave(wave_name='EndRF').get_single()
        StartRF = single_wave(wave_name='StartRF').get_single()
        imfreq = single_wave(wave_name='ImageFreq').get_single()
        rc2 = single_wave(wave_name='ReduceCurrentF2').get_single()
        rc = single_wave(wave_name='ReduceCurrentF').get_single()
        rc1 = single_wave(wave_name='ReduceCurrentF1').get_single()
        TpyShimTrap = single_wave(wave_name='TpyShimTrap').get_single()
        TpxShim = single_wave(wave_name='TpxShim').get_single()
        TpxShimTrap = single_wave(wave_name='TpxShimTrap').get_single()
        TpzShim = single_wave(wave_name='TpzShim').get_single()
        dint = single_wave(wave_name='DipoleInt').get_single()
        dint2 = single_wave(wave_name='DipoleInt2').get_single()
        dint3 = single_wave(wave_name='DipoleInt3').get_single()
        fd = single_wave(wave_name='FinalDipole').get_single()
        evaptime = single_wave(wave_name='RFevapTime').get_single()
        evaptime1 = single_wave(wave_name='EvapTime').get_single()
        evaptime2 = single_wave(wave_name='EvapTime_2').get_single()
        evaptime3 = single_wave(wave_name='EvapTime_3').get_single()
        evaptimeTotal = single_wave(wave_name='EvapTimeTotal').get_single()
        dipoleEvapTau = single_wave(wave_name='dipoleEvapTau').get_single()
        DecompTime2 = single_wave(wave_name='DecompTime2').get_single()
        DecompTime1 = single_wave(wave_name='DecompTime1').get_single()
        fdint = single_wave(wave_name='FinalDipole').get_single()
        TpzShimEnd = single_wave(wave_name='TpzShimEnd').get_single()
        TpyShimTrap1 = single_wave(wave_name='TpyShimTrap1').get_single()
        TpyShimTrap2 = single_wave(wave_name='TpyShimTrap2').get_single()
        TpyShimTrap3 = single_wave(wave_name='TpyShimTrap3').get_single()
        delay_time = single_wave(wave_name='set_freq_delay_time').get_single()
        imagefreq = single_wave(wave_name='ImageFreq').get_single()
        UVxShim = single_wave(wave_name='UVxShim').get_single()
        UVyShim = single_wave(wave_name='UVyShim').get_single()
        UVzShim = single_wave(wave_name='UVzShim').get_single()
        MOTxShim = single_wave(wave_name='MOTxShim').get_single()
        MOTyShim = single_wave(wave_name='MOTyShim').get_single()
        MOTzShim = single_wave(wave_name='MOTzShim').get_single()
        MolXBias = single_wave(wave_name='MolXBias').get_single()
        MolYBias = single_wave(wave_name='MolYBias').get_single()
        MolZBias = single_wave(wave_name='MolZBias').get_single()
        OpxShim = single_wave(wave_name='OpxShim').get_single()
        OpyShim = single_wave(wave_name='OpyShim').get_single()
        OpzShim = single_wave(wave_name='OpzShim').get_single()
        CMOTFreq = single_wave(wave_name='CMOTFreq').get_single()
        RepumpCMOT  = single_wave(wave_name='RepumpCMOT').get_single()
        MOTAmplitude = single_wave(wave_name='MOTAmplitude').get_single()
        UVMOTRepump = single_wave(wave_name='UVMOTRepump').get_single()
        UVMOTAmplitude = single_wave(wave_name='UVMOTAmplitude').get_single()
        RepumpMOT  = single_wave(wave_name='RepumpMOT').get_single()
        UVMOTFreq   = single_wave(wave_name='UVMOTFreq').get_single()
        MOTFreq   = single_wave(wave_name='MOTFreq').get_single()
        CMOTFreq   = single_wave(wave_name='CMOTFreq').get_single()
        MOT_time = single_wave(wave_name='MOT_time').get_single()
        StartFreqMol = single_wave(wave_name='StartFreqMol').get_single()
        EndFreqMol = single_wave(wave_name='EndFreqMol').get_single()
        TimeMol = single_wave(wave_name='TimeMol').get_single()
        TauMol = single_wave(wave_name='TauMol').get_single()
        OpPumpAmplitude = single_wave(wave_name='OpPumpAmplitude').get_single()
        RepumpOpPump = single_wave(wave_name='RepumpOpPump').get_single()
        TimePump = single_wave(wave_name='TimePump').get_single()
        CapxShim = single_wave(wave_name='CapxShim').get_single()
        CapyShim = single_wave(wave_name='CapyShim').get_single()
        CapzShim = single_wave(wave_name='CapzShim').get_single()
        CaptureWidth  = single_wave(wave_name='CaptureWidth').get_single()
        EndFreqOpPump   = single_wave(wave_name='EndFreqOpPump').get_single()
        TrapTime    = single_wave(wave_name='TrapTime').get_single()
        MOTCaptureCurrent   = single_wave(wave_name='MOTCaptureCurrent').get_single()
        TimePrePump    = single_wave(wave_name='TimePrePump').get_single()
        UVMOTCurrent     = single_wave(wave_name='UVMOTCurrent').get_single()
        Vel9  = single_wave(wave_name='Vel9').get_single()
        #x_amp = single_wave(wave_name='x_gauss_amp', single_shot_routine='fit_slice').get_single()
        
        tof = single_wave(wave_name='TimeShTOF').get_single()
        #tof_squared = tof
        #tof_squared['wave_data'] = tof_squared['wave_data']*tof_squared['wave_data']
        #tof_squared['wave_name'] = 'TimeShTOF squared'
        probeInt = single_wave(wave_name='ProbeInt').get_single()
        #motlocktime = single_wave(wave_name='MOTLockFreqTimeOffset').get_single()
        
        holdtime = single_wave(wave_name='MagTrapTimeAfterMoveUp').get_single()
        #reduce_cur = single_wave(wave_name='ReduceCurrentF').get_single()
        N = single_wave(wave_name='N', single_shot_routine='fit_slice').get_single()
        max_OD = single_wave(wave_name='max_OD', single_shot_routine='insitu&TOF_OD').get_single()
        x_width = single_wave(wave_name='x_gauss_width', single_shot_routine='fit_slice').get_single()
        #x_width_squared = x_width
        #x_width_squared['wave_data'] = x_width_squared['wave_data']*x_width_squared['wave_data']
        #x_width_squared['wave_name'] = 'x_gauss_width squared'
        
        
        y_width = single_wave(wave_name='y_gauss_width', single_shot_routine='fit_slice').get_single()
        #y_width_squared = y_width
        #y_width_squared['wave_data'] = y_width_squared['wave_data']*y_width_squared['wave_data']
        #y_width_squared['wave_name'] = 'y_gauss_width squared'
        y_center = single_wave(wave_name='y_gauss_center', single_shot_routine='fit_slice').get_single()
        
    # 1D Fourier analysis
    if False:
        number_psd = single_wave(wave_name='N_(2,2)', single_shot_routine='get_OD').psd()
        f_space = single_wave(wave_name='TimeHold').reciprocal_space()
        
    # Spin projection analysis
    if False:
        pulse_time = single_wave(wave_name='time_microwave_pulse').get_single()
        #number_m1 = single_wave(wave_name='N_(1,-1)', single_shot_routine='get_OD').get_single()
        #number_0 = single_wave(wave_name='N_(1,0)', single_shot_routine='get_OD').get_single()
        number_2 = single_wave(wave_name='N_(2,2)', single_shot_routine='get_OD').get_single()
        
        oneplustwo = multi_wave(number_2, number_1).add_waves()
        #total_number = multi_wave(onepluszero, number_m1).add_waves()
        #pop_ratio_0 = multi_wave(total_number, number_0).take_ratio()
        pop_ratio_2 = multi_wave(oneplustwo, number_2).take_ratio()
        pop_ratio_1 = multi_wave(oneplustwo, number_1).take_ratio()
    
    # Fits 
    if False:
        try:
            fit_x1, fit_y1= multi_wave(tof ,x_width).fit_wave(trace='parabola')
        except:
            fit_x1, fit_y1 = magtime, counts
    # Plot multi-plots on demand
    
    plot_triplet = True
    plotShotNumber = False
    scanned_variable = rc1
    if plot_triplet:
        fig1 = plt.figure(1, figsize=(40, 40))
           
        
        fig1.add_subplot(131) 
        plot_waves(scanned_variable, N, show_errorbar=False, is_data_fit=False)
        fig1.add_subplot(132) 
        plot_waves(scanned_variable, x_width, show_errorbar=False, is_data_fit=False)
        fig1.add_subplot(133) 
        psd = N['wave_data']/(x_width['wave_data']*y_width['wave_data']/1e12)**3
        # # psd = psd/psd[0]
        plt.scatter(scanned_variable['wave_data'],psd)
        plt.title('relative psd vs ' + scanned_variable['wave_name'])
        plt.xlabel(scanned_variable['wave_name'])
        plt.ylabel('relative psd')
        
        
        
        #plt.legend(['two stages decomp','one stage decomp'])
        #plt.xlabel('run number')
        #plt.ylabel('atom number')
        #plt.scatter(NN[5:],max_OD['wave_data'][5:],'b')
        #t = linspace(0.0,0.014,100)
        #yy = parabola(t,fit_y1['wave_fit_params'][0],fit_y1['wave_fit_params'][1])
        #plt.plot(t,yy)
        #plot_waves(fit_x1, fit_y1, show_errorbar=False, is_data_fit=True)
        #plot_waves(fit_freq, fit_number_1, show_errorbar=False, is_data_fit=True)  
    if plotShotNumber:
        fig2 = plt.figure(2)
        NN = range(N['wave_data'].shape[0])
        plt.plot(NN[0:5],N['wave_data'][0:5],'o',NN[5:],N['wave_data'][5:],'v')
        plt.show()
        plt.tight_layout()
    
    # Show
    fig1 = plt.gcf()
    fig2 = plt.gcf()
    
    
    print '\n ********** Successful **********\n\n'
except Exception as e:
    print '%s' %e 
    print '\n ********** Not Successful **********\n\n'

    
""" Wishlist """

    #avg_wave_x, avg_wave_y = average_wave(('None', 'get_OD', 'None'), ('sequence_index', 'y_gauss_center', 'short_TOF'),
    #                                      as_wave_x='short_TOF', wave_avg='sequence_index', common='y_gauss_center')
    
# - Parse fit_params guess all the way into traces.py
