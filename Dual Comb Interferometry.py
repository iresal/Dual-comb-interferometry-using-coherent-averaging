"""
- Dual-Comb Spectroscopy - Post-processing
    Script dedicated to post-process the data measured in the DCI setup.
    The system requires a calibration step using a first measurement as reference, this data is processed
    and the complex amplitude is retrieved. The subsecuent measurements are called SIGNAL, the same 
    processing is performed but with the aditional step of substractinf the reference complex amplitude. 
    Import MAT files for both cases
    
    1. REFERENCE. Output of the balanced photodetector (BPD) and 25MHz signal without a DUT.
    2. SIGNAL. BPD output and 25MHz signal with the DUT (WS, fiber, etc).  

Created on 2020-06-26
@autor: Israel Rebolledo (israels@chalmers.se)

"""
import os
import numpy as np
import scipy.io as sio
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.signal import find_peaks,peak_widths
from scipy.fftpack import fft, ifft,ifftshift,fftshift
import time
start_time = time.time()

os.chdir('C:\\Users\\israels\\Box\\Ultrafast Photonics\\Isra\\Measurements\\DCI\\') # Office path 
#os.chdir('X:\\Measurements\\DCI\\0831\\') # Lab path 

"""
- Functions Description 
    - Clock Extraction. Extracts the zero crossing index of a signal for the ascending edge, we use
    the second derivative of the signal to find the zero-crossing.
"""
def clock_extr (clock_ref, delta_t):
    clock = clock_ref.copy()   #just copy to another array 
    f = np.fft.fftfreq(N,delta_t) # Frequency array for the full length
    fft_clock = fft(clock)   
    filter_clock = (abs(f)< 27.5e6) & (abs(f)> 22.5e6) #Just a BP filter  
    filtered = ifft(fft_clock*filter_clock)            
    diff_clock = np.diff(np.sign(np.diff(np.real(filtered))))  #Second derivative 
    out = np.where(diff_clock >0)

    return out

def load_data (filename):
    """
    Load the data 
    We will be using MAT files to record data, the purpose is to have data files that can be used
    in Matlab and Python. """
    data_reference = sio.loadmat(filename)
    clock_ref = data_reference['clock']
    signal_ref = data_reference['signal']
    clock_ref = np.array(clock_ref).ravel()  #Convert to 1D np array
    signal_ref = np.array(signal_ref).ravel()

    return [clock_ref,signal_ref]

def interferogram_slicing(clock_ref,signal_ref):
    ### Spectral Processing of the data
    # Here we extract individual interferograms of 40ns duration using the 25MHz as reference
    ref_index = np.array(clock_extr(clock_ref,delta_t)).ravel() # Make the index extraction 
    interferograms = np.zeros(((n_interferograms-start_interf),np.int(samples))) #Create an empy array with the actual number of samples

    for i in range (0,len(interferograms)): #Indexing the interferograms 
        interferograms[i,:] = signal_ref[ref_index[i+start_interf]:(ref_index[i+start_interf]+np.int(samples))]

    return interferograms

def find_ideal_peaks ():
    #Ideal peaks for RF Comb
    ideal_peaks = f_AOM + delta_f/2 +np.arange(-n_lines/2,n_lines/2) *delta_f 
    pks_idx = np.zeros(len(ideal_peaks))
    
    for i in range (0,len(ideal_peaks)):
        pks_idx[i] = np.where (freq ==ideal_peaks[i])[0]
    pks_idx = pks_idx.astype(int)   # Convert the index into INT array for indexing

    return pks_idx

def spectral_processing (fft_interf):
    ## Spectral Amplitude
    amp_array = np.abs(fft_interf)       # Amplitude Spectrum 
    amp_array = amp_array[:,pks_idx]     # Mapping to select only where the RF notes are
    ## Spectral Phase
    angle = np.angle(fft_interf)
    angle = angle[:,pks_idx]
    phase_array = np.unwrap(angle)  # Contain the phase from peaks recognized for the 'n' interferograms processed
    # Removing the linear term of the Phase profile 
    linear_term = np.zeros((n_interferograms,n_lines))
    for jj in range(0,n_interferograms):
        p = np.polyfit(freq_peaks,phase_array[jj,:],1)
        fit_curve = np.poly1d(p)
        linear_term[jj,:] = fit_curve(freq_peaks) 
    
    return amp_array,phase_array,linear_term

if __name__ == "__main__": 
    
    #####################   Parameters of the system    #####################
    fr =  50e9          #Sampling frequency
    delta_t = 1/fr
    delta_f = 100e6     #Repetition rate of the RF Comb
    f_AOM = 25e6        #Frequency of Acousto Optic Modulator
    n_lines = 45        #Number of spectral lines
    N = int(5e6)        #Number of samples recorded
    n_interferograms = 1000   #Number of interferograms to process
    start_interf = 0         #Start in the second interferogram
    samples = 4/(delta_f*delta_t)   # Number of samples for individual interferograms 
    t = np.arange(0,N*delta_t,delta_t)  #time array 
    freq = np.fft.fftfreq(np.int(samples),delta_t) # Frequency vector
    print('Number of comb lines = ', n_lines)
    #####################################################################
    #####################   Reference Measurement   #####################
    filename = '0504_ZeroPhase_50Gs_100us_compensated.mat' #This is our reference measurement file
    #filename = '1020_50Gs_100us_1545_25_ref.mat'
    [clock_ref, signal_ref] = load_data(filename)
    ###############
    fft_signal = ifftshift(fft(fftshift(signal_ref)))
    f_grid = np.linspace(-fr/2,fr/2,fft_signal.size, endpoint=False)
    #plt.figure()
    #plt.semilogy(f_grid,abs(fft_signal))
    #plt.show()
    filename=('Y:\\Lab Members Data\\Isra\\data_all.mat')
    sio.savemat(filename, {'f_grid': f_grid,'fft_signal':fft_signal})
    #plt.plot(f_grid,abs(fft_signal))
    


    ####################

    ### First we chop the recorded signal to extract individual interferograms
    interferograms = interferogram_slicing(clock_ref,signal_ref)
    pks_idx = find_ideal_peaks() # Find the vector position of the beating frequencies according the ideal peaks
    freq_peaks = freq[pks_idx]   # Creates a new RF vector according the peaks recognized
    ### Fourier spectral processing
    fft_interf = fft(interferograms)   
    [amp_array,phase_array,linear_term] = spectral_processing(fft_interf)
    phase_RF = phase_array - linear_term
    mean_amp = np.mean(amp_array, axis = 0)
    mean_phase = np.mean(phase_RF, axis = 0)  # Recovered phase from reference measurement 
    
    ##################################################################
    #####################   Signal Measurement   #####################
    filename_dtu = '0504_2Disp_50Gs_100us.mat'  #This is our reference measurement file
    #filename_dtu = '1020_50Gs_100us_1545_25_fiber.mat'
    [clock_dtu, signal_dtu] = load_data(filename_dtu)
    ### First we chop the recorded signal to extract individual interferograms
    interferograms_dtu = interferogram_slicing(clock_dtu,signal_dtu)
    ### Fourier spectral processing
    fft_interf_dtu = fft(interferograms_dtu)   
    [amp_array_dtu,phase_array_dtu,linear_term_dtu] = spectral_processing(fft_interf_dtu)
    phase_RF_dtu = phase_array_dtu - linear_term_dtu
    mean_amp_dtu = np.mean(amp_array_dtu, axis = 0)
    mean_phase_dtu = np.mean(phase_RF_dtu, axis = 0)   

    phase_RF_corrected = mean_phase_dtu - mean_phase # Substract the reference measurement
 
    ############################################################################
    #####################   Figures Reference Measurement  #####################
    # Show the first chopped interferogram
    fig_interf = plt.figure() 
    plt.plot(interferograms[0,:]) 
    plt.title('First interferogram')
    # Showing the RF Comb with the Ideal Peaks mapped for a single interferogram
    fig_RFComb = plt.figure()
    plt.plot(freq,abs(fft_interf[1,:]))
    plt.xlim(-2.5e9,2.5e9)
    plt.title('RF Spectrum')
    # Show the RF combs obtained (for each interferogram)
    fig_allRFC = plt.figure()
    plt.plot(freq,abs(np.transpose(fft_interf[:,:])))    
    plt.title('All RF combs')  
    # Show the first interferogram with the peaks recognized
    power_RFSpectrum = 10*np.log10(abs(fft_interf[1,:])**2) #Power Spectrum 
    power_peaks = power_RFSpectrum[pks_idx] # Only for ONE interferograms
    fig_peaks_recogn = plt.figure()
    plt.plot(freq,power_RFSpectrum)
    plt.plot(freq_peaks,power_peaks,'-o')
    plt.xlim(-2.5e9,2.5e9)#,plt.ylim(20,38)
    plt.title('Peaks detected in a single interferogram')
    plt.xlabel('Frequency (THz)')
    plt.ylabel('Power(dB)')
    # Show the Amplitude recovered    
    fig_meanAmp = plt.figure()
    plt.plot(freq_peaks, 20*np.log10(mean_amp),'-o')
    #plt.plot(freq_peaks, 20*np.log10(mean_amp_dtu),'-o')
    plt.title('Mean Amplitude Spectrum')
    # Show the phase sets of each interferogram processed      
    fig_allPhases = plt.figure()
    plt.plot(freq[pks_idx],np.transpose(phase_array[:,:]))
    plt.title('Phase set retrieved for each interferogram')       
    plt.xlabel('Frequency [Hz]')
    plt.ylabel('Phase [rad]')
    # Show the phase sets without the linear term    
    fig_noLinear = plt.figure()
    plt.plot(freq_peaks, np.transpose(phase_RF[:,:]))
    plt.title('Phase removing the linear term'), plt.xlabel('Frequency'), plt.ylabel('Phase(rad)')
    # Show the Phase recovered
    fig_MeanPhase = plt.figure()
    plt.plot(freq_peaks, mean_phase, 'go',label ='Mean Phase reference')
    plt.plot(freq_peaks, phase_RF_corrected,label='Mean Phase signal',linestyle=':', marker='o')
    plt.legend()
    plt.title('Mean phase removing the linear term')
    plt.xlabel('Frequency [Hz]'), plt.ylabel('Phase [rad]')
    ############## you are done
    print("--- %s seconds ---" % (time.time() - start_time))
    ## Figures I want to show
    fig_interf.show()       # Show the first interferogram extracted
    fig_RFComb.show()       # RF Comb of the first interferogram (linear)
    fig_allRFC.show()       # All RF Combs processed (linear)
    fig_peaks_recogn.show() # RF Spectrum with ideal peaks selected
    fig_allPhases.show()    # All Phase profiles of interferograms processed
    fig_meanAmp.show()      # Mean Amplitude
    fig_noLinear.show()     # All phase profiles removing linear term
    fig_MeanPhase.show()    # Mean Phase

    input()
    










