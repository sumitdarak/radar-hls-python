import numpy as np
import os
import matplotlib.pyplot as plt



# RMS quantization level
Nq      = 15
SNR_q   = 6.02*Nq - 1.76


def revBits(n, no_of_bits):
    result = 0
    for i in range(no_of_bits):
        result <<= 1
        result |= n & 1
        n >>= 1
    return result

    
def butter_time_int16(x, y, w):
    """
    Radix 2 butterfly implementation (decimation-in-time) with rounding to int16 and scaling factor for output
    @param complex int16 x:  FFT complex point (sample)
    @param complex int16 y:  FFT complex point (sample)
    @param complex int16 w:  Complex coefficient
    @return x_t, y_t complex  int16 samples
    """
    y_w = np.round((y * w)/2**15) # rounding back to 16 bits after multiplication
    y_t = x - y_w
    x_t = x + y_w
    return x_t/2, y_t/2

## twiddling coefficients
def coef_init(Npoints):
    """
    Twiddling coefficients generation 
    @param int Npoints:  length of FFT
    @return wk_16 complex int16 coefficient array
    """
    wk_16 = np.zeros(Npoints, dtype='complex')
    for k in range(Npoints):
        wk_16[k] = np.round((2**15-1) * np.exp(-1j*2*np.pi*k / Npoints))
    return wk_16
    
def fft_dit(x, w):
    """
    FFT decimation-in-time implementation 
    @param  complex int16 x: input data
    @param  complex int16 w: twiddling coefficients
    @return complex int16 y: output data
    """
    Np = np.size(x)
    Ns = int(np.log2(Np))
    Y_int16 = np.zeros((Ns + 1, Np), dtype=complex)
    
    # reverse input
    for k in range(Np):
        Y_int16[0, k] = x[revBits(k, Ns)]
        
        
    for casc in range(Ns):
        d = 0
        for k in range(Np // 2):
            idx_w = int(np.mod(k, 2**casc))*2**(Ns - 1 - casc)
            if np.mod(k, 2**casc) == 0:
                d = 2*k
            idx_1 = d 
            idx_2 = d + 2**casc
            d = d + 1
            Y_int16[casc + 1, idx_1], Y_int16[casc + 1, idx_2] = butter_time_int16(Y_int16[casc, idx_1], Y_int16[casc, idx_2], w[idx_w])
   
    return np.round(Y_int16[Ns, :])



def main():
    print('<< DFT / FFT math modelling')
    print('<< aleksei.rostov@protonmail.com')
    
    curr_path = os.getcwd()
    if (curr_path[-16:] != 'radar-hls-python'):
        print("<< Error! Please change directory!")
        exit()
    
    if not (os.path.exists(curr_path + '/sim_files')):
        print("<< Error! Please run signal_generator.py first!")
        exit()
    
    hw      = np.loadtxt(os.getcwd() + "/sim_files/cmpx_hls.txt", dtype=int)
    hw_re   = hw[0::2]
    hw_im   = hw[1::2]
    hw_cmpx = hw_re + 1j*hw_im
        
    x_re    = np.loadtxt(os.getcwd() + "/sim_files/scaled_re.txt", dtype=int)
    x_im    = np.loadtxt(os.getcwd() + "/sim_files/scaled_im.txt", dtype=int)
    xcmpx   = x_re + 1j*x_im
    Np = np.size(xcmpx)
    Nstages = int(np.log2(Np))
    Nb      = Nstages
        
    u_re    = np.loadtxt(os.getcwd() + "/sim_files/nonscaled_re.txt", dtype=float)
    u_im    = np.loadtxt(os.getcwd() + "/sim_files/nonscaled_im.txt", dtype=float)
    u       = u_re + 1j*u_im
    
    print("<< Coefficients initialization")
    wk_16 = coef_init(Np)
    print("<< FFT computing")
    # print("""
    # """)
    py_cmpx  = fft_dit(xcmpx, wk_16)     # normal output order
    uFFT     = np.fft.fft(u) / Np
    # exit()
    print("<< Plotting results")
    plt.figure(num=1, figsize=(10,10))
    
    plt.subplot(211)
    plt.plot(20*np.log10(np.abs(uFFT))        , 'x-g', label='Numpy FFT')
    plt.plot(20*np.log10(np.abs(py_cmpx/2**15)), 'o-r', label='SW    FFT')
    plt.plot(20*np.log10(np.abs(hw_cmpx/2**15)), '.-b', label='HLS   FFT')
    plt.plot(-SNR_q*np.ones(Np), '--m', label='RMS quantization level')
    plt.plot((-SNR_q - 10*np.log10(Np))*np.ones(Np), '--k', label='FFT noise floor')
    plt.legend()
    plt.title("FFT {} points, RMS quantization level is 92.06 dB, FFT noise floor is {:3.2f} dB".format(Np, SNR_q + 10*np.log10(Np)), fontweight="bold", fontsize=14)
    plt.xlabel('bin')
    plt.ylabel('power, dB')
    plt.grid()
    

    
    plt.subplot(212)
    plt.plot(np.abs(py_cmpx) - np.abs(hw_cmpx), '.-r')
    plt.title("SW FFT - HLS FFT", fontweight="bold", fontsize=14)
    plt.xlabel('bin')
    plt.ylabel('Error')
    plt.grid()
    
    plt.tight_layout()
    plt.show()
    print('<< End Modelling')


if __name__ == "__main__":
    main()