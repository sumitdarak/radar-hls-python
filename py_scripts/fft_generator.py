import numpy as np
import argparse
import os
import sys


def coef_init(Npoints):
    """
    Twiddling coefficients generation 
    @param int Npoints:  length of FFT
    @return wk_16 complex int16 coefficient array
    """
    wk_16 = np.zeros(Npoints, dtype='complex')
    for k in range(Npoints):
        wk_16[k] = np.round(32767 * np.exp(-1j*2*np.pi*k / Npoints))
    return wk_16

Nq      = 15
# RMS quantization level
SNR_q   = 6.02*Nq + 1.76


parser = argparse.ArgumentParser()
parser.add_argument('Npoints' ,default=None, type=int)
parser.add_argument('Nsignals',default=None, type=int)
parser.add_argument('SNR_dB'  ,default=None, type=int)
args = parser.parse_args()

# example
# python data_generator.py 1024 4 40

def main():
    print('<< Generating Input Data for FFT')
    print('<< aleksei.rostov@protonmail.com')
    
       
    curr_path = os.getcwd()
    if (curr_path[-16:] != 'radar-hls-python'):
        print("<< Error! Please change directory!")
        exit()
        
    if not (os.path.exists(curr_path + '/sim_files')):
        os.makedirs(curr_path + '/sim_files')
   
   
    if args.Npoints is None:
        Npoints = 32
    else:
        Npoints = args.Npoints
        
    if args.Nsignals is None:
        Nsignals = 4
    else:
        Nsignals = args.Nsignals
        
    if args.SNR_dB is None:
        SNR_dB = 0
    else:
        SNR_dB = args.SNR_dB


    Ampl     = 2**15
    x        = np.zeros(Npoints, dtype=complex)
    for k in range(Nsignals):
        nbin = Npoints // 8 * (k + 1)
        x   += 10**((SNR_dB - 5*k)/20) * np.exp(2*1j*np.pi*nbin*np.arange(Npoints)/Npoints) + np.random.randn(Npoints) # 
    
    x       *= 10**(-SNR_q/20)
    # saving non-scaled input signal
    np.savetxt(curr_path + '/sim_files/nonscaled_re.txt', np.real(x),fmt='%f')
    np.savetxt(curr_path + '/sim_files/nonscaled_im.txt', np.imag(x),fmt='%f')
    
    
    x_16     = np.round(Ampl*x)
    
    np.savetxt(curr_path + '/sim_files/scaled_re.txt', np.real(x_16),fmt='%d')
    np.savetxt(curr_path + '/sim_files/scaled_im.txt', np.imag(x_16),fmt='%d')
    
    
    
    w_cmpx16 = coef_init(Npoints)
    
    
    with open(curr_path + '/hls_src/fft/parameters.h', 'w') as fp:
        fp.write("\n")
        fp.write("#define NPOINTS   ")
        fp.write(str(Npoints))
        fp.write("\n#define FFTRADIX  ")
        fp.write(str(int(np.log2(Npoints))))
        fp.write("\n#define FFTRAD_1  ")
        fp.write(str(int(np.log2(Npoints) + 1)))
        fp.write("\n")
        fp.write("\n")
        fp.write("\n")
        fp.write('static uint32_t wcoe[] = {')
        for idx in range(Npoints // 2):
            w_re = int(np.real(w_cmpx16[idx]))
            w_im = int(np.imag(w_cmpx16[idx]))
            fp.write("0x")
            if w_im >= 0:
                h_16 = hex(w_im)
                if   len(h_16) == 2:fp.write('0000')
                elif len(h_16) == 3:fp.write('000')
                elif len(h_16) == 4:fp.write('00')
                elif len(h_16) == 5:fp.write('0')
                
                fp.write(h_16[2:])
            else:
                w_im = 2**16 - abs(w_im)
                h_16 = hex(w_im)
                fp.write(h_16[2:])
                
                
            if w_re >= 0:
                h_16 = hex(w_re)
                if   len(h_16) == 2:fp.write('0000')
                elif len(h_16) == 3:fp.write('000')
                elif len(h_16) == 4:fp.write('00')
                elif len(h_16) == 5:fp.write('0')
                
                fp.write(h_16[2:])
            else:
                w_re = 2**16 - abs(w_re)
                h_16 = hex(w_re)
                fp.write(h_16[2:])
            if(idx < Npoints-1):fp.write(',')
            if(Npoints-1 > idx > 0):
                if((np.mod(idx + 1, 8) == 0)):
                    fp.write("\n")
                    fp.write("                          ")
            
        fp.write("};")
    
    print('<< Successfully Done')
    
    
if __name__ == "__main__":
    main()
