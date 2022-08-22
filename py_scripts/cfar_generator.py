import numpy as np
import argparse
import os
import sys

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx], idx

def cfar_din_generator(N):

    fn_1        = N // 100 * 83
    SNR_dB      = 3
    
    # signal and 1st clutter zone
    s       = 10**(SNR_dB/20) * np.exp(2*1j*np.pi*fn_1*np.arange(N)/N) + np.random.randn(N)
    
    fn_1    = N // 100 * 6
    SNR_dB  = -6
    s       += 10**(SNR_dB/20) * np.exp(2*1j*np.pi*fn_1*np.arange(N)/N) + np.random.randn(N)
    sF      = np.fft.fft(s) / N
    
    # 2nd clutter zone (30 %)
    n       = 10**(20/20) * np.random.randn(N)
    nF      = np.fft.fft(n) / N
    nF[:(N // 100) * 70] = 0
    
    sF     += nF
    
    # 3rd clutter zone (30 %)
    n       = 10**(30/20) * np.random.randn(N)
    nF      = np.fft.fft(n) / N
    nF[:(N // 100) * 30] = 0
    nF[(N // 100) * 45:] = 0
    
    sF     += nF
    xF      = np.abs(sF)
    # output norming
    xF     /= np.max(xF)
  
    
    return xF

parser = argparse.ArgumentParser(description='Input signal size and CFAR parameters')
parser.add_argument('NPOINTS' , help='Number of points', default=None, type=int)
parser.add_argument('REFWIND', help='Number of cells in CFAR window',default=None, type=int)
parser.add_argument('PFA' ,help='Probability of false alarm' ,default=None, type=float)
args = parser.parse_args()
    
def main():
    print('<< Generating Input Data for CFAR')
    print('<< aleksei.rostov@protonmail.com')
    
       
    if args.NPOINTS is None:
        NPOINTS = 32
    else:
        NPOINTS = args.NPOINTS
        
    if args.REFWIND is None:
        REFWIND = 4
    else:
        REFWIND = args.REFWIND
        
    if args.PFA is None:
        PFA = 1e-1
    else:
        PFA = args.PFA
    
    curr_path = os.getcwd()
    if (curr_path[-16:] != 'radar-hls-python'):
        print("<< Error! Please change directory!")
        exit()
        
    if not (os.path.exists(curr_path + '/sim_files')):
        os.makedirs(curr_path + '/sim_files')
        
        
    sF      = cfar_din_generator(NPOINTS)
    
    # k-th cell is 75 % from size of sliding window
    KTH_CELL    = (REFWIND*75) // 100
    
    # Pfa 	    = PFA
    
    N           = REFWIND
    # scaling factor calculating for Pfa
    dPfa_0 = np.zeros(REFWIND)
    for k in range(REFWIND):
        alpha = k + 1
        dPfa_0[k] = (np.math.factorial(N) * np.math.factorial(alpha + N - KTH_CELL)) / (np.math.factorial(N - KTH_CELL) * np.math.factorial(alpha + N))
    
    val, SCALING = find_nearest(dPfa_0, PFA)
    
    py_param = np.zeros(5)
    py_param[0] = NPOINTS
    py_param[1] = REFWIND
    py_param[2] = KTH_CELL
    py_param[3] = SCALING
    py_param[4] = PFA
        
    T_u16       = np.round(SCALING*2**10)
    
    np.savetxt(curr_path + '/sim_files/cfarIn_u16.txt', np.round(2**16*sF),fmt='%d')
    np.savetxt(curr_path + '/sim_files/cfarIn_float.txt', sF,fmt='%f')
    np.savetxt(curr_path + '/sim_files/cfarPy_param.txt', py_param,fmt='%f')
    
    
    with open(curr_path + '/hls_src/cfar/parameters.h', 'w') as fp:
        fp.write("\n")
        fp.write("#define NPOINTS   ")
        fp.write(str(NPOINTS))
        fp.write("\n#define REFWIND   ")
        fp.write(str(REFWIND))
        fp.write("\n#define KTH_CELL  ")
        fp.write(str(KTH_CELL))
        fp.write("\n#define Z_COEF    ")
        fp.write(str(int(T_u16)))
        fp.write("\n")
        fp.write("\n")
        fp.write("\n")
                
                
    print('<< Successfully Done')
     
    
    
if __name__ == "__main__":
    main()
