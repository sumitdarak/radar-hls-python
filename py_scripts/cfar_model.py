import numpy as np
import matplotlib.pyplot as plt
from decimal import Decimal
import sys
import os

# greater or equal with bool output 
comp_geb = lambda a,b: 0 if a >= b else 1

MAX_VALUE = 2**32


def PE(din, staticVar, k):
    """
    Processing Element - part of systolic array
    @param uint32_t din: new input 
    @param uint32_t staticVar: array for static var
    @param uint32_t k: index
    @return dout uint32_t
    """   
    if din > staticVar[k]: # 
        dout = staticVar[k] # 
        staticVar[k] = din
    else:
        dout = din
    return dout

def sort_insertion(x):
    """
    Sorting network (insertion algorithm)
    @param uint32_t x: input array
    @return y uint32_t: output array
    """  
    # length of input array
    N = np.size(x)
    y = np.zeros(N)
    d = np.zeros(N)
    stVar = np.zeros(N)
    
    for k in range(2*N):
        if k < N:
            d[0] = PE(x[k], stVar, 0)
        else:
            d[0] = PE(MAX_VALUE, stVar, 0)

        for k_th in range(N-1):
            d[k_th + 1] = PE(d[k_th],   stVar, k_th + 1)

        if k > N - 1:
            y[k - N] = d[N-1]
            
    return y


def cfar_1d(s, GrdHalf, RefHalf, T, k_th, Type):
    """CFAR for one-dimentional array.

    Args:
    s (numpy.ndarray)	: Input array.
            GrdHalf(int): Half of guard window.
            RefHalf(int): Half of reference window.
                T(float): Scaling factor for thresholding.
               k_th(int): cell for ordered-statistic
            Type(String): OS - ordered statistic
                          CA - cell-averaging
                          GO - greatest-of
                          SO - smallest-of

    Returns:
        s_cfar(numpy.ndarray): output array of thresholding.
        t_cfar(numpy.ndarray): output array of detection.
    """
    if s.ndim > 1:
        print('numpy.ndarray error')
        return -1
    N 			= np.size(s)
    s_cfar 		= np.empty(N)
    t_cfar 		= np.zeros(N)
    s_cfar[:] 	= np.nan 
    ref_win 	= np.zeros(2*RefHalf)
    
    for idx in range(GrdHalf+RefHalf, N-(GrdHalf+RefHalf)):
        ref_win[0:RefHalf] = s[idx-(GrdHalf + RefHalf):idx-GrdHalf]
        ref_win[RefHalf :] = s[idx+ GrdHalf + 1: idx + GrdHalf + RefHalf + 1]
        ref_win 	= sort_insertion(ref_win)
        Z 			= T*ref_win[k_th]
            
        s_cfar[idx + 0] = Z
        t_cfar[idx + 1] = comp_geb(Z, s[idx + 1])
    return s_cfar, t_cfar
    
    
    
def main():
    print('<< Ordered-Statistic CFAR  math modelling')
    print('<< aleksei.rostov@protonmail.com')
    
    GrdHalf     = 0
    Pfa 	    = 1e-4
    
    curr_path   = os.getcwd()
    if (curr_path[-16:] != 'radar-hls-python'):
        print("<< Error! Please change directory!")
        exit()
    
    if not (os.path.exists(curr_path + '/sim_files')):
        print("<< Error! Please run signal_generator.py first!")
        exit()
        
        
   
    py_param    = np.loadtxt(os.getcwd() + "/sim_files/cfarPy_param.txt", dtype=float)
    RefHalf     =   int(py_param[1] // 2)
    k_th        =   int(py_param[2])
    T           =   int(py_param[3])
    PFA         = float(py_param[4])
    T_u16       = T*2**10
    print("<< Coefficient is {:3.2f}".format(T))
    
    hls_in      = np.loadtxt(os.getcwd() + "/sim_files/cfarIn_u16.txt", dtype=int)
    N           = np.size(hls_in)
    
    hls_out     = np.loadtxt(os.getcwd() + "/sim_files/cfarOut_u32.txt", dtype=int)
    
    x_in        = np.loadtxt(os.getcwd() + "/sim_files/cfarIn_float.txt", dtype=float)
    

    z_cfarOS, t_cfarOS 	= cfar_1d(x_in, GrdHalf, RefHalf, T, k_th, "OS")
    indx        = np.array(np.where(t_cfarOS == 1))
    indx        = indx.T
    
    z0_cfarOS, t0_cfarOS 	= cfar_1d(hls_in, GrdHalf, RefHalf, T_u16, k_th, "OS")
    
    nt = 0
    for k in range(np.size(z0_cfarOS)):
        if(hls_out[k] < hls_in[k] and hls_out[k] > 0):
            nt += 1
            
    print(PFA)
    plt.figure(2)
    plt.plot(10*np.log10(x_in),     '.-r', label="Input Data")
    plt.plot(10*np.log10(z_cfarOS), 'b', linestyle='--', label="Adaptive Threshold")
    plt.plot(-8 * np.ones(N), 'k', linestyle='--', label="Constant Threshold")
    plt.plot(indx, 10*np.log10(x_in[indx[:]]), 'xg', label="Targets")
    plt.legend()
    plt.title("Python OS CFAR, window is {} cells, PFA is {}".format(RefHalf*2, PFA))
    plt.xlabel('range, bins')
    plt.ylabel('level, dB')
    plt.ylim([-30, 10])
    plt.grid()
    
    plt.show()
    
    return
    

    plt.figure(num=1, figsize=(10,10))
    plt.subplot(211)
    plt.plot(10*np.log10(x_in),     '.-r', label="Input Data")
    plt.plot(10*np.log10(z_cfarOS), '.-b', label="OS CFAR Threshold")
    plt.plot(indx, 10*np.log10(x_in[indx[:]]), 'xg', label="detected")
    plt.legend()
    plt.title("Python OS CFAR, {} targets detected".format(np.size(indx)))
    plt.xlabel('bins')
    plt.ylabel('dB')
    plt.ylim([-30, 10])
    plt.grid()
    
    plt.subplot(212)
    plt.plot(10*np.log10(hls_in), '.-r', label="Uint16 Input Data")
    plt.plot(10*np.log10(z0_cfarOS / 2**10), '.-b', label="Python OS CFAR Threshold")
    plt.title("HLS Simulation Output - OS CFAR, {} targets detected".format(nt))
    plt.plot(10*np.log10(hls_out), '.-g', label="HLS OS CFAR Threshold")
    plt.legend()
    plt.grid()
    
    plt.show()
    print('<< End Modelling')





if __name__ == "__main__":
    main()
