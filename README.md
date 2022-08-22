# High-Level Synthesis implementation of Radar SIgnal and Data Processing Algorithms

Some common Radar Signal and Data Processing Algorithms Implementation in Vivado HLS with detailed explanation and math modelling (Python).
The lessons show how to do for some basic algorithms:

- mathematical modelling in Python for verification of simulation output and hardware design consideration
- Vivado High-Level Synthesis simulation, implementation and hardware design optimization

The lessons are very useful for understanding Radar Signal and Data processing fundamentals, 
FPGA hardware design and implementation by using advantage of high-level synthesis. 

1. [Ordered-Statistic CFAR (Sorting Algorithm)](./doc/cfar-hls-python.md)
2. [FFT decimation-in-time (Linear Algorithm)](./doc/fft-hls-python.md)


## Structure of the repository
Folders and files of the repository
```
+-- doc\   
|   |
|   +-- fft-hls-python.md   -- main documentation 
|   +-- cfar-hls-python.md  -- main documentation
|	+-- images\
|
+-- py_scripts\   
|   |
|   +-- fft_generator.py 	-- complex data generator for RTL simulation and fft_model.py
|   +-- fft_model.py        -- math model of FFT DIT algorithm for int16
|   +-- cfar_generator.py 	-- data generator for RTL simulation and cfar_model.py
|   +-- cfar_model.py       -- math model of OS-CFAR algorithm
|
+-- hls_src\                -- Vivado HLS sources
|   |
|	+-- fft\                
|   |	|
|   |	+-- fft_hls.cpp
|   |	+-- fft_hls.h
|   |	+-- fft_hls_tb.cpp
|   |	+-- fft_hls_tb.h
|	|
|	+-- cfar\
|   |   |
|   |	+-- cfar_hls.cpp
|   |	+-- cfar_hls.h
|   |	+-- cfar_hls_tb.cpp
|   |	+-- cfar_hls_tb.h
|
+-- tcl\              		-- tcl script for compiling Vivado HLS
|   |
|   +-- run_cfar_hls.tcl
|   +-- run_fft_hls.tcl
|
+-- Makefile                -- Makefile for building Vivado HLS project
|
+-- README.md
|
+-- LICENSE
```

How to work with the repository

1. Create synthetic signal with the <i>fft_generator.py</i> script

```sh
python3 py_scripts/fft_generator.py 1024 3 40
```

2. In the repository root directory run <i>make</i> command for building and launching Vivado HLS project

```sh
make fft_accel
```

3. Run <i>fft_model.py</i> script for cheking Vivado HLS simulation result

```sh
python3 py_scripts/fft_model.py
```