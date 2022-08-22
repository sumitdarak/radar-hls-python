############################################################
## This file is generated automatically by Vitis HLS.
## Please DO NOT edit it.
## Copyright 1986-2022 Xilinx, Inc. All Rights Reserved.
############################################################
open_project -reset cfar_accel
set_top CFAR_TOP
add_files ../hls_src/cfar/parameters.h
add_files ../hls_src/cfar/cfar_hls.cpp
add_files ../hls_src/cfar/cfar_hls.h
add_files -tb ../hls_src/cfar/cfar_hls_tb.cpp
add_files -tb ../hls_src/cfar/cfar_hls_tb.h


open_solution -reset "solution1"
set_part {xczu2cg-sfvc784-1-e}
create_clock -period 10 -name default
#source "./fft_accel/solution1/directives.tcl"
csim_design

exit

csynth_design
cosim_design
export_design -format ip_catalog
