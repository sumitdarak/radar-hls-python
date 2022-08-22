all: create_proj


create_proj: fft_accel

fft_accel:
	@echo "Create build folder"
	mkdir -p build
	@echo $(PWD)
	@echo "Copy files"
	cp -f $(PWD)/tcl/run_fft_hls.tcl $(PWD)/build
	@echo "go to build dir"
	cd $(PWD)/build	&&	vitis_hls -f run_fft_hls.tcl


create_proj: cfar_accel


cfar_accel:
	@echo "Create build folder"
	mkdir -p build
	@echo $(PWD)
	@echo "Copy files"
	cp -f $(PWD)/tcl/run_cfar_hls.tcl $(PWD)/build
	@echo "go to build dir"
	cd $(PWD)/build	&&	vitis_hls -f run_cfar_hls.tcl



create_proj: conv2d_accel

conv2d_accel:
	@echo "Create build folder"
	mkdir -p build
	@echo $(PWD)
	@echo "Copy files"
	cp -f $(PWD)/tcl/run_conv2d_hls.tcl $(PWD)/build
	@echo "go to build dir"
	cd $(PWD)/build	&&	vitis_hls -f run_conv2d_hls.tcl

clean:
	rm -rf build/ sim_files/

.PHONY: all create_proj clean