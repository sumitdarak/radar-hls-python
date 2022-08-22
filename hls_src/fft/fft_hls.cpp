#include "fft_hls.h"


void FFT_TOP(stream<stream_1ch> &in_stream, stream<stream_1ch> &out_stream)
{
#pragma HLS DATAFLOW
#pragma HLS INTERFACE ap_ctrl_none port=return
#pragma HLS INTERFACE axis port=in_stream
#pragma HLS INTERFACE axis port=out_stream

	wrapped_fft_hw <uint32_t, int32_t, int16_t>(in_stream, out_stream);
}


void BUTTERFLY_TOP(uint32_t x0, uint32_t y0, uint32_t w0, uint32_t *x1, uint32_t *y1)
{
#pragma HLS PIPELINE
#pragma HLS INTERFACE ap_ctrl_none port=return
	butter_dit<uint32_t, int32_t, int16_t, 15>(x0, y0, w0, x1, y1);

}
