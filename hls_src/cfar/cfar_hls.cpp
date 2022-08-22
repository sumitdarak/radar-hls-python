#include "cfar_hls.h"




void CFAR_TOP(stream<stream_1ch> &in_stream, stream<stream_1ch> &out_stream, uint32_t z_coef)
{
	wrapper_cfar<uint32_t>(in_stream, out_stream, z_coef);

}



void SORT_TOP(uint32_t x[REFWIND], uint32_t y[REFWIND])
{
	sort_insert<uint32_t>(x, y);

}
