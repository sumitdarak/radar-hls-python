#include "ap_axi_sdata.h"
#include "hls_stream.h"
#include "parameters.h"

using namespace hls;

/* ****************************** DEFINES ************************************** */

#define MAX_VALUE      		0x1FFFFFFF


typedef ap_axiu<32, 0, 0, 0> stream_1ch;

/* ****************************** STRUCTURES *********************************** */


/* ****************************** FUNCTIONS ************************************ */


/* ****************************** TOP FUNCTIONS DECLARATION ******************** */

void CFAR_TOP(stream<stream_1ch> &in_stream, stream<stream_1ch> &out_stream, uint32_t z_coef);

void SORT_TOP(uint32_t x[REFWIND], uint32_t y[REFWIND]);

/* ****************************** C++ TEMPLATES ********************************* */



/*
 * 	reading a sample from axis interface
 */
template <typename T>
T read_stream(stream_1ch const &e)
{
#pragma HLS INLINE
	union
	{
		int ival;
		T oval;
	} converter;
	converter.ival = e.data;
	T ret = converter.oval;

	return ret;
} // read_stream

template <typename T>
void pop_input(stream<stream_1ch> &in_stream,T y[NPOINTS])
{
// #pragma HLS INLINE
	pin_L:for(uint16_t idx = 0; idx <  NPOINTS; idx ++)
		y[idx] = read_stream<T>(in_stream.read());
} // pop_input


/*
 * 	writing a sample to axis interface
 */
template <typename T>
stream_1ch write_stream(T const &v, bool last = false)
{
#pragma HLS INLINE
	stream_1ch e;

	union
	{
		int oval;
		T ival;
	} converter;
	converter.ival = v;
	e.data = converter.oval;
	e.last = last ? 1 : 0;

	return e;
}

template <typename T>
void push_output(stream<stream_1ch> &out_stream,T y[NPOINTS])
{
// #pragma HLS INLINE
	pout_L:for(uint16_t idx = 0; idx <  NPOINTS; idx ++)
		out_stream.write(write_stream<T>(y[idx], (idx == NPOINTS - 1)));
} // push_output


template <typename T>
void pe(T din, T *staticvar, T *dout)
{

	if(din > *staticvar)
	{
		*dout = *staticvar;
		*staticvar = din;
	}
	else
		*dout = din;

}

template <typename T>
void sort_insert(T x[REFWIND], T y[REFWIND])
{
	static T d[REFWIND];
	static T s[REFWIND];
	for(char idx = 0; idx < REFWIND; idx ++)
	{
		d[idx] = 0;
		s[idx] = 0;
	}

	for(T idx = 0; idx < 2*REFWIND; idx ++)
	{
		if(idx < REFWIND)
			pe<T>(x[idx], &s[0], &d[0]);
		else
			pe<T>(MAX_VALUE, &s[0], &d[0]);

		for (T idx_2 = 0; idx_2 < REFWIND-1; idx_2 ++)
			pe<T>(d[idx_2], &s[idx_2 + 1], &d[idx_2 + 1]);

		if(idx >= REFWIND)
			y[idx - REFWIND] = d[REFWIND-1];
	} // for
}



template <typename T>
void cfar(T x[NPOINTS], T y[NPOINTS], T z_coef)
{
	static T x_sort[REFWIND];
	static T y_sort[REFWIND];

	for (uint16_t idx =  0; idx < NPOINTS; idx ++)
	{

		if(idx >=  REFWIND / 2 && idx < (NPOINTS - REFWIND / 2))
		{
			L1_cfar:for(char t1 = 0; t1 < REFWIND / 2; t1 ++) 			// before CUT
				x_sort[t1] = x[idx - t1 - 1];

			L2_cfar:for(char t1 = 0; t1 < REFWIND / 2; t1 ++) 			// after CUT
				x_sort[t1 + REFWIND / 2] = x[idx + t1];

			sort_insert<T>(x_sort, y_sort);

			T tmp_0  	= z_coef  * y_sort[KTH_CELL];
			y[idx] 		= tmp_0 >> 10;
		}
		else
			y[idx] = 0;


	}

}


template <typename T>
void wrapper_cfar(stream<stream_1ch> &in_stream, stream<stream_1ch> &out_stream, T z_coef)
{
	T x[NPOINTS], y[NPOINTS];
	pop_input  	<T>( in_stream, x);
	cfar		<T>(x, y, z_coef);
	push_output <T>(out_stream, y);
}
