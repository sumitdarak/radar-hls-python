#include "fft_hls.h"


/* ****************************** DEFINES ************************************** */

// #define CONSOLE // output to console (print out to console for FFT points less or equal  16 )
#define DIN_RE 	"..\\..\\..\\..\\..\\sim_files\\scaled_re.txt"
#define DIN_IM 	"..\\..\\..\\..\\..\\sim_files\\scaled_im.txt"
#define DOUT    "..\\..\\..\\..\\..\\sim_files\\cmpx_hls.txt"

using namespace std;

/* ****************************** C++ TEMPLATES ************************************** */

void read_txt(uint32_t *dout)
{
	// array
	std::cout << "READING INPUT" << std::endl;
	cmpx_t<int16_t> _dout[NPOINTS];
	uint16_t idx = 0;

	ifstream fpReIn (DIN_RE);
	ifstream fpImIn (DIN_IM);

	if (fpReIn.is_open())
	{
		idx = 0;
		while ( !fpReIn.eof())
		{
			fpReIn >> _dout[idx].re_;
			if(abs(_dout[idx].re_) > 32767)
			{
				std::cout << "ERROR REAL INPUT \n";
				throw std::exception();
			}

			idx++;
		}
		fpReIn.close();
	}

	else cout << "Unable to open file" << endl;


	if (fpImIn.is_open())
	{
		idx = 0;
		while ( !fpImIn.eof())
		{
			fpImIn >> _dout[idx].im_;
			if(abs(_dout[idx].im_) > 32767)
			{
				std::cout << "ERROR IMAG INPUT \n";
				throw std::exception();
			}

			idx++;
		}
		fpImIn.close();
	}

	else std::cout << "Unable to open file" << std::endl;

	memcpy(&dout[0], &_dout, NPOINTS*sizeof(uint32_t));

	std::cout  << "DONE!" << std::endl;

}



/**
 *
 *  TOP BUTTERFLY TESTING
 *
 */

void test_butterfly()
{
	uint32_t x0 = 0x03458755;
	uint32_t y0 = 0x0234F435;
	uint32_t w0 = 0xF4359985;

	uint32_t x1 = 0, y1 = 0;

	cout  << endl << "START BUTTERFLY SIMULATION" << endl << endl;
	BUTTERFLY_TOP(x0, y0, w0, &x1, &y1);

	uint2cmpx.uint = x0;
	cmpx_t<int16_t> x_0 = uint2cmpx.cmpx;

	uint2cmpx.uint = y0;
	cmpx_t<int16_t> y_0 = uint2cmpx.cmpx;

	uint2cmpx.uint = w0;
	cmpx_t<int16_t> w_0 = uint2cmpx.cmpx;

	cout << setw(6) << "x0.re = " << x_0.re_ << "\t";
	cout << setw(6) << "x0.im = " << x_0.im_ << "\n";
	cout << setw(6) << "y0.re = " << y_0.re_ << "\t";
	cout << setw(6) << "y0.im = " << y_0.im_ << "\n";
	cout << setw(6) << "w0.re = " << w_0.re_ << "\t";
	cout << setw(6) << "w0.im = " << w_0.im_ << "\n";

	uint2cmpx.uint = x1;
	cmpx_t<int16_t> x_1 = uint2cmpx.cmpx;

	uint2cmpx.uint = y1;
	cmpx_t<int16_t> y_1 = uint2cmpx.cmpx;

	cout  << endl << "HARDWARE OUTPUT" << endl << endl;

	cout << setw(6) << "x1.re = " << x_1.re_ << "\t";
	cout << setw(6) << "x1.im = " << x_1.im_ << "\n";
	cout << setw(6) << "y1.re = " << y_1.re_ << "\t";
	cout << setw(6) << "y1.im = " << y_1.im_ << "\n";

	/*
	 *
	 *   x_1 = x_0 + y_0 * w_0
	 *   y_1 = x_0 - y_0 * w_0
	 *
	 */

	x_1 = {0, 0};
	y_1 = {0, 0};

	cmpx_t<int32_t> yw = {0, 0};
	yw.re_ = ((int32_t)(y_0.re_ * w_0.re_) - (int32_t)(y_0.im_ * w_0.im_)) >> 15;
	yw.im_ = ((int32_t)(y_0.im_ * w_0.re_) + (int32_t)(y_0.re_ * w_0.im_)) >> 15;

	x_1.re_ = (int16_t)(((int32_t)x_0.re_ + yw.re_) >> 1);
	x_1.im_ = (int16_t)(((int32_t)x_0.im_ + yw.im_) >> 1);

	y_1.re_ = (int16_t)(((int32_t)x_0.re_ - yw.re_) >> 1);
	y_1.im_ = (int16_t)(((int32_t)x_0.im_ - yw.im_) >> 1);

	cout  << endl << "GOLD OUTPUT" << endl << endl;

	cout << setw(6) << "x1.re = " << x_1.re_ << "\t";
	cout << setw(6) << "x1.im = " << x_1.im_ << "\n";
	cout << setw(6) << "y1.re = " << y_1.re_ << "\t";
	cout << setw(6) << "y1.im = " << y_1.im_ << "\n";


	cout << endl << "END BUTTERFLY SIMULATION" << endl << endl;


}


/**
 *
 * 	TOP FFT testing
 *
 */

void test_fft()
{
	cout  << endl << "START FFT SIMULATION" << endl;

	uint32_t din_uint[NPOINTS];
	stream<stream_1ch> 	src, dst;

	read_txt(&din_uint[0]);

#ifdef CONSOLE
	cmpx_t<int16_t> din_points = {0, 0};
	cout << endl << "Input Points " << endl;
	cout << setw(3) << "NUM" << "\t";
	cout << setw(6) << "REAL" << "\t";
	cout << setw(6) << "IMAG" << "\n";
	for(int idx = 0; idx < NPOINTS; idx++)
	{
		uint2cmpx.uint = din_uint[idx];
		din_points     = uint2cmpx.cmpx;

		cout << setw(3) << idx << ":" << "\t";
		cout << setw(6) << din_points.re_ << "\t";
		cout << setw(6) << din_points.im_ << "\n";

	}
#endif

	for(unsigned int idx = 0; idx < NPOINTS; idx ++)
	{
	   stream_1ch VallIn;
	   VallIn.data = din_uint[idx];
	   VallIn.keep = -1;
	   VallIn.last = (idx == NPOINTS - 1);
	   src << VallIn;
	}


	cout  << endl << "RUN HARDWARE TEST FOR " << NPOINTS << " FFT points"  << endl;
	FFT_TOP(src, dst);
	cout  << "DONE!" << endl;

	cmpx_t<int16_t> dout_points[NPOINTS];

	for(unsigned int idx = 0; idx < NPOINTS; idx ++)
	{
		stream_1ch VallOut;
		dst >> VallOut;
		uint2cmpx.uint = VallOut.data;
		dout_points[idx].re_ = uint2cmpx.cmpx.re_;
		dout_points[idx].im_ = uint2cmpx.cmpx.im_;
	}

#ifdef CONSOLE
	std::cout << endl << "Output Points " << endl;
	cout << setw(3) << "NUM" << "\t";
	cout << setw(6) << "REAL" << "\t";
	cout << setw(6) << "IMAG" << "\n";
	for(int idx = 0; idx < NPOINTS; idx++)
	{
		cout << setw(3) << idx << ":" << "\t";
		cout << setw(6) << dout_points[idx].re_ << "\t";
		cout << setw(6) << dout_points[idx].im_ << "\n";
	}


#else
	cout << endl << "WRITING RESULT TO TXT FILE" << endl;

	ofstream myfile;
	myfile.open (DOUT);
	for(int k = 0; k < NPOINTS; k++)
	{
		myfile << dout_points[k].re_  << "\n";
		myfile << dout_points[k].im_  << "\n";
	}
	myfile.close();
	cout  << "DONE!" << endl;

#endif

	cout << endl << "END FFT SIMULATION" << endl << endl;
}
