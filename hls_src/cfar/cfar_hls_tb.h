#include "cfar_hls.h"
#include <stdint.h>
#include <algorithm>
#include <math.h>
#include <stdlib.h>
#include <memory>
#include <cstdlib>
#include <fstream>

using namespace std;

#define DIN 	"..\\..\\..\\..\\..\\sim_files\\cfarIn_u16.txt"
#define DOUT    "..\\..\\..\\..\\..\\sim_files\\cfarOut_u32.txt"


void read_txt(uint32_t *dout)
{
	// array
	std::cout << "READING INPUT" << std::endl;
	uint32_t _dout[NPOINTS];
	uint16_t idx = 0;

	ifstream fpIn (DIN);

	if (fpIn.is_open())
	{
		idx = 0;
		while ( !fpIn.eof())
		{
			fpIn >> _dout[idx];
			idx++;
		}
		fpIn.close();
	}

	else cout << "Unable to open file" << endl;

	memcpy(&dout[0], &_dout, NPOINTS*sizeof(uint32_t));

	std::cout  << "DONE!" << std::endl;

}



void TEST_SORT()
{

	cout << "\n<< Start TEST_SORT \n";

	uint32_t      x[REFWIND];
	uint32_t   y_hw[REFWIND];
	uint32_t y_gold[REFWIND];

	srand(20);
	for(unsigned int k = 0; k < REFWIND; k ++)
	{
		     x[k] = rand() % 1024;
		  y_hw[k] = 0;
		y_gold[k] = 0;
	}

	memcpy(&y_gold[0], &x[0], REFWIND*sizeof(uint32_t));
    // C++ std::sort() algorithm
	sort(std::begin(y_gold), std::end(y_gold));

	SORT_TOP(x, y_hw);

	cout << "\n<< Comparing output with sort() algorithm \n";
	uint32_t error = 0;
	for(uint16_t idx = 0; idx < REFWIND; idx ++)
		error += abs(y_gold[idx] - y_hw[idx]);


	cout << "\n<< Error is " << error << endl;

	cout << "\n<< End TEST_SORT \n" << endl;

}



void TEST_CFAR()
{

	cout << "\n<< Start TEST_CFAR \n" << endl;

	uint32_t      x[NPOINTS];
	uint32_t      y[NPOINTS];
	read_txt(&x[0]);

	stream<stream_1ch> 	src, dst;

	for(unsigned int idx = 0; idx < NPOINTS; idx ++)
	{
	   stream_1ch VallIn;
	   VallIn.data = x[idx];
	   VallIn.keep = -1;
	   VallIn.last = (idx == NPOINTS - 1);
	   src << VallIn;
	}


	CFAR_TOP(src, dst, (uint32_t)Z_COEF);


	for(unsigned int idx = 0; idx < NPOINTS; idx ++)
	{
		stream_1ch VallOut;
		dst >> VallOut;
		y[idx] = VallOut.data;
	}

	cout << endl << "WRITING RESULT TO TXT FILE" << endl;

	ofstream myfile;
	myfile.open (DOUT);
	for(int k = 0; k < NPOINTS; k++)
	{
		myfile << y[k]  << "\n";

	}
	myfile.close();
	cout  << "DONE!" << endl;

	cout << "\n<< End TEST_CFAR \n" << endl;

}
