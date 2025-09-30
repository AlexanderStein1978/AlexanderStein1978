//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//


#ifndef DOUBLEVECTOR
#define DOUBLEVECTOR


#include "flexvector.h"


class DoubleVector : public FlexVector<double>
{
	public:
		DoubleVector(int start = -499);
		double &operator[](const int &i); 
};

#endif
