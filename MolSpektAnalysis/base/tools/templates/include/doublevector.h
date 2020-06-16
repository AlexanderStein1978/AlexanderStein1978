//
// C++ Interface: DoubleVector
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2006 - 2019
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
