//
// C++ Implementation: utils
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2006 - 2019
//
// Copyright: See README file that comes with this source code
//
//

#include "utils.h"


QString **CreateQString(const int &I, const int &J)
{
	QString **R = new QString*[I];
	int i;
	for (i=0; i<I; i++) R[i] = new QString[J];
	return R;
}


void Destroy(QString **v, const int &I)
{
	int i;
	for (i=0; i<I; i++) delete[] v[i];
	delete[] v;
}
