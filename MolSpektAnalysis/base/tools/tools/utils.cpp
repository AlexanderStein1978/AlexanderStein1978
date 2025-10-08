//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//


#include "utils.h"

#include <QString>


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
