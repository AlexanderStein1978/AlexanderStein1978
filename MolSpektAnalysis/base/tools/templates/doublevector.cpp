//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//


#include "doublevector.h"


DoubleVector::DoubleVector(int start) : FlexVector<double>(start)
{
	int i;
	for (i=0; i< 1000; i++) CurElement->Data[i] = 0.0;
}

double &DoubleVector::operator[](const int &i)
{
	while (CurElement->startI > i)
	{
		if (CurElement->prev == 0)
		{
			int j;
			CurElement->prev = new VElement<double>;
			CurElement->prev->next = CurElement;
			CurElement->prev->prev = 0;
			CurElement->prev->startI = CurElement->startI - 1000;
			CurElement->prev->endI = CurElement->startI - 1;
			for (j=0; j < 1000; j++) CurElement->prev->Data[j] = 0.0;
		}
		CurElement = CurElement->prev;
	}
	while (CurElement->endI < i)
	{
		if (CurElement->next == 0)
		{
			int j;
			CurElement->next = new VElement<double>;
			CurElement->next->prev = CurElement;
			CurElement->next->next = 0;
			CurElement->next->startI = CurElement->endI + 1;
			CurElement->next->endI = CurElement->endI + 1000;
			for (j=0; j < 1000; j++) CurElement->next->Data[j] = 0.0;
		}
		CurElement = CurElement->next;
	}
	return CurElement->Data[i - CurElement->startI];
}
