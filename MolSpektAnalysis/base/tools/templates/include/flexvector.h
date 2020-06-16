//
// C++ Interface: FlexVector
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2006 - 2019
//
// Copyright: See README file that comes with this source code
//
//


#ifndef FLEXVECTOR
#define FLEXVECTOR


#include "velement.h"


template <class T> class FlexVector
{
	public:
		FlexVector(int start = -499);
		~FlexVector();
		T &operator[](const int &i);
		T *getData();
		void clear();
		int getFirstIndex();
		int getLastIndex();
		
	protected:
		VElement<T> *CurElement;
		int FI, LI;
};

template <class T> FlexVector<T>::FlexVector(int start)
{
	CurElement = new VElement<T>;
	CurElement->next = CurElement->prev = 0;
	CurElement->startI = start;
	CurElement->endI = start + 999;
	FI = LI = start;
}

template <class T> FlexVector<T>::~FlexVector()
{
	clear();
}

template <class T> void FlexVector<T>::clear()
{
	while (CurElement->next != 0) CurElement = CurElement->next;
	while (CurElement->prev != 0) 
	{
		CurElement = CurElement->prev;
		delete CurElement->next;
	}
	delete CurElement;
}

template <class T> T *FlexVector<T>::getData()
{
	int i, j, L = FI - LI + 1;
	T *R = new T[L];
	while (CurElement->startI > FI) CurElement = CurElement->prev;
	for (i=0, j = FI - CurElement->startI; j <= LI; i++, j++) 
	{
		if (j == 1000) 
		{
			CurElement = CurElement->next;
			j=0;
		}
		R[i] = CurElement->Data[j];
	}
	return R;
}

template <class T> int FlexVector<T>::getFirstIndex()
{
	return FI;
}

template <class T> int FlexVector<T>::getLastIndex()
{
	return LI;
}

template <class T> T &FlexVector<T>::operator[](const int &i)
{
	while (CurElement->startI > i)
	{
		if (CurElement->prev == 0)
		{
			CurElement->prev = new VElement<T>;
			CurElement->prev->next = CurElement;
			CurElement->prev->prev = 0;
			CurElement->prev->startI = CurElement->startI - 1000;
			CurElement->prev->endI = CurElement->startI - 1;
		}
		CurElement = CurElement->prev;
	}
	while (CurElement->endI < i)
	{
		if (CurElement->next == 0)
		{
			CurElement->next = new VElement<T>;
			CurElement->next->prev = CurElement;
			CurElement->next->next = 0;
			CurElement->next->startI = CurElement->endI + 1;
			CurElement->next->endI = CurElement->endI + 1000;
		}
		CurElement = CurElement->next;
	}
	return CurElement->Data[i - CurElement->startI];
}

#endif
