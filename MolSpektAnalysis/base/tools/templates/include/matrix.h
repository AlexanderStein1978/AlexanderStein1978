//
// C++ Interface: Matrix
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2006 - 2019
//
// Copyright: See README file that comes with this source code
//
//

#ifndef MATRIX
#define MATRIX


template <class T> struct Matrix
{
	Matrix(const int &i, const int &j)
	{
		R = new T*[i];
		int n;
		for (n=0; n<i; n++) R[n] = new T[j];
		nr = i;
	}
	~Matrix()
	{
		int n;
		for (n=0; n<nr; n++) delete[] R[n];
		delete[] R;
	}
	int nr;
	T **R;
};

#endif
