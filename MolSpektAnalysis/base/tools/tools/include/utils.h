//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//

#ifndef UTILS
#define UTILS

#include <cstring>

#define sqr(x) ((x)*(x))

class QString;


QString **CreateQString(const int &i, const int &j);
void Destroy(QString **v, const int &i);

inline double **Create(const int &I, const int &J, bool initialize = false)
{
	double **R = new double*[I];
	int i;
	for (i=0; i<I; i++)
    {
        R[i] = new double[J];
        if (initialize) memset(R[i], 0, J * sizeof(double));
    }
	return R;
}

inline int **CreateInt(const int &I, const int &J)
{
	int **R = new int*[I];
	int i;
	for (i=0; i<I; i++) R[i] = new int[J];
	return R;
}

inline int ***CreateInt(const int &I, const int &J, const int &K)
{
	int ***R = new int**[I];
	int i, j;
	for (i=0; i<I; i++)
	{
		R[i] = new int*[J];
		for (j=0; j<J; j++) R[i][j] = new int[K];
	}
	return R;
}

inline bool **CreateBool(const int &I, const int &J)
{
	bool **R = new bool*[I];
	int i;
	for (i=0; i<I; i++) R[i] = new bool[J];
	return R;
}

inline double **Create1(const int &I, const int &J)
{
	double **R = new double*[I];
	int i;
	for (i=0; i<I; i++) R[i] = (new double[J]) - 1;
	return R - 1;
}

inline double ***Create(const int &I, const int &J, const int &K)
{
	double ***R = new double**[I];
	int i;
	for (i=0; i<I; i++) R[i] = Create(J, K);
	return R;
}

inline bool ***CreateBool(const int &I, const int &J, const int &K)
{
	bool ***R = new bool **[I];
	int i;
	for (i=0; i<I; i++) R[i] = CreateBool(J, K);
	return R;
}

inline double ****Create(const int &I, const int &J, const int &K, const int &L)
{
	double ****R = new double ***[I];
	int i;
	for (i=0; i<I; i++) R[i] = Create(J, K, L);
	return R;
}

inline bool ****CreateBool(const int &I, const int &J, const int &K, const int &L)
{
	bool ****R = new bool ***[I];
	int i;
	for (i=0; i<I; i++) R[i] = CreateBool(J, K, L);
	return R;
}

inline double***** Create(const int& i, const int& j, const int& k, const int& l, const int& m)
{
	double *****R = new double****[i];
	int n;
	for (n=0; n<i; n++) R[n] = Create(j, k, l, m);
	return R;
}

inline bool***** CreateBool(const int& i, const int& j, const int& k, const int& l, const int& m)
{
	bool *****R = new bool****[i];
	int n;
	for (n=0; n<i; n++) R[n] = CreateBool(j, k, l, m);
	return R;
}

inline void Destroy(double **v, const int &I)
{
	int i;
    //printf("Vor delete v[i]\n");
	for (i=0; i<I; i++) {
	    //printf("i=%d\n", i);
		delete[] v[i];
	}
    //printf("Vor delete v\n");
	delete[] v;
}

inline void Destroy(int **v, const int &I)
{
	int i;
	for (i=0; i<I; i++) delete[] v[i];
	delete[] v;
}

inline void Destroy(bool **v, const int &I)
{
	int i;
	for (i=0; i<I; i++) delete[] v[i];
	delete[] v;
}

inline void Destroy1(double **v, const int &I)
{
	int i;
	for (i=1; i<=I; i++) delete[] (v[i] + 1);
	delete[] (v + 1);
}

inline void Destroy(double ***v, const int &I, const int &J)
{
	int i;
	for (i=0; i<I; i++) Destroy (v[i], J);
	delete[] v;
}

inline void Destroy(bool ***v, const int &I, const int &J)
{
	int i;
	for (i=0; i<I; i++) Destroy(v[i], J);
	delete[] v;
}

inline void Destroy(int ***v, const int &I, const int &J)
{
	int i, j;
	for (i=0; i<I; i++) 
	{
		for (j=0; j<J; j++) delete[] v[i][j];
		delete[] v[i];
	}
	delete[] v;
}

inline void Destroy(double ****v, const int &I, const int &J, const int &K)
{
	int i;
	for (i=0; i<I; i++) Destroy(v[i], J, K);
	delete[] v;
}

inline void Destroy(bool ****v, const int &I, const int &J, const int &K)
{
	int i;
	for (i=0; i<I; i++) Destroy(v[i], J, K);
	delete[] v;
}

inline void Destroy(double***** v, const int& i, const int& j, const int& k, const int& l)
{
	int n;
	for (n=0; n<i; n++) Destroy(v[n], j, k, l);
	delete[] v;
}

inline void Destroy(bool***** v, const int& i, const int& j, const int& k, const int& l)
{
	int n;
	for (n=0; n<i; n++) Destroy(v[n], j, k, l);
	delete[] v;
}

inline int rund(double V)
{
	int R = int(V);
	double r = V - double(R);
	if (r >= 0.5) return ++R;
	if (r <= -0.5) return --R;
	return R;
}


#endif //UTILS
