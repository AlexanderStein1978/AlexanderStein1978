//
// C++ Interface: potential
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2007 - 2019
//
// Copyright: See README file that comes with this source code
//
//


#ifndef TANGTOENNIESPOT_H
#define TANGTOENNIESPOT_H


class TangToenniesPot : public PotWorker
{
public:
	TangToenniesPot(PotFit *Fit);
	TangToenniesPot(const TangToenniesPot &C);
	~TangToenniesPot();
	void getMinimum(double &R, double &T);
	void setPotCoeff(double *A, int NA, double b, double *C, int N, double Uinf);
	void getPotCoeff(double *&A, int &NA, double &b, double *&C, int &N, double &Uinf);
	bool calcPotCoeff(double *CC, int NCC, int FNC, double Re, double De);
	double calcPotCoeff(int numPoints, double *R, double *U, double *Sig, int numCoeff, double b);
	double improvePotential(int numPoints, double *R, double *U, double *Sig, double nb = 0.0);
	void updateCoefficients(double *CD);
	PotentialData *getPotentialData();
	
	inline double GetUinf()
	{
		return Uinf;
	}
	
	inline double getb()
	{
		return b;
	}
	
	inline void setb(double nb)
	{
		b = nb;
	}
	
	inline int getNFitCoefficients()
	{
		return NC + NA + (A_ex == 0.0 ? 0 : 1);
	}
	
	inline void setExcInt(double A, double Beta, double Gamma)
	{
		A_ex = A;
		beta = Beta;
		gamma = Gamma;
	}
	
	inline void getExcInt(double &A, double &Beta, double &Gamma)
	{
		A = A_ex;
		Beta = beta;
		Gamma = gamma;
	}
	
	inline void setNA(int nA)
	{
		int n;
		if (A != 0) delete[] A;
		if (nA > 0) for (n=0, A = new double[NA = nA]; n < nA; n++) A[n] = 0.0;
		else A=0;
	}
	
protected:
	void getLRCoeff(double &R, int &numCoefficients, int *&Exponents, double *&Coefficients);
	void calcS(int numPoints, double Rmin, double Rmax);
	double Point(double R, int FC = 0);
	
private:
	void getPreFakt(double R, double *C);
	bool calcAb(double Re, double De);
	void AB(double BI, double &Y, double &AA, double *D, double RM);
	
	double Uinf, *A, b, *C, *SWF, A_ex, beta, gamma, RMin, UMax;
	int NC, NA;
};

#endif
