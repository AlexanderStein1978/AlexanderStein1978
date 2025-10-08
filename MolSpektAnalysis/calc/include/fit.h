//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//

#ifndef FIT_H
#define FIT_H
/**
	@author Alexander Stein <AlexanderStein@t-online.de>
*/

#include "../../../Numerical-Recipes/src/include/nr3.h"

class Transition;
class MLRPot;
class IsoTab;
class ResidualFit;
struct SplinePoint;
class Spline;
class LocalPerturbation;


void InitDebugLogStream();
void DestroyDebugLogStream();

void FitBandConstants(double *o_BandConstants, double *o_err, int i_NBandC, double i_IsoF, double i_Omega, double *i_J, double *i_E, double *i_unc, int i_NData);
void AutoFitLocalPerturbation(LocalPerturbation& io_perturbation, const int i_maxBandC);
void FitLocalPerturbation(LocalPerturbation& io_perturbation, double i_pertBe);
void LocalPerturbationFitFuncs(const Doub i_J, VecDoub_I& i_par, Doub& o_E_calc, VecDoub_O& o_dE_dParam);
void FitToStraitLine(int N, double *x, double *y, double *unc, double x1, double x2,
					 double &y1, double &y2);
VecDoub PotFuncs(const Doub x);
double fitPotCoeff(int N, double *R, double *U, int n, double nRm, double nb, 
				   double *PotCoeff, double &Tm, double *Sig);
				
//VecDoub DunFitFuncs(const VecDoub_I &C);
double fitDun(int N, int nC, int **DunList, double ***vF, double ***JF, double *Res, double *err, int xx[][5],
			  double *yy, double *ssig, double *SS, double *dev, IsoTab *Iso, double *V = 0);
double **cEQS(int N, int &nC, int **DunList, double ***vF, double ***JF, int xx[][5], double *yy, IsoTab *Iso,
			  double *V = 0);
double cFQS(int N, int nC, double **EQS, int xx[][5], double *SS, double *x, double *ssig, double *dev);

void SvdFit(double **EQS, double *x, double *y, double *sig, int nC, int nL, double t);

bool fitMLRPot(double Rmin, double Rmax, double *U, int N, MLRPot *Pot, bool *fLRC);
void MLRFuncs(const Doub r, VecDoub_I &a, Doub &y, VecDoub_O &dyda);

void CalcLMatrix(double **L, SplinePoint *points, double iExp, int N, bool nat = false);
void SplineFuncs(const Doub x, VecDoub_I &a, Doub &y, VecDoub_O &dyda);
/*double fitSplinePot(Potential *Pot, double ****E, double *****K, int *mJ, int **mv, int Mv, int N, int NL, 
				    int **QN, double *SSig, double **EQS, double *EO, double *DE, 
				    double &minR, double &maxR, int NPar, double *Par);*/

void fitBeta(double *PotData, double Rmin, double Rmax, int NP, double EMin, double *beta, int Nb, int p, int q, double De, double Re,
			 double Rref, double *LRC, int *PLRC, int NLRC);
VecDoub betaFuncs(const Doub R);

void calcEK(double ****E, double *****K, double *U, double el_S, IsoTab *IsoT, double *IsoF, int NumFreePar,
			int *SMap, double **S, int NAdCorr, int TAdCorr, int PAdCorr, double RIso1AdCorr, double RIso2AdCorr, double AdCorr_Rm,
			double AdCorr_b, double *adCorr, int NumFreeAdCorr, int Omega, double q_e, double q_f, double SpinRRm, double SpinRb, 
			int NSpinRGamma, double *SpinRGamma, double SpinRR1, double SpinRR2, int *mJ, int **mv, int Mv, int N, int NL, 
            int **QN, double *Sig, double **EQS, double *EO, double *DE, double &minR, double &maxR, int NumWFPoints, double ****MaxErr = 0);
void calcAdCorrPot(double* InPot, double* OutPot, double Rmin, double Step, int NPoints, double MIso1, double MIso2,
	               int NAdCorr, int TAdCorr, int PAdCorr, double RIso1AdCorr, double RIso2AdCorr, double AdCorr_Rm, 
				   double AdCorr_b, double *adCorr);
void cRotPot(int J, int Iso, double *U, double *R, double minR, double h, int n, int Omega, double q_e, double q_f,
	         double *IsoF, double SpinRRm, double SpinRb, int NSpinRGamma, double *SpinRGamma, double SpinRR1, double SpinRR2,
			 int c = 0, int FC = 0);
void NumerovCooley(double *&Y, double M, int &nv, double *E, double **WF, 
                              int *ri, int *ra, int *CP, double *F, double *SQ, const int n, double *U, double *MinErr = 0);
double WFinoutMaxRatio(double* WF, int bi, int bo);
void getMinMaxInt(double *U, double SE, double M, const int n, int &xi, int &xo);
void Numerov(double M, double E, double &DE, int &v, const int n, double *U,
						double *WF, int &CP, double &F, double &SQ, double *Y);
void scatWaveFunc(double* U, double E, double M, int pi, int NP, double* WF);
void CalcLMEQS(int NEQ, double** EQS, double* Sig, double* EDiff, double** LMEQS, double hCi_cP, int NumFreePar);
void solveLinEQSbyLU(double ** EQS, int N);
void solveLinEQSbySVD(double** EQS, int N, double t);
double calcEQSFQS(MatDoub a, VecDoub b, VecDoub x);

void calcFFT(const double *const inputData, const int N, const double delta, double **const realOutput, double **const imaginaryOutput);
void backtransformFFT(const double *const realInput, const double *const imaginaryInput, const int N, const double delta, double **const realOutput, double **const imaginaryOutput);

inline double SpinRotation(double R, double SpinRRm, double SpinRb, int NSpinRGamma, double *SpinRGamma, 
						   double SpinRR1, double SpinRR2)
{
	double x = (R - SpinRRm) / (R + SpinRb * SpinRRm), E = 0.0;
	int n;
	for (n = NSpinRGamma - 1; n>=0; n--)
	{
		E *= x;
		E += SpinRGamma[n];
	}
	E /= (exp((R - SpinRR1) / SpinRR2) + 1.0);
	return E;
}


class DAddMaximizeFunctor
{
public:
	DAddMaximizeFunctor(Spline *spline, double **Data, int NJ, 
						int JStart, int JStep, int NE, double EStart, double ERes, bool absorption);
	Doub operator()(VecDoub x);
	
private:
	Spline *spline;
	int NJ, JStart, JStep, NE;
	double **Data, EStart, ERes;
	bool isAbsorption;
};

void DAddMaximize(ResidualFit *resFit, double **Data, int NJ, int JStart, int JStep,
                  int NE, double EStart, double ERes, double Diff, bool absorption);

#endif
