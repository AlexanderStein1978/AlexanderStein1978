//
// C++ Implementation: fit
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2008 - 2020
//
// Copyright: See README file that comes with this source code
//
//


#define UseDebugLogStream

#ifdef UseDebugLogStream
#include <QFile>
#include <QTextStream>
#include <QString>

QFile* DebugFile;
QTextStream* DebugStream;
#endif

#include "fit.h"
#include "utils.h"
#include "constants.h"
#include "isotab.h"
#include "SplinePoint.h"
#include "Spline.h"
#include "ResidualFit.h"
#include "LocalPerturbation.h"

#include "../../NR_C301/code/gamma.h"
#include "../../NR_C301/code/incgammabeta.h"
#include "../../NR_C301/code/fitab.h"
#include "../../NR_C301/code/svd.h"
#include "../../NR_C301/code/fitsvd.h"
#include "../../NR_C301/code/ludcmp.h"
#include "../../NR_C301/code/gaussj.h"
#include "../../NR_C301/code/fitmrq.h"
#include "../../NR_C301/code/amoeba.h"
#include "../../NR_C301/code/fourier.h"

#include <limits>
#include <cmath>

int nPotCoeff, NBeta, p, q;
double Rref, Re;
Doub Rm, b;
MLRPot *mlrPot;
//Potential *Pot;
double ****E, *****K, *Par, *SSig, **EQS, *EO, *DE, minR, maxR;
int **QN, NPar, *mJ, **mv, Mv, N, NL;


void InitDebugLogStream()
{
#ifdef UseDebugLogStream
	DebugFile = new QFile("DebugLog.txt");
	DebugFile->open(QIODevice::WriteOnly);
	DebugStream = new QTextStream(DebugFile);
#endif
}

void DestroyDebugLogStream()
{
#ifdef UseDebugLogStream
	delete DebugStream;
	DebugStream = 0;
	delete DebugFile;
	DebugFile = 0;
#endif
}

void FitToStraitLine(int N, double* x, double* y, double* unc, double x1, double x2, double& y1, double& y2)
{
	VecDoub xx(N), yy(N), ssig(N);
	int n;
	for (n=0; n<N; n++)
	{
		xx[n] = x[n];
		yy[n] = y[n];
		ssig[n] = unc[n];
	}
	Fitab Fit(xx, yy, ssig);
	y1 = Fit.a + x1 * Fit.b;
	y2 = Fit.a + x2 * Fit.b;
}

VecDoub PotFuncs(const Doub R)
{
	VecDoub C(nPotCoeff+1);
	int i;
	Doub x((R - Rm) / (R + b * Rm));
	for (i=1, C[0] = 1.0; i <= nPotCoeff; i++) C[i] = x * C[i-1];
	return C;
}

double fitPotCoeff(int N, double *R, double *U, int n, double nRm, double nb, 
				   double *PotCoeff, double &Tm, double *Sig)
{
	b = nb;
	Rm = nRm;
	nPotCoeff = n;
	VecDoub xx(N), yy(N), ssig(N);
	int i;
	for (i=0; i<N; i++)
	{
		xx[i] = R[i];
		yy[i] = U[i];
		ssig[i] = Sig[i];
	}
	Fitsvd SVD(xx, yy, ssig, PotFuncs);
	SVD.ma = n + 1;
	SVD.fit();
	Tm = SVD.a[0];
	for (i=0; i<n; i++) PotCoeff[i] = SVD.a[i+1];
	printf("FQS=%g\n", SVD.chisq);
	return SVD.chisq;
}

/*VecDoub DunFitFuncs(const VecDoub_I &x)
{
	//printf("Beginn DunFitFuncs, nP=%d, nC=%d, x[0]=%f, x[1]=%f, x[2]=%f\n", nP, nC, x[0], x[1], x[2]);
	VecDoub R(nC);
	int n;
	for (n=0; n < nC; n++) 
		R[n] = vF[int(x[0])][int(x[1])][DunList[n][0]] * JF[int(x[0])][int(x[2])][DunList[n][1]];
	//for (; n < nC + nP; n++) R[n] = 0.0;
	//printf("nP=%d, x[3]=%f\n", nP, x[3]);
	//if (x[3] >= 0.0) R[nC + int(x[3])] = 1.0;
	//printf("Ende DunFitFuncs\n");
	//for (n=0; n < nC; n++) printf("R[%d]=%f\n", n, R[n]);
	return R;
}*/

double **cEQS(int N, int &nC, int **DunList, double ***vF, double ***JF, int xx[][5], double *yy, IsoTab *Iso, double *V)
{
	int i, j, n, NC = nC;
	if (V!=0) for (nC = i = 0; i < NC; i++) if (V[i] == 0.0) nC++;
		//printf("V[%d]=%f\n", i, V[i]);}
	double **EQS = Create(N, nC + 1);
	//for (i=0; i < NC; i++) 
		//printf("NC=%d, DunList[%d][0]=%d, DunList[%d][1]=%d\n", NC, i, DunList[i][0], i, DunList[i][1]);
	//printf("cEQS: N=%d\n", N);
	for (n=0; n<N; n++)
	{
		//printf("xx[%d][0]=%d, xx[%d][1]=%d, xx[%d][2]=%d, xx[%d][3]=%d\n", 
			//   n, xx[n][0], n, xx[n][1], n, xx[n][2], n, xx[n][3]);
		for (i=j=0, EQS[n][nC] = 0.0; i < NC; i++)
		{
			switch (DunList[i][2])
			{
				case 0:
					if (V!=0 ? V[i]!=0.0 : false) 
						EQS[n][nC] 
								-= V[i] * vF[xx[n][0]][xx[n][1]][DunList[i][0]] * JF[xx[n][0]][xx[n][2]][DunList[i][1]];
					else
						EQS[n][j++] = vF[xx[n][0]][xx[n][1]][DunList[i][0]] * JF[xx[n][0]][xx[n][2]][DunList[i][1]];
					break;
				case 1:
					if (2 * (xx[n][4] / 2) == xx[n][4])
					{
						if (V!= 0 ? V[i] != 0.0 : false)
							EQS[n][nC] -= V[i] * vF[xx[n][0]][xx[n][1]][DunList[i][0]] * JF[xx[n][0]][xx[n][2]][DunList[i][1]];
						else EQS[n][j++] = vF[xx[n][0]][xx[n][1]][DunList[i][0]] * JF[xx[n][0]][xx[n][2]][DunList[i][1]];
					}
					else if (V!=0 ? V[i] == 0 : true) EQS[n][j++] = 0.0;
					break;
				case 2:
					if (xx[n][4] >= 2)
					{
						//printf("n=%d, F1\n", n);
						if (V!= 0 ? V[i] != 0.0 : false)
							EQS[n][nC] += V[i] * vF[xx[n][0]][xx[n][1]][DunList[i][0]] 
											* (xx[n][2] + 1) * JF[xx[n][0]][xx[n][2]][DunList[i][1]-1];
						else EQS[n][j++] = -vF[xx[n][0]][xx[n][1]][DunList[i][0]] 
											* (xx[n][2] + 1) * JF[xx[n][0]][xx[n][2]][DunList[i][1]-1];
					}
					else
					{
						//printf("n=%d, F0\n", n);
						if (V!= 0 ? V[i] != 0.0 : false)
							EQS[n][nC] -= V[i] * vF[xx[n][0]][xx[n][1]][DunList[i][0]] 
											* xx[n][2] * JF[xx[n][0]][xx[n][2]][DunList[i][1]-1];
						else 
							EQS[n][j++] = vF[xx[n][0]][xx[n][1]][DunList[i][0]] 
											* xx[n][2] * JF[xx[n][0]][xx[n][2]][DunList[i][1]-1];
					}
					break;
				case 3:
					if (V!=0 ? V[i] != 0.0 : false)
						EQS[n][nC] -= V[i] * (1.0 - Iso->relRedMass[xx[n][0]]) 
										* vF[xx[n][0]][xx[n][1]][DunList[i][0]] * JF[xx[n][0]][xx[n][2]][DunList[i][1]];
					else
						EQS[n][j++] = (1.0 - Iso->relRedMass[xx[n][0]]) 
										* vF[xx[n][0]][xx[n][1]][DunList[i][0]] * JF[xx[n][0]][xx[n][2]][DunList[i][1]];
					break;
			}
			//printf(" EQS[%d][%d]=%f,", n, i, EQS[n][i]);
		}
		EQS[n][nC] += (xx[n][3]==-1 ? yy[n] : -yy[n]);
		//printf(" EQS[%d][nC=%d]=%f\n", n, nC, EQS[n][nC]);
	}
	//for (i=0; i<nC; i++) 
		//printf("DunList[%d][0]=%d, DunList[%d][1]=%d\n", i, DunList[i][0], i, DunList[i][1]);
	return EQS;
}

double cFQS(int N, int nC, double **EQS, int xx[][5], double *SS, double *x, double *ssig, 
			double *dev)
{
	int i, l, n, m;
	double B = 0.0, FQS, T[N], Max, stDev = 0.0;
	for (l=-1, n=m=0; n<N; n++) 
	{
		for (T[n]=0.0, i=0; i < nC; i++) T[n] += EQS[n][i] * x[i];
		if (xx[n][3] != l)
		{
			if (l != -1)
			{
				//printf("B=%g ", B);
				B /= SS[l];
				//printf("SS[%d]=%g, B=%g\n", l, SS[l], B);
				for (i=m; i<n; i++) 
				{
					//printf("EQS[%d][nC]=%g, yy[%d]=%g, B-yy[%d]=%g\n", 
						//   i, EQS[i][nC], i, yy[i], i, B-yy[i]);
					EQS[i][nC] += B;
					//if (!(EQS[i][nC] - B < 1.0)) 
						//printf("error: i=%d, l=%d, n=%d, m=%d\n", i, l, n, m);
				}
			}
			l = xx[n][3];
			B = 0.0;
			m = n;
		}
		if (xx[n][3] != -1) B += (T[n] - EQS[n][nC]) / ssig[n];
	}
	if (l!=-1) 
	{
		B /= SS[l];
		for (i=m; i<n; i++) EQS[i][nC] += B; 
	}
	//for (n=0; n<nC; n++) for (i=0, b[n] = 0.0; i<N; i++) b[n] += EQS[i][n] * EQS[i][nC] * ssig[i];
		//printf("b[%d]=%g\n", n, b[n]);}
	for (FQS = 0.0, Max = 0.0, n=0; n<N; n++) 
	{
		dev[n] = B =  EQS[n][nC] - T[n];
		stDev += (B*=B);
		FQS += B / ssig[n];
		if (fabs(dev[n]) > fabs(Max)) Max = dev[n]; 
	}
    //printf("Max dev=%f, stdDev=%f, DevSqS=%f, N=%d\n", Max, sqrt(stDev / (N - nC)), stDev, N);
	//printf("FQS = %f\n", FQS);
	return FQS;
}

double fitDun(int N, int nC, int **DunList, double ***vF, double ***JF, double *Res, double *err, 
			  int xx[][5], double *yy, double *ssig, double *SS, double *dev, IsoTab *Iso, double *V)
{
	//printf("Beginn fitDun\n");
	int i, j, r, s, n, l, m;
	double **EQS = cEQS(N, nC, DunList, vF, JF, xx, yy, Iso, V), sig;
	MatDoub a(nC, nC), C(nC, nC);
	VecDoub b(nC), x(nC);
	/*for (n=0; n<nC; n++) 
	{
		for (m=0; m<nC; m++) for (i=0, a[n][m] = 0.0; i<N; i++) 
				a[n][m] += EQS[i][n] * EQS[i][m] * ssig[i];
			
		for (i=0, b[n] = 0.0; i<N; i++) b[n] += EQS[i][n] * EQS[i][nC] * ssig[i];
	
	}*/
	for (n=0; n<nC; n++)
	{
		b[n] = 0.0;
		for (m=0; m<nC; m++) a[n][m] = 0.0;
	}
	for (l=-1, n=m=0; n<=N; n++)
	{
		if (n<N ? xx[n][3] != l : true)
		{
			if (l != -1) for (i=m+1; i<n; i++) for (j=m; j<i; j++) 
						for (r=0, sig = 1.0 / (ssig[i] + ssig[j]); r<nC; r++)
			{
				for (s=0; s<nC; s++) a[s][r] += (EQS[i][s] - EQS[j][s]) * (EQS[i][r] - EQS[j][r]) * sig;
				b[r] += (EQS[i][r] - EQS[j][r]) * (EQS[i][nC] - EQS[j][nC]) * sig;
			}
			if (n<N) 
			{
				l = xx[n][3];
				m=n;
			}
		}
		if (n<N ? xx[n][3] == -1 : false) for (r=0, sig = 1.0 / ssig[n]; r<nC; r++)
		{
			for (s=0; s<nC; s++) a[s][r] += EQS[n][s] * EQS[n][r] * sig;
			b[r] += EQS[n][r] * EQS[n][nC] * sig;
		}
	}	
	LUdcmp LU(a);
	double FQS = -1.0, aFQS = -1.0;
	while (aFQS == -1.0 || (aFQS - FQS) / FQS > 1e-7)
	{
		aFQS = FQS;
		if (FQS == -1.0) LU.solve(b, x);
		else LU.mprove(b, x);
		for (n=0; n<nC; n++) Res[n] = x[n];
		FQS = cFQS(N, nC, EQS, xx, SS, Res, ssig, dev);
	}
	
	//for (n=0; n<nC; n++) printf("P[%d]=%g\n", n, x[n]);
	LU.inverse(C);
	for (n=0; n < nC; n++) err[n] = sqrt(C[n][n]);
	/*for (n=0; n < nC; n++) 
	{
		Res[n] = FS.a[n];
		err[n] = sqrt(FS.covar[n][n]);
	}
	for (n=0; n < nC; n++) printf("P%d=%g\n", n, FS.a[n]);
	printf("FQS=%g\n", FS.chisq);*/
	Destroy(EQS, N);
	return FQS;
}

void solveLinEQSbyLU(double** EQS, int N)
{
	int n, m;
	double FQS, lFQS;
	MatDoub a(N, N);
	VecDoub b(N), x(N);
	for (n=0; n<N; n++)
	{
		for (m=0; m<N; m++) a[n][m] = EQS[n][m];
		b[n] = EQS[n][N];
	}
	LUdcmp LU(a);
	LU.solve(b, x);
	lFQS = calcEQSFQS(a, b, x);
	LU.mprove(b, x);
	FQS = calcEQSFQS(a, b, x);
	while ((lFQS - FQS) / FQS > 1e-7)
	{
		lFQS = FQS;
		LU.mprove(b, x);
		FQS = calcEQSFQS(a, b, x);
	}
	for (n=0; n<N; n++) EQS[n][N] = x[n];
}

double calcEQSFQS(MatDoub a, VecDoub b, VecDoub x)
{
	int n, m, N = b.size();
	double FQS = 0.0, S;
	for (n=0; n < N; n++)
	{
		for (m=0, S = 0.0; m<N; m++) S += a[n][m] * x[m];
		S -= b[n];
		FQS += S*S;
	}
	printf("FQS=%g\n", FQS);
	return FQS;
}

void solveLinEQSbySVD(double** EQS, int N, double t)
{
	int n, m;
	MatDoub a(N, N);
	VecDoub b(N), x(N);
	for (n=0; n<N; n++)
	{
		for (m=0; m<N; m++) a[n][m] = EQS[n][m];
		b[n] = EQS[n][N];
	}
	SVD svd(a);
	svd.solve(b, x, t);
	calcEQSFQS(a, b, x);
	for (n=0; n<N; n++) EQS[n][N] = x[n];
}

void SvdFit(double **EQS, double *x, double *y, double *sig, int nC, int nL, double t)
{
	MatDoub aa(nL, nC);
	VecDoub b(nL), a(nC);
	int n, c;
	//double KS, W;
	//printf("ssig[0]=%f\n", sig[0]);
	for (n=0; n<nL; n++) 
	{
		//for (c=0, KS = 0.0; c < nC; c++) KS += fabs(EQS[n][c]);
		for (c=0; c < nC; c++) aa[n][c] = sig[n] * EQS[n][c];
		b[n] = sig[n] * y[n];
	}
	SVD svd(aa);
	svd.solve(b, a, t);
	for (c=0; c < nC; c++) 
	{
		x[c] = a[c];
		//printf("x[%d]=%g\n", c, x[c]);
	}
}

/*void MLRFuncs(const Doub r, VecDoub_I& a, Doub& E, VecDoub_O& dyda)
{
	int p, q, Nbeta, NLRC, *pLRC, n, m;
	double TAs, De, Re, RRef, *beta, *LRC;
	mlrPot->getCoefficients(p, q, TAs, De, Re, RRef, Nbeta, beta, NLRC, pLRC, LRC);
	TAs = a[0];
	De = a[1];
	Re = a[2];
	for (n=0, m=3; n < Nbeta; n++, m++) beta[n] = a[m];
	for (n=0; n < NLRC; n++, m++) LRC[n] = a[m];
	mlrPot->setCoefficients(p, q, TAs, De, Re, RRef, Nbeta, beta, NLRC, pLRC, LRC);
	double dbeta[Nbeta], dC[NLRC];
	mlrPot->getDerivatives(r, E, dbeta, dC, Re, De);
	dyda[0] = TAs;
	dyda[1] = De;
	dyda[2] = Re;
	for (n=0, m=3; n < Nbeta; n++, m++) dyda[m] = dbeta[n];
	for (n=0; n < NLRC; n++, m++) dyda[m] = dC[n];
}*/

/*bool fitMLRPot(double Rmin, double Rmax, double *U, int N, MLRPot* Pot, bool *fLRC)
{
	int p, q, Nbeta, NLRC, *pLRC, n, m;
	double cTAs, TAs, De, Re, RRef, *beta, *LRC, r, h = (Rmax - Rmin) / double(N - 1);
	bool R = true;
	mlrPot = Pot;
	mlrPot->getCoefficients(p, q, TAs, De, Re, RRef, Nbeta, beta, NLRC, pLRC, LRC);
	VecDoub xx(N), yy(N), ssig(N), aa(Nbeta + NLRC + 3);
	for (n=0, r = Rmin; n<N; n++, r+=h)
	{
		xx[n] = r;
		yy[n] = U[n];
		ssig[n] = 1e-3;
	}
	aa[0] = cTAs = TAs;
	aa[1] = De;
	aa[2] = Re;
	for (n=0, m=3; n < Nbeta; n++) aa[m] = beta[n];
	for (n=0; n < NLRC; n++, m++) aa[m] = LRC[n];
	Fitmrq FLM(xx, yy, ssig, aa, MLRFuncs, 1e-3);
	if (TAs > 0.0) FLM.hold(0, TAs);
	for (n=0, m = Nbeta + 3; n < NLRC; n++, m++) if (fLRC[n]) FLM.hold(m, LRC[n]);
	try
	{
		FLM.fit();
	}
	catch (...)
	{
		//printf(err);
		R = false;
	}
	delete[] beta;
	delete[] pLRC;
	delete[] LRC;
	return R;
}*/

VecDoub betaFuncs(const Doub R)
{
	double Rp = pow(R, p), Rep = pow(Re, p), Rrefp = pow(Rref, p), Rq = pow(R, q), Rrefq = pow(Rref, q);
	double yrefq = (Rq - Rrefq) / (Rq + Rrefq);
	int n;
	VecDoub Rv(NBeta);
	Rv[0] = (Rep - Rp) * 2.0 * Rrefp / ((Rp + Rep) * (Rp + Rrefp));
	for (n=1; n < NBeta; n++) Rv[n] = Rv[n-1] * yrefq;
	return Rv;
}

void fitBeta(double *Data, double Rmin, double Rmax, int NP, double Offs, double* beta, int Nb, int np, int nq, double De, double nRe, 
			 double nRref, double* LRC, int* PLRC, int NLRC)
{
	int n, m;
	double ULR, ULRRe = 0.0, betaInf, R, h = (Rmax - Rmin) / double(NP - 1);// a = sqrt(mu * C_u * 2e-18 * C_c/ (C_h * De)) * we * M_PI;
	double Rep = pow(Re = nRe, p = np), Rp, dDe = 1.0 / De;
	double Rrefp = pow(Rref = nRref, p);
	VecDoub xx(NP), yy(NP), sig(NP);
	NBeta = Nb;
	q = nq;
	for (m=0; m < NLRC; m++) ULRRe += pow(Re, -PLRC[m]) * LRC[m];
	betaInf = log(2.0 * De / ULRRe);
	for (R = Rmin, n=0; n < NP; n++, R+=h)
	{
		xx[n] = R;
		for (ULR = 0.0, m=0; m < NLRC; m++) ULR += pow(R, -PLRC[m]) * LRC[m];
		Rp = pow(R, p);
		//yy[n] = a * (Re - R) - log(ULR * ULRRe) + (Rp - Rep) * (Rp - Rrefp) / ((Rp + Rep) * (Rp + Rrefp)) * betaInf;
		if (R < Re) yy[n] = log((sqrt((Data[n] - Offs) * dDe) + 1.0) / ULR * ULRRe);
		else yy[n] = log((1.0 - sqrt((Data[n] - Offs) * dDe)) / ULR * ULRRe);
		yy[n] += (Rp - Rep) * (Rp - Rrefp) / ((Rp + Rep) * (Rp + Rrefp)) * betaInf;
		sig[n] = 0.01;
	}
	Fitsvd SVD(xx, yy, sig, betaFuncs);
	SVD.fit();
	for (n=0; n < NBeta; n++) beta[n] = SVD.a[n];
	printf("Betafit: chisq=%e\n", SVD.chisq);
	/*VecDoub bF;
	for (R = Rmin, n=0; R <= Rmax; R += 100 * h, n+=100)
	{
		bF = betaFuncs(R);
		for (ULR = 0.0, m=0; m < NLRC; m++) ULR += pow(R, -PLRC[m]) * LRC[m];
		for (Rp = 0.0, m=0; m < Nb; m++) Rp += bF[m] * SVD.a[m];
		printf ("R=%f, a=%e, b=%e, exp=%e\n", R, a * (Re - R), Rp - yy[n] + a * (Re - R), Rp - yy[n] + a * (Re - R) - log(ULR * ULRRe));
	}*/
}

/*void SplineFuncs(const Doub x, VecDoub_I& a, Doub& y, VecDoub_O& dyda)
{
	int n, i = int(x);
	for (n=0; (n < NPar ? a[n] == Par[n] : false); n++) ;
	if (n < NPar)
	{
		for (n=0; n < NPar; n++) Par[n] = a[n];
		Pot->SetSplinePar(Par, E, K, mJ, mv, Mv, N, NL, QN, SSig, EQS, EO, DE, minR, maxR);
	}
	y = E[QN[i][5]][QN[i][0]][QN[i][1]][QN[i][2]];
	for (n=0; n < NPar; n++) dyda[n] = K[QN[i][5]][QN[i][0]][QN[i][1]][QN[i][2]][n];
}*/

/*double fitSplinePot(Potential* nPot, double**** nE, double***** nK, int* nmJ, int** nmv, int nMv, int nN, int nNL,
				    int** nQN, double* nSSig, double** nEQS, double* nEO, double* nDE, 
				    double& nminR, double& nmaxR, int nNPar, double* nPar)
{
	int n, l, m, i, done;
	double wt, ochisq, chisq = 0.0, alamda = 1e-3;
	Pot = nPot;
	E = nE;
	K = nK;
	mJ = nmJ;
	mv = nmv;
	Mv = nMv;
	N = nN;
	NL = nNL;
	QN = nQN;
	SSig = nSSig;
	EQS = nEQS;
	EO = nEO;
	DE = nDE;
	minR = nminR;
	maxR = nmaxR;
	NPar = nNPar;
	Par = nPar;
	VecDoub /x(NL), y(NL),/ sig2i(NL), a(NPar), beta(NPar), da(NPar);
	MatDoub alpha(NPar, NPar), covar(NPar, NPar);
	for (n=0; n < NPar; n++) for (l=0, beta[n] = 0.0; l < NPar; l++) alpha[n][l] = 0.0;
	for (n=0; n < NL; n++)
	{
		//x[n] = n;
		//y[n] = EO[n];
		sig2i[n] = SSig[n] * SSig[n];
		for (l=0; l < NPar; l++)
		{
			wt = EQS[n][l] * sig2i[n];
			for (m=0; m<=l; m++) alpha[l][m] += wt * (EQS[n][m]);
			beta[l] += DE[n] * wt;
		}
		chisq += DE[n] * DE[n] * sig2i[n];
	}
	ochisq = chisq;
	for (l=1; l < NPar; l++) for (m=0; m<l; m++) alpha[m][l] = alpha[l][m];
	for (n=0; n < NPar; n++) a[n] = Par[n];
	//Fitmrq MRQ(xx, yy, ssig, aa, SplineFuncs);
	//MRQ.fit();
	while (true)
	{
		if (done == 4) alamda = 0.0;
		for (l=0; l < NPar; l++) 
		{
			for (m=0; m < NPar; m++) covar[l][m] = alpha[l][m];
			covar[l][l] *= (1.0 + alamda);
			for (m=0; m < NPar; m++) printf("covar[%d][%d]=%f\t", l, m, covar[l][m]);
		}
		SVD svd(covar);
		svd.solve(beta, da, 1e-7);
		if (done == 4) break;
		for (l=0, chisq = 0.0; l < NPar; l++) 
		{
			Par[l] = a[l] + da[l];
			for (m=0, da[l] = 0.0; m < NPar; m++) covar[l][m] = 0.0;
		}
		Pot->SetSplinePar(Par, E, K, mJ, mv, Mv, N, NL, QN, SSig, EQS, EO, DE, minR, maxR);
		for (n=0; n < NL; n++)
		{
			for (l=0; l < NPar; l++)
			{
				wt = EQS[n][l] * sig2i[n];
				for (m=0; m<=l; m++) covar[l][m] += wt * EQS[n][m];
				da[l] += DE[n] * wt;
			}
			chisq += DE[n] * DE[n] * sig2i[n];
		}
		for (l=1; l < NPar; l++) for (m=0; m<l; m++) covar[m][l] = covar[l][m];
		if (abs(chisq - ochisq) < MAX(1e-8, 1e-8 * chisq)) done++;
		if (chisq < ochisq)
		{
			alamda *= 0.1;
			ochisq = chisq;
			printf("chisq=%g\n", chisq);
			for (l=0; l < NPar; l++) for (m=0, beta[l] = da[l], a[l] = Par[l]; m < NPar; m++) 
				alpha[l][m] = covar[l][m];
		}
		else
		{
			alamda *= 10.0;
			chisq = ochisq;
		}
	}
	for (n=0; n < NPar; n++) Par[n] = a[n];
	return chisq;
}*/

void CalcLMatrix(double** L, SplinePoint* points, double iExp, int N, bool nat)
{
	int n, m, s = (nat ? 1 : 0);
	double d3 = 1.0 / 3, d6 = 1.0 / 6, dx[N-1];
	MatDoub Li(N, N+2-2*s), a(N, N);
	MatDoub Lo(N, N+2-2*s);
	for (n=0; n<N; n++)
	{
		for (m=0; m<N; m++) a[n][m] = 0.0;
		for (m=0; m<=N+1-2*s; m++) Li[n][m] = 0.0;
	}
	for (n=0; n<N-1; n++) dx[n] = points[n+1].x - points[n].x;
	Li[0][1] = -1.0 / dx[0];
	Li[0][2] = -Li[0][1];
	if (nat) Li[0][1] = Li[0][2];
	a[0][0] = d3 * dx[0];
	a[0][1] = 0.5 * a[0][0];
	for (n=1; n<N-1; n++)
	{
		a[n][n-1] = a[n-1][n];
		a[n][n] = d3 * (points[n+1].x - points[n-1].x);
		a[n][n+1] = d6 * dx[n];
		Li[n][n-s] = Li[n-1][n+1-s];
		Li[n][n+2-s] = 1.0 / dx[n];
		Li[n][n+1-s] = -Li[n][n-s] - Li[n][n+2-s];
	}
	if (nat)
	{
		a[0][0] = a[N-1][N-1] = 1.0;
		Li[0][1] = Li[0][2] = a[0][1] = 0.0;
	}
	else
	{
		Li[0][0] = (iExp > 0.0 ? iExp * pow(points[0].x, -iExp - 1.0) : -1);
		Li[N-1][N-1] = 1.0 / dx[N-2];
		Li[N-1][N] = - Li[N-1][N-1];
		Li[N-1][N+1] = 1.0;
		a[N-1][N-2] = d6 * dx[N-2];
		a[N-1][N-1] = 2.0 * a[N-1][N-2];
	}
	SVD svd(a);
	svd.solve(Li, Lo, 1e-12);
	for (n=0; n<N; n++) for (m=0; m<=N+1-2*s; m++) L[n][m] = Lo[n][m];
}

void calcEK(double ****E, double *****K, double *U, double el_S, IsoTab *IsoT, double *IsoF, int NumFreePar,
			int *SMap, double **S, int NAdCorr, int TAdCorr, int PAdCorr, double RIso1AdCorr, double RIso2AdCorr, double AdCorr_Rm,
			double AdCorr_b, double *adCorr, int NumFreeAdCorr, int Omega, double q_e, double q_f, double SpinRRm, double SpinRb, 
			int NSpinRGamma, double *SpinRGamma, double SpinRR1, double SpinRR2, int *mJ, int **mv, int Mv, int N, int NL, 
            int **QN, double *Sig, double **EQS, double *EO, double *DE, double &minR, double &maxR, int NumWFPoints, double ****MaxErr)
{
	int I, J, v, n, m, k, nv, minn = NumWFPoints, maxn = 0, FC, NFC = int(2.0 * el_S) + 1, NIso = IsoT->numIso;
	double h = (rmax - rmin) / (NumWFPoints - 1), *RPot = new double[NumWFPoints], Sum, CE1/*, CE2 = 1.0*/, M, aM;
	double MF, *Y = new double[NumWFPoints];
	double **WF = Create(Mv + 1, NumWFPoints), *AdCorrPot = (NAdCorr > 0 ? new double[NumWFPoints]: U), uE = 0.0, ErrS;
	//printf("Nach Create WF\n");
	//aPot->GetMin(R, M);
	double *F = new double[Mv+1], MFAdCorr = 1.0;
	double *SQ = new double[Mv+1], uK[N], Sq;
	//double *U2, *A2, *R2, *E2;
	int *CP = new int[Mv+1];
	int *ri = new int[Mv+1];
	int *ra = new int[Mv+1];
	/*if (Type == SplinePotential && L==0) calcLMatrix();
	else if (Type != SplinePotential && AdCorrWF == 0) AdCorrWF = new double[NAdCorr];
	if (S==0) calcS(NumWFPoints, rmin, rmax);*/
	/*if (WFSource != 0) 
	{
		U2 = WFSource->getPoints(rmin, rmax, NumWFPoints);
		E2 = new double[Mv + 1];
		if (NAdCorr > 0) A2 = new double[NumWFPoints];
		else A2 = U2;
		R2 = new double[NumWFPoints];
	}*/
	//printf("Nach calcS\n");
	for (FC=0; FC < NFC; FC++) for (I=0; I < NIso; I++) 
	{
		//printf("I=%d\n", I);
		MF = IsoF[I] * h*h;
		if (NAdCorr > 0) 
		{
			calcAdCorrPot(U, AdCorrPot, rmin, h, NumWFPoints, IsoT->mIso1[I], IsoT->mIso2[I], NAdCorr, TAdCorr, PAdCorr, RIso1AdCorr,
				          RIso2AdCorr, AdCorr_Rm, AdCorr_b, adCorr);
			//if (WFSource != 0) 
				//WFSource->calcAdCorrPot(U2, A2, rmin, h, NumWFPoints, IsoT->mIso1[I], IsoT->mIso2[I]);
		}
		if (NumFreeAdCorr > 0) switch (TAdCorr)
		{
			case 1:
				MFAdCorr = 1.0 - RIso1AdCorr / IsoT->mIso1[I];
				break;
			case 2:
				MFAdCorr = 1.0 - RIso2AdCorr / IsoT->mIso2[I];
				break;
			case 3:
				MFAdCorr = 1.0 - (RIso1AdCorr * RIso2AdCorr) / ((RIso1AdCorr + RIso2AdCorr) * IsoT->redMass[I]);
				break;
		}
		for (J=0; J <= mJ[I]; J++) if (mv[I][J] > -1)
		{
			cRotPot(J, I, AdCorrPot, RPot, rmin, h, NumWFPoints, Omega, q_e, q_f, IsoF, SpinRRm, SpinRb,
				    NSpinRGamma, SpinRGamma, SpinRR1, SpinRR2, 0, (NFC > 1 ? FC + 1 : 0));
			NumerovCooley(Y, MF, nv = mv[I][J] + 1, E[FC][I][J], WF, ri, ra, CP, F, SQ, NumWFPoints, RPot, 
						  (MaxErr != 0 ? MaxErr[FC][I][J] : 0));
			/*if (WFSource != 0) 
			{
				WFSource->CRotPot(J, I, A2, R2, rmin, h, NumWFPoints, 0, (NFC > 1 ? FC + 1 : 0));
				WFSource->NumerovCooley(Y, MF, nv, E2, WF, ri, ra, CP, F, SQ, NumWFPoints, R2,
										(MaxErr != 0 ? MaxErr[FC][I][J] : 0));
			}*/
			if (nv < mv[I][J] + 1)
			{
				printf("For J=%d and Iso%d v=%d is not existing, highest v is %d\n", 
					   J, I, mv[I][J], nv - 1);
				for (v = nv; v <= mv[I][J]; v++) 
				{
					E[FC][I][J][v] = 0.0;
					for (k=0; k<N; k++) K[FC][I][J][v][k] = 0.0;
				}
			}
			if (nv == 0) continue;
			for (n = ri[v = nv - 1], M = 0.0, aM = -1.0; M > aM && n <= ra[v]; n++)
			{
				aM = M;
				M = fabs(WF[v][n]);
			}
			for (M *= 0.05; fabs(WF[v][n]) > M; n--) ;
			if (n < minn) minn = n;
			for (M = 0.05 * fabs(WF[v][n = CP[v] + 1]); (n < NumWFPoints ? fabs(WF[v][n]) > M : false); n++) ;
			if (n > maxn) maxn = n;
			for (v=0; v < nv; v++) for (k=0; k<N; k++)
			{
				Sum = K[FC][I][J][v][k] = 0.0;
				for (n = ri[v]; n <= CP[v]; n++) 
				  K[FC][I][J][v][k] += WF[v][n] * WF[v][n] * S[SMap[k]][n];
				for ( ; n < ra[v]; n++) Sum += WF[v][n] * WF[v][n] * S[SMap[k]][n];
				K[FC][I][J][v][k] = (K[FC][I][J][v][k] + F[v] * F[v] * Sum) * SQ[v];
				if (SMap[k] >= NumFreePar - NumFreeAdCorr) K[FC][I][J][v][k] *= MFAdCorr;
				//printf("I=%d, v=%d, J=%d, k=%d, R=%e, WF=%g, S=%g\n", 
					//   I, v, J, k, rmin + h * CP[v], WF[v][CP[v]], S[k][CP[v]]);
			}
		}
	}
	//printf("Vor delete\n");
	delete[] F;
	delete[] U;
	delete[] Y;
	Destroy(WF, Mv + 1);
	delete[] SQ;
	delete[] CP;
	delete[] ri;
	delete[] ra;
	delete[] RPot;
	/*if (WFSource != 0)
	{
		delete[] E2;
	    if (NAdCorr > 0) delete[] A2;
		delete[] U2;
		delete[] R2;
	}*/
	if (NAdCorr > 0) delete[] AdCorrPot;
	//printf("Nach delete\n");
	minR = rmin + h * minn;
	maxR = rmin + h * maxn;
	for (n=0; n < NL; n++) 
	{
		//printf("NL=%d, Iso=%d, v0=%d, J0=%d, J1=%d, v1=%d\n", 
			//   NL, QN[n][0], QN[n][1], QN[n][2], QN[n][3], QN[n][4]);
		if (QN[n][3] >= 0)
		{
			if (n>0 ? QN[n][3] != QN[n-1][3] : true) 
			{
				for (k=0, uE = ErrS = 0.0; k<N; k++) uK[k] = 0.0;
				for (m=n; m < n + QN[n][4]; m++)
				{
					for (k=0, Sq = Sig[m] * Sig[m]; k<N; k++) 
						uK[k] += Sq * K[QN[m][5]][QN[m][0]][QN[m][1]][QN[m][2]][k];
					uE += Sq * (EO[m] + E[QN[m][5]][QN[m][0]][QN[m][1]][QN[m][2]]);
					ErrS += Sq;
				}
				for (k=0, Sq = 1.0 / ErrS; k<N; k++) uK[k] *= Sq;
				uE *= Sq;
			}
			DE[n] = uE - EO[n] - (CE1 = E[QN[n][5]][QN[n][0]][QN[n][1]][QN[n][2]]);
		}
		else DE[n] = EO[n] - (CE1 = E[QN[n][5]][QN[n][0]][QN[n][1]][QN[n][2]]);
		//if (QN[n][3] != -1) DE[n] += (CE2 = E[QN[n][5]][QN[n][0]][QN[n][3]][QN[n][4]]);
		if (CE1 == 0.0) // || (QN[n][3] != -1 && CE2 == 0.0))
		{
			DE[n] = 0.0;
			for (k=0; k<N; k++) EQS[n][k] = 0.0;
		}
		else for (k=0; k<N; k++) 
		{
			EQS[n][k] = K[QN[n][5]][QN[n][0]][QN[n][1]][QN[n][2]][k] 
						- (QN[n][3] != -1 ? uK[k] : 0.0);
			//if (QN[n][3] != -1) EQS[n][k] -= K[QN[n][5]][QN[n][0]][QN[n][3]][QN[n][4]][k];
			//if (QN[n][2] == 0) printf(", EQS[%d][%d]=%g", n, k, EQS[n][k]);
		}
		//printf("DE[%d]=%g, Eo=%g, Ec=%d\n", n, DE[n], EO[n], CE1 - CE2);
		//printf("N=%d\n", N);
	}
	//printf("Ende calcEK\n");
}

void calcAdCorrPot(double* InPot, double* OutPot, double Rmin, double Step, int NPoints, double MIso1, double MIso2,
	               int NAdCorr, int TAdCorr, int PAdCorr, double RIso1AdCorr, double RIso2AdCorr, double AdCorr_Rm, 
				   double AdCorr_b, double *adCorr)
{
	int n, m;
	double R, x, T;
	double MR = 1.0 - (TAdCorr == 1 ? RIso1AdCorr / MIso1 : (TAdCorr == 2 ? RIso2AdCorr / MIso2 
							: RIso1AdCorr * RIso2AdCorr * (MIso1 + MIso2)  / (MIso1 * MIso2 * (RIso2AdCorr + RIso1AdCorr))));
	if (NAdCorr == 0)
	{
		for (n=0; n < NPoints; n++) OutPot[n] = InPot[n];
		return;
	}
	for (n=0, R = Rmin; n < NPoints; n++, R += Step)
	{
		T = 0.0;
		if (NAdCorr > 1)
		{
			x = (R - AdCorr_Rm) / (R + AdCorr_b * AdCorr_Rm);
			for (m = NAdCorr -1; m > 0; m--)
			{
				T += adCorr[m];
				T *= x;
			}
		}
		T += adCorr[0];
		OutPot[n] = InPot[n] + T * MR * pow(2.0 * AdCorr_Rm / (R + AdCorr_Rm), PAdCorr);
	}
}

void cRotPot(int J, int Iso, double *U, double *R, double minR, double h, int n, int Omega, double q_e, double q_f, double *IsoF, 
			 double SpinRRm, double SpinRb, int NSpinRGamma, double *SpinRGamma, double SpinRR1, double SpinRR2, int c, int FC)
{
	//printf("Beginn CRotPot\n");
	int r;
	double p, O = Omega * Omega, q = (c==0 ? q_e : q_f);
	if (FC == 0) for (r=0, p=minR; r < n; r++, p += h) R[r] = U[r] + double(J * (J + 1) - O) * (1.0 + q) / (p * p * IsoF[Iso]);
	else if (FC == 1) for (r=0, p=minR; r < n; r++, p += h) 
		R[r] = U[r] + (double(J) * (double(J + 1) + SpinRotation(p, SpinRRm, SpinRb, NSpinRGamma, SpinRGamma, SpinRR1, SpinRR2)) 
		     - double(O)) / (p * p * IsoF[Iso]);
	else for (r=0, p=minR; r < n; r++, p += h) 
		R[r] = U[r] + (double(J+1) * (double(J) - SpinRotation(p, SpinRRm, SpinRb, NSpinRGamma, SpinRGamma, SpinRR1, SpinRR2)) 
		     - double(O)) / (p * p * IsoF[Iso]);
}

void NumerovCooley(double *&Y, double M, int &nv, double *E, double **WF, 
                              int *ri, int *ra, int *CP, double *F, double *SQ, const int n, double *U, double *MinErr)
{
	//printf("Beginn NumerovCooley\n");
	int i, hv, mv, v, fv, xm = n, lxm, lxM = 0, xM = n-2, Mv = nv - 1;//, raM;
	double Min = U[0], lMin = Min, SE, DE, wE, lwE, lDE, Max = U[n-1], lMax = Max, uDE, *tY/*, *x, *y, *sig*/;
	//bool isShapeRes;
	//Lorentzian ShapeRes;
	for (i=1; i<n-1; i++) 
	{
		if (U[i] < U[i-1] && U[i] <= U[i+1] && U[i] < Min) 
		{
			if (Max - Min > 1.0)
			{
				lMin = Min;
				lMax = Max; 
				lxm = xm;
				lxM = xM;
			}
			Max = Min = U[i];
			xm = i;
		}
		if (i > xm && U[i] > Max)
		{
			xM = i;
			Max = U[i];
		}
	}
	if (Max - Min <= 1.0 && lMax - lMin > 1.0)
	{
		Min = lMin;
		Max = lMax;
		xm = lxm;
		xM = lxM;
	}
	if (Min >= Max)
	{
		nv = 0;
		return;
	}
	SE = 0.5 * (Max + Min);
	//printf("Min=%f, Max=%f, SE=%f\n", Min, Max, SE);
	getMinMaxInt(U, SE, M, n, ri[0], ra[0]);
	if (SE > U[n-1]) ra[0] = xM;
	Numerov(M, SE, DE, hv, ra[0] - ri[0], U + ri[0], WF[0] + ri[0], CP[0], F[0], SQ[0], Y);
	//for (i=17000; i<n; i++) printf("E=%f, U=%f, WF=%f\n", SE, U[i], WF[0][i]);
	if (hv == 0) hv = 1;
	lwE = wE = (SE - Min) / hv;
	getMinMaxInt(U, Max, M, xM + 1, ri[0], ra[0]);
	Numerov(M, Max, DE, mv, ra[0] - ri[0], U + ri[0], WF[0], CP[0], F[0], SQ[0], Y);
	if (DE > 1e-7) mv--;
	//printf("hv=%d, mv=%d\n", hv, mv);
	//printf("ri=%d, ra=%d, n=%d\n", ri[0], ra[0], n);
	if (mv > Mv) 
	{
		//printf("mv=%d > Maxv=%d\n", mv, Mv);
		mv = Mv;
	}
	SE = Min - 0.5 * wE;
	//printf("Min=%f, wE=%f, M=%g\n", Min, wE, M);
	for (v=0; v<=mv; v++) 
	{
		//printf("v=%d, mv=%d\n", v, mv);
		SE += 2.0 * wE - lwE;
		if (v>0) if (SE < E[v-1]) SE = 0.5 * (Max + E[v-1]);
		lwE = wE;
		for (i=0; i<n; i++) WF[v][i] = 0.0;
		getMinMaxInt(U, SE, M, n, ri[v], ra[v]);
		if (ra[v] > xM) ra[v] = xM + 1;
		//isShapeRes = U[raM = ra[v]] <= SE;
		for (DE = 500.0, lDE = 1000.0, i=0; (fabs(DE) > (MinErr != 0 ? MinErr[v] : 1e-3) || fabs(lDE) > 0.1) && i < 30; i++)
		{
			if (SE > Max) SE = Max - 1e-5;
			lDE = DE;
			/*if (v==63)
			{
				delete[] WF[63];
				WF[63] = new double[NumWFPoints];
			}*/
			/*if (isShapeRes)
			{
				while (U[ra[v]] > SE) ra[v]++;
				while (U[ra[v]] <= SE) ra[v]--;
				if (ra[v] <= xM) ra[v] = xM + 1;
			}
			else if (U[ra[v]] <= SE) while (U[ra[v]] <= SE) ra[v]--;
			else if (ra[v] < raM)
			{
				while(U[ra[v]] > SE && ra[v] < raM) ra[v]++;
				if (U[ra[v]] <= SE) ra[v]--; 
			}*/
			//if (v==61) *DebugStream << "v''=61, SE=" << QString::number(SE, 'f', 10);
			Numerov(M, SE, DE, fv, ra[v] - ri[v], U + ri[v], WF[v] + ri[v], CP[v], F[v], SQ[v], Y);
			//if (v==61) *DebugStream << ", DE=" << QString::number(DE, 'f', 10) << '\n';
			CP[v] += ri[v];
			//printf("SE=%f, U[n-1]=%f, xM=%d, xo=%d, xi=%d\n", SE, U[n-1], xM, ra[v], ri[v]);
			/*if (v==63)
			{
				
				delete[] WF[63];
				WF[63] = new double[NumWFPoints];
						//break;
			}*/
			//printf("SE=%f, DE=%f, lDE=%f, wE=%f, lwE=%f, v=%d, fv=%d\n", SE, DE, lDE, wE, lwE, v, fv);
			if (fv < v)
			{
				lDE = 1000.0;
				DE = 500.0;
				SE += wE;
				while (SE >= Max) SE -= (wE *= 0.5);
				continue;
			}
			if (fv > v)
			{
				lDE = 1000.0;
				DE = 500.0;
				SE -= (wE *= 0.5);
				continue;
			}
			if (fabs(DE) > fabs(lDE))
			{
				if (DE > 0.0) SE += 0.5 * wE;
				else SE -= 0.5 * wE;
				DE = wE;
			}
			else if (SE + DE > Max)
			{
				printf("v=%d, (SE=%f) + (DE=%f) > Max=%f\n", v, SE, DE, Max); 
				for (uDE = 0.5 * DE; SE + uDE > Max; uDE *= 0.5) ;
				SE += uDE;
			} 
			else 
			{
				SE += DE;
				wE = (v > 0 ? SE - E[v-1] : 2.0 * (SE - Min));
			}
			//if (v==63) printf("SE=%f, DE=%f\n", SE, DE);
			if (isnan(SE) || isinf(SE) || isnan(DE) || isinf(DE))
			{
				SE = 0.0;
				printf("Invalid numbers, v=%d\n", v);
				break;
			}
		}
		if (i==30 || SE == 0.0) 
		{
			printf("Maximale Anzahl Iterationen erreicht, v=%d\n", v);
			mv = v - 1;
			break;
		}
		E[v] = SE;
		
		/*if (v >= 60)
		{
			double h = (rmax - rmin) / (NumWFPoints - 1);
			*DebugStream << QString("v=%1: R_min=%2, R_max=%3\n").arg(v).arg(rmin + h * double(ri[v]), 0, 'f', 5).arg(rmin + h * double(ra[v]), 0, 'f', 5);
		}*/
		
		if (SE > U[n-1]) 
		{
			if (Max - SE < 10.0)
			{
				double Err = (MinErr != 0 ? MinErr[v] : 1e-3), MR = 0.0, lMR = -1.0;
				for(DE = SE; MR > lMR; DE -= Err)
				{
					lMR = MR;
					scatWaveFunc(U, DE, M, ri[v], n, Y);
					MR = WFinoutMaxRatio(Y, ri[v], n);
					if (MR > lMR && lMR > 0.0)
					{
						tY = WF[v];
						WF[v] = Y;
						Y = tY;
					}
				}
				if (SE - DE < 2.1 * Err)
				{
					for (MR = lMR, DE = SE + Err; MR >= lMR; SE += Err)
					{
						lMR = MR;
						scatWaveFunc(U, SE, M, ri[v], n, WF[v]);
						MR = WFinoutMaxRatio(WF[v], ri[v], n);
						if (MR > lMR)
						{
							tY = WF[v];
							WF[v] = Y;
							Y = tY;
						}
					}
					E[v] = SE - 2.0 * Err;
				}
				else E[v] = DE + 2.0 * Err;
				F[v] = 1.0;
				while (ra[v] < n ? fabs(WF[v][ra[v]-1]) > fabs(WF[v][ra[v]]) : false) ra[v]++;
                for (SQ[v] = 0.0, i = ri[v]; i < ra[v]; i++) SQ[v] += WF[v][i] * WF[v][i];
				SQ[v] = 1.0 / SQ[v];
				//printf("ShapeRes: E[%d]=%f, SE=%f, DE=%f, Err=%f, DE/Err = %f\n", v, E[v], SE, SE - E[v], Err, (SE - E[v]) / Err);
				
			}
			/*double W, FQS, yf;
			if ((W = WF[v][CP[v] + 1] / WF[v][xM]) < 1e4)
			{
				for (i=0, x = new double[5], y = new double[5], sig = new double[5]; i<5; i++) sig[i] = 1e4;
				scatWaveFunc(U, x[2] = SE, M, ri[v], n, WF[v]);
				y[2] = WFinoutMaxRatio(WF[v], ri[v], n);
				scatWaveFunc(U, x[0] = SE - 0.01, M, ri[v], n, WF[v]);
				y[0] = WFinoutMaxRatio(WF[v], ri[v], n);
				scatWaveFunc(U, x[1] = SE - 0.005, M, ri[v], n, WF[v]);
				y[1] = WFinoutMaxRatio(WF[v], ri[v], n);
				while (y[1] < 0.25 * y[2])
				{
					y[0] = y[1];
					x[0] = x[1];
					scatWaveFunc(U, x[1] = 0.5 * (SE + x[0]), M, ri[v], n, WF[v]);
					y[1] = WFinoutMaxRatio(WF[v], ri[v], n);
				}
				scatWaveFunc(U, x[3] = 2 * SE - x[1], M, ri[v], n, WF[v]);
				y[3] = WFinoutMaxRatio(WF[v], ri[v], n);
				scatWaveFunc(U, x[4] = 2 * x[3] - SE, M, ri[v], n, WF[v]);
				y[4] = WFinoutMaxRatio(WF[v], ri[v], n);
				//printf("y[0]=%g, y[1]=%g, y[2]=%g, y[3]=%g, y[4]=%g\n", y[0], y[1], y[2], y[3], y[4]);
				ShapeRes.setData(x, y, sig, 5);
				FQS = ShapeRes.LevenbergMarquardt(100, 0.001);
				DE = ShapeRes.getE();
				while (fabs(DE - SE) > 1e-3)
				{
					scatWaveFunc(U, SE = DE, M, ri[v], n, WF[v]);
					ShapeRes.addDataPoint(DE, WFinoutMaxRatio(WF[v], ri[v], n), 1e4);
					FQS = ShapeRes.LevenbergMarquardt(100, 0.001);
					DE = ShapeRes.getE();
				}
				scatWaveFunc(U, DE, M, ri[v], n, WF[v]);
				yf = WFinoutMaxRatio(WF[v], ri[v], n);
				printf("ShapeRes: W=%g, DeltaE=%g, FQS=%f, y_fin=%g\n", W, DE - E[v], FQS, yf);
				E[v] = DE;
				F[v] = 1.0;
				ra[v] = n-1;
			}
			else ra[v] = xM + 1;*/
			
		}
		//printf("SQ[%d]=%g, ri=%d, ra=%d\n", v, 1.0/SQ[v], ri[v], ra[v]);
		//if (v==9) for (i=ri[v]; i<ri[v]+10; i++) printf("WF[9][%d]=%g\n", i, WF[9][i]);
	}
	//printf("WF[9][2208]=%g\n", WF[9][2208]);
	nv = mv + 1;
	//printf("Ende NumerovCooley, CP[0]=%d, nv=%d\n", CP[0], nv);
}

double WFinoutMaxRatio(double* WF, int bi, int bo)
{
	while (WF[bi + 1] > WF[bi]) bi++;
	if (WF[bo - 2] < WF[bo - 1])
	{
		for (bo--; WF[bo - 1] < WF[bo]; bo--) ;
		return sqr(WF[bi] / WF[bo]);
	}
	if (WF[bo - 2] > WF[bo -2])
	{
		for (bo--; WF[bo - 1] > WF[bo]; bo--) ;
		return sqr(WF[bi] / WF[bo]);
	}
	return sqr(WF[bi] / WF[bo - 1]);
}

void getMinMaxInt(double *U, double SE, double M, const int n, int &xi, int &xo)
{
	double LE, AE, NE;
	for (xi=1; (xi < n ? U[xi] > SE : false); xi++) ;
	for (LE = 0.0, AE = 1e-5, xi--; xi > 0 && AE < 1e7; xi--)
	{
		NE = AE * (2.0 + M * (U[xi] - SE)) - LE;
		LE = AE;
		AE = NE;
	}
	for (xo = n-2; (xo > xi ? U[xo] < SE : false); xo--) ;
	while (xo > xi ? U[xo] > SE : false) xo--;
	if (xo == xi)
	{
		xo = n;
		return;
	}
	//printf("xo=%d\n", xo);
	for (LE = 0.0, AE = 1e-5, xo++; xo < n && AE < 1e17; xo++)
	{
		NE = AE * (2.0 + M * (U[xo] - SE)) - LE;
		LE = AE;
		AE = NE;
	}
	//printf("getMinMaxInt: SE=%f, xi=%d, xo=%d\n", SE, xi, xo);
}

void Numerov(double M, double E, double &DE, int &v, const int n, double *U,
						double *WF, int &CP, double &F, double &SQ, double *Y)
{
	//printf("Beginn Numerov, n=%d, Iso=%d\n", n, Iso);
	int i, m;
	double N, eins = 1.0, PM = M * 0.0833333333333, ED;
	//printf("M=%f, h=%f, I=%f\n", M, h, IsoF[Iso]);
	bool gn = true;
	SQ = 0.0;
	Y[1] = Y[n-2] = 1e-5;
	Y[0] = Y[1] * exp(-sqrt(M*(U[0]-E)));
	Y[n-1] = Y[n-2] * exp(-sqrt(M*(U[n-1]-E)));
	WF[0] = Y[0] / (eins - PM * (U[0] - E));
	WF[n-1] = Y[n-1] / (eins - PM * (U[n-1] - E));
	//printf("WF[n-1]=%f, Y[n-1]=%f, M=%f, E=%f, U[n-1]=%f\n", WF[n-1], Y[n-1], M, E, U[n-1]);
	for (m=n-3; Y[m+1] > Y[m+2] && m>=0; m--)
	{
		WF[m+1] = Y[m+1] / (eins - PM * (ED = U[m+1] - E));
		Y[m] = M * ED * WF[m+1] + 2.0 * Y[m+1] - Y[m+2];
		//if (U[m+1] > E && U[m] <= E) printf("Ende Y[&d]=%g\n", m, Y[m]);
		//if ((isinf(Y[m]) || isnan(Y[m])) && !(isinf(Y[m+1]) || isnan(Y[m+1]))) 
			//printf("Y[%d]=%g, Y[%d]=%g\n", m, Y[m], m+1, Y[m+1]);
	}
	//printf("WF[n-2]=%f, n=%d, m=%d\n", WF[n-2], n, m);
	/*QFile DDatei("Debug.WF");
	DDatei.open(QIODevice::WriteOnly);
	QTextStream S(&DDatei);
	for (i=n-1; i>m; i--) S << "i=" << i << " Y=" << QString::number(Y[i], 'g', 10) << " WF=" 
				<< QString::number(WF[i], 'g', 10) 
				<< " ED= " << QString::number(U[i]-E, 'g', 10) << "\n";
	if (m >= n-2) printf("m=%d, n=%d\n", m, n);*/
	F = eins / Y[CP = ++m];
	for (i=2, v=0; i<=m; i++)
	{
		WF[i-1] = Y[i-1] / (eins - PM * (ED = U[i-1] - E));
		Y[i] = M * ED * WF[i-1] + 2.0 * Y[i-1] - Y[i-2];
		//if (U[i-1] > E && U[i] <= E) printf("Beginn Y[%d]=%g\n", i, Y[i]);
		//if ((isinf(Y[i]) || isnan(Y[i])) && !(isinf(Y[i-1]) || isnan(Y[i-1]))) 
			//printf("Y[%d]=%g, Y[%d]=%g\n", i, Y[i], i-1, Y[i-1]);
		if (gn && Y[i] < 0.0)
		{
			gn = false;
			v++;
		}
		else if (!gn && Y[i] > 0.0)
		{
			gn = true;
			v++;
		}
	}
	WF[m] = Y[m] / (eins - PM * (U[m] - E));
	for (i=0; i<=m; i++) SQ += WF[i] * WF[i];
	//if (v==9) printf("E=%g, WF[%d]=%g, SQ=%g\n", E, m, WF[m], SQ);
	F *= Y[m];
	for (i=m+1, N=0.0; i<n; i++) N += WF[i] * WF[i];
	SQ += N*F*F;
	SQ = eins / SQ;
	DE = ((2.0 * Y[m] - Y[m-1] - Y[m+1] * F) / M + (U[m] - E) * WF[m]) * SQ * Y[m];
	//if (isnan(DE)) 
		//printf("Y[%d]=%f, Y[%d]=%f, Y[%d]=%f, F=%f, M=%f, U[%d]=%f, E=%f, WF[%d]=%f, N=%f, n=%d\n",
			//   m-1, Y[m-1], m, Y[m], m+1, Y[m+1], F, M, m, U[m], E, m, WF[m], SQ, n);
	//printf("Ende Numerov, WF[n-1]=%f, CP=%d\n", WF[n-1], CP);
}


void scatWaveFunc(double* U, double E, double M, int pi, int NP, double* WF)
{
	double Y1 = 1e-5, PM = M * 0.0833333333333, Y0 = Y1 * exp(-sqrt(M*(U[pi]-E))), eins = 1.0, ED, Y2;
	int n;
	WF[n = pi] = Y0 / (eins - PM * (U[pi] - E));
	n++;
	WF[n] = Y1 / (eins - PM * (ED = U[n] - E));
	for (n++; n < NP; n++)
	{
		Y2 = M * ED * WF[n-1] + 2.0 * Y1 - Y0;
		WF[n] = Y2 / (eins - PM * (ED = U[n] - E));
		Y0 = Y1;
		Y1 = Y2;
	}
}

void CalcLMEQS(int NEQ, double** EQS, double* Sig, double* EDiff, double** LMEQS, double, int NumFreePar)
{
	int n, i, j;
	double wt;
	for (n=0; n < NumFreePar; n++) for (i=0; i <= NumFreePar; i++) LMEQS[n][i] = 0.0; 
	for (n=0; n < NEQ; n++) 
	{
		for (i=0; i < NumFreePar; i++)
		{
			wt = EQS[n][i] * Sig[n] * Sig[n];
			for (j=0; j < NumFreePar; j++) LMEQS[i][j] += wt * EQS[n][j];
			LMEQS[i][NumFreePar] += wt * EDiff[n];
		}
	}
}


DAddMaximizeFunctor::DAddMaximizeFunctor(Spline* ispline, double** iData, int iNJ, 
										 int iJStart, int iJStep, int iNE, 
										 double iEStart, double iERes, bool absorption)
{
	spline = ispline;
	Data = iData;
	NJ = iNJ;
	JStart = iJStart;
	JStep = iJStep;
	NE = iNE;
	EStart = iEStart;
	ERes = iERes;
	isAbsorption = absorption;
}

Doub DAddMaximizeFunctor::operator()(VecDoub x)
{
	int n, N = x.size(), J, JM = spline->getxN(), J0 = spline->getx0();
	int m = (J0 - JStart) / JStep, d;
	double R = 0.0, E, F1, F2;
	for (n=0; n<N; n++) spline->setyn(n, x[n]);
	spline->calcYss();
	for (J = J0; J <= JM; J += JStep, m++) 
	{
		E = 0.5 * NE + spline->gety(J) * ERes;
		d = int(E);
		F2 = E - double(d);
		F1 = 1.0 - F2;
		R += F1 * Data[m][d] + F2 * Data[m][d+1];
	}
	printf("R=%g\n", R);
	return (isAbsorption ? R : -R);
}

void DAddMaximize(ResidualFit* resFit, double** Data, int NJ, int JStart, int JStep, 
				  int NE, double EStart, double ERes, double Diff, bool absorption)
{
	int n, N = resFit->getNumSplines();
	for (n=0; n<N; n++)
	{
		Spline *spline = resFit->getSpline(n);
		spline->shifty(Diff);
		printf("Diff=%g\n", Diff);
		DAddMaximizeFunctor func(spline, Data, NJ, JStart, JStep, NE, EStart, ERes, absorption);
		Amoeba amoeba(1e-12);
		int p, NP = spline->getNPar();
		VecDoub point(NP);
		for (p=0; p < NP; p++) point[p] = spline->getyn(p);
		VecDoub R = amoeba.minimize<DAddMaximizeFunctor>(point, Diff, func);
		for (p=0; p < NP; p++) spline->setyn(p, R[p]);
		spline->calcYss();
	}
}

void LocalPerturbationFitFuncs(const Doub i_J, VecDoub_I &i_par, Doub &o_E_calc, VecDoub_O &o_dE_dParam)
{
    double J = abs(i_J), JStep = i_par[3];
    int NBandC = static_cast<int>(i_par[0]), Jstart = static_cast<int>(i_par[NBandC + 5]), i, Ji = static_cast<int>(J);
    double Hsq = i_par[4] * i_par[4], E1_calc = i_par[NBandC + 6 + (Ji - Jstart) / JStep], E2_calc = i_par[NBandC + 4];
    double JtJp1 = i_par[2] * (J * (J + 1.0) - i_par[1] * i_par[1]);
    for (i = NBandC + 3; i>4; --i)
    {
        E2_calc *= JtJp1;
        E2_calc += i_par[i];
    }
    double DeltaE_calc = E1_calc - E2_calc, DeltaE_calc_haSq = 0.25 * DeltaE_calc * DeltaE_calc, lRoot = sqrt(DeltaE_calc_haSq + Hsq), throughlRoot = 1.0 / lRoot;
    double signProd = i_J * DeltaE_calc;
    o_E_calc = 0.5 * (E1_calc + E2_calc) + (signProd > 0.0 ? lRoot : -lRoot);
#ifdef UseDebugLogStream
	*DebugStream << QString("H=%1\tT=%2\tB=%3\tJ=%4\tE=%5\n").arg(i_par[4], 0, 'f', 4).arg(i_par[5], 0, 'f', 4).arg(i_par[6], 0, 'f', 6)
						.arg(i_J, 0, 'f', 0).arg(o_E_calc, 0, 'f', 4);
#endif
    o_dE_dParam[4] = i_par[4] * (signProd > 0.0 ? throughlRoot : -throughlRoot);
    o_dE_dParam[5] = 0.5 + (signProd > 0.0 ? -0.25 : 0.25) * DeltaE_calc * throughlRoot;
    for (i=1; i < NBandC; ++i) o_dE_dParam[5+i] = o_dE_dParam[4+i] * JtJp1;
}

void FitLocalPerturbation(LocalPerturbation &io_perturbation, double i_pertBe)
{
    int NData = io_perturbation.GetNData(), minJ = -1, maxJ = -1, JStep = io_perturbation.GetJStep(), NBandC = io_perturbation.GetNBandC();
    VecDoub xx(NData), yy(NData), ssig(NData);
    for (int n=0; n < NData; ++n)
    {
        xx[n] = io_perturbation.GetDataJValue(n);
        if (xx[n] < minJ || minJ < 0) minJ = static_cast<int>(xx[n]);
        if (xx[n] > maxJ || maxJ < 0) maxJ = static_cast<int>(xx[n]);
        yy[n] = io_perturbation.GetDataEobs(n);
        ssig[n] = io_perturbation.GetDataUnc(n);
    }
    int NJ = (maxJ - minJ) / JStep + 1;
    double data[NJ], unc;
    for (int n=0; n < NJ; ++n) data[n] = 0.0;
    for (int n=0; n < NData; ++n) data[(static_cast<int>(xx[n]) - minJ) / JStep] = io_perturbation.GetDataECalc(n);
    VecDoub aa(6 + NBandC + NJ);
    io_perturbation.GetH12(aa[4], unc);
    for (int n=0; n < NBandC; ++n) io_perturbation.GetBandC(n, aa[5+n], unc);
    if (i_pertBe > 0.0)
    {
        double JPert = io_perturbation.GetCenter(), EPert = io_perturbation.GetPointCalcBand(JPert);
        io_perturbation.SetBandC(1, i_pertBe, 0.0);
        aa[5] += EPert - io_perturbation.GetPointCalcBand(JPert);
    }
    Fitmrq Lmq(xx, yy, ssig, aa, LocalPerturbationFitFuncs);
    Lmq.hold(0, NBandC);
    Lmq.hold(1, io_perturbation.GetOmega());
    Lmq.hold(2, io_perturbation.GetIsoF());
    Lmq.hold(3, JStep);
    if (i_pertBe > 0.0) Lmq.hold(6, i_pertBe);
    Lmq.hold(5 + NBandC, minJ);
    for (int n=0; n < NJ; ++n) Lmq.hold(6 + NBandC + n, data[n]);
    Lmq.fit();
#ifdef UseDebugLogStream
	*DebugStream << QString("FQS=%1\n").arg(Lmq.chisq, 0, 'f', 3);
#endif
    io_perturbation.SetChisq(Lmq.chisq);
    io_perturbation.SetH12(Lmq.a[4], Lmq.covar[4][4] * Lmq.covar[4][4]);
    for (int n=0; n < NBandC; ++n) io_perturbation.SetBandC(n, Lmq.a[5+n], Lmq.covar[5+n][5+n] * Lmq.covar[5+n][5+n]);
}

void AutoFitLocalPerturbation(LocalPerturbation &io_perturbation, const int i_maxBandC)
{
    bool NotTooMuch;
    LocalPerturbation fitPerturbation(io_perturbation);
    int m, n, J, Jl = fitPerturbation.GetDataJValue(0), nm = 1, nM = 1, minJ = Jl, maxJ = fitPerturbation.GetDataJValue(fitPerturbation.GetNData() - 1);
    double value, stdDev, maxDevNeighJ = 0.0, minDevSameJ = 9e99, unc, uncl = fitPerturbation.GetDataUnc(0), cdev, dev, devl = fitPerturbation.GetDataEobs(0) - fitPerturbation.GetDataECalc(0);
    double JPert, EPert, H, B, T, O = fitPerturbation.GetOmega(), IsoF = fitPerturbation.GetIsoF();
    io_perturbation.SetChisq(9e99);
    while (fitPerturbation.GetNBandC() > 2) fitPerturbation.RemoveLastBandC();
    for (n=1; n < fitPerturbation.GetNData(); ++n)
    {
        J = fitPerturbation.GetDataJValue(n);
        unc = fitPerturbation.GetDataUnc(n);
        dev = fitPerturbation.GetDataEobs(n) - fitPerturbation.GetDataECalc(n);
        cdev = abs(dev - devl);
        if (J == Jl)
        {
            if (cdev > 4.0 * unc && cdev > 4.0 * uncl && cdev < minDevSameJ)
            {
                minDevSameJ = cdev;
                nm = n;
            }
        }
        else
        {
            if (cdev > maxDevNeighJ)
            {
                maxDevNeighJ = cdev;
                nM = n;
            }
        }
        Jl = J;
        uncl = unc;
        devl = dev;
    }
    if (minDevSameJ < 9e99) nM = nm;
    JPert = 0.5 * (fitPerturbation.GetDataJValue(nM) + fitPerturbation.GetDataJValue(nM - 1));
    EPert = 0.5 * (fitPerturbation.GetDataEobs(nM) + fitPerturbation.GetDataEobs(nM - 1));
    double EDiffLeft = fitPerturbation.GetDataEobs(nM - 1) - fitPerturbation.GetDataEobs(0) - fitPerturbation.GetDataECalc(nM - 1) + fitPerturbation.GetDataECalc(0);
    if (minDevSameJ < 9e99)
    {
        H = 0.5 * minDevSameJ;
        double E2 = 0.0, J2 = 0.0;
        if (nm >= 3 && fitPerturbation.GetDataJValue(nm - 3) == fitPerturbation.GetDataJValue(nm - 2))
        {
            J2 = fitPerturbation.GetDataJValue(nm - 3);
            E2 = fitPerturbation.GetDataEobs(nm - 3) + fitPerturbation.GetDataEobs(nm - 2) - fitPerturbation.GetDataECalc(nm - 3);
        }
        else if (nm <= fitPerturbation.GetNData() - 3 && fitPerturbation.GetDataJValue(nm + 1) == fitPerturbation.GetDataJValue(nm + 2))
        {
            E2 = fitPerturbation.GetDataEobs(nm + 1) + fitPerturbation.GetDataEobs(nm + 2) - fitPerturbation.GetDataECalc(nm + 1);
            J2 = fitPerturbation.GetDataJValue(nm + 1);
        }
        if (E2 != 0.0 && J2 != 0.0)
        {
            B = (2.0 * EPert - fitPerturbation.GetDataECalc(nm) - E2) / (IsoF * (JPert * (JPert + 1.0) - J2 * (J2 + 1.0)));
            T = E2 - IsoF * (J2 * (J2 + 1.0) - O * O) * B;
        }
        else
        {
            double E1 = 2.0 * EPert - fitPerturbation.GetDataECalc(nm), E2a = fitPerturbation.GetDataEobs(nm - 2) + (EDiffLeft > 0.0 ? -H : H);
            double E2b = fitPerturbation.GetDataEobs(nm + 1) + (EDiffLeft > 0 ? H : -H);
            double J2a = fitPerturbation.GetDataJValue(nm - 2), J2b = fitPerturbation.GetDataJValue(nm + 1);
            double Ba = (E1 - E2a) / (IsoF * (JPert * (JPert + 1.0) - J2a * (J2a + 1.0)));
            double Bb = (E1 - E2b) / (IsoF * (JPert * (JPert + 1.0) - J2b * (J2b + 1.0)));
            B = 0.5 * (Ba + Bb);
            T = E1 - IsoF * (JPert * (JPert + 1.0) - O * O) * B;
        }
        int lastUp = -1, lastDown = -1;
        if (EDiffLeft > 0.0)
        {
            for (n = nm; n >= 0; --n)
            {
                double Jval = fitPerturbation.GetDataJValue(n);
                bool isDown = fitPerturbation.GetDataEobs(n) < T + IsoF * (Jval * (Jval + 1.0) - O * O) * B;
                if (isDown) lastDown = n;
                fitPerturbation.SetDataIsUp(n, !isDown);
            }
            for (n = nm - 1; n < fitPerturbation.GetNData(); ++n)
            {
                double Jval = fitPerturbation.GetDataJValue(n);
                bool isUp = fitPerturbation.GetDataEobs(n) > T + IsoF * (Jval * (Jval + 1.0) - O * O) * B;
                if (isUp) lastUp = n;
                fitPerturbation.SetDataIsUp(n, isUp);
            }
        }
        else
        {
            for (n = nm; n >= 0; --n)
            {
                double Jval = fitPerturbation.GetDataJValue(n);
                bool isUp = fitPerturbation.GetDataEobs(n) > T + IsoF * (Jval * (Jval + 1.0) - O * O) * B;
                if (isUp) lastUp = n;
                fitPerturbation.SetDataIsUp(n, isUp);
            }
            for (n = nm - 1; n < fitPerturbation.GetNData(); ++n)
            {
                double Jval = fitPerturbation.GetDataJValue(n);
                bool isDown = fitPerturbation.GetDataEobs(n) < T + IsoF * (Jval * (Jval + 1.0) - O * O) * B;
                if (isDown) lastDown = n;
                fitPerturbation.SetDataIsUp(n, !isDown);
            }
        }
        fitPerturbation.SetLastLevels(lastUp, lastDown);
    }
    else
    {
        double E1 = fitPerturbation.GetDataEobs(nM - 1), E2 = fitPerturbation.GetDataEobs(nM);
        double J1 = fitPerturbation.GetDataJValue(nM - 1), J2 = fitPerturbation.GetDataJValue(nM);
        B = (EDiffLeft > 0.0 ? 1.2 * fitPerturbation.GetBeCurState() : 0.85 * fitPerturbation.GetBeCurState());
        T = 0.5 * (E1 + E2) - (0.5 * (J1 * (J1 + 1) + J2 * (J2 + 1)) - O*O) * IsoF * B;
        H = abs(T - E1 + (J1 * (J1 + 1) - O*O) * IsoF * B);
        if (EDiffLeft > 0.0) fitPerturbation.SetLastLevels(nm - 1, nm);
        else fitPerturbation.SetLastLevels(nm, nm - 1);
        for (n=0; n < fitPerturbation.GetNData(); ++n)
        {
            double Jval = fitPerturbation.GetDataJValue(n);
            fitPerturbation.SetDataIsUp(n, fitPerturbation.GetDataEobs(n) > T + IsoF * (Jval * (Jval + 1.0) - O * O) * B);
        }
    }
    fitPerturbation.SetH12(H, 0.0);
    fitPerturbation.SetBandC(0, T, 0.0);
    fitPerturbation.SetBandC(1, B, 0.0);
    fitPerturbation.SetFirstUp(EDiffLeft > 0.0);
    for (NotTooMuch = true, n=1; NotTooMuch && n <= i_maxBandC; ++n)
    {
        if (n>1) fitPerturbation.SetBandC(n, 0.0, 0.0);
        //FitLocalPerturbation(fitPerturbation, -1.0);
        fitPerturbation.LevenbergMarquardt(100, 1e-6);
        if (io_perturbation.GetFitSigma() <= fitPerturbation.GetFitSigma()) NotTooMuch = false;
        fitPerturbation.GetH12(value, stdDev);
        if (abs(stdDev) > 0.1 * abs(value)) NotTooMuch = false;
        for (m=0; m <= n; ++m)
        {
            fitPerturbation.GetBandC(n, value, stdDev);
            if (abs(stdDev) > 0.1 * abs(value)) NotTooMuch = false;
        }
        if (NotTooMuch || io_perturbation.GetNBandC() < 2) io_perturbation = fitPerturbation;
    }
    while (JPert <= maxJ && (EDiffLeft > 0 ? io_perturbation.GetPointCalcBand(JPert) < io_perturbation.GetPointUpOBCa(JPert) : io_perturbation.GetPointCalcBand(JPert) > io_perturbation.GetPointUpOBCa(JPert)))
        JPert += 1.0;
    while (JPert >= minJ && (EDiffLeft > 0 ? io_perturbation.GetPointCalcBand(JPert) > io_perturbation.GetPointUpOBCa(JPert) : io_perturbation.GetPointCalcBand(JPert) < io_perturbation.GetPointUpOBCa(JPert)))
        JPert -= 1.0;
    JPert += 0.5;
    if (JPert < io_perturbation.GetDataJValue(0)) JPert = io_perturbation.GetDataJValue(0) + 0.5;
    else if (JPert > io_perturbation.GetDataJValue(io_perturbation.GetNData() - 1)) JPert = io_perturbation.GetDataJValue(io_perturbation.GetNData() - 1) - 0.5;
    io_perturbation.SetCenter(JPert);
}

void FitBandConstants(double *o_BandConstants, double *o_err, int i_NBandC, double i_IsoF, double i_Omega, double *i_J, double *i_E, double *i_unc, int i_NData)
{
    MatDoub a(i_NBandC, i_NBandC), C(i_NBandC, i_NBandC);
    VecDoub b(i_NBandC), x(i_NBandC);
    for (int n=0; n < i_NBandC; ++n)
    {
        for (int m=0; m < i_NBandC; ++m) a[n][m] = 0.0;
        b[n] = 0.0;
    }
    double OmegaSq = i_Omega * i_Omega;
    for (int i=0; i < i_NData; ++i)
    {
        double sig = 1.0 / (i_unc[i] * i_unc[i]), JF = i_IsoF * (i_J[i] * (i_J[i] + 1.0) - OmegaSq);
        VecDoub X(i_NBandC);
        X[0] = 1.0;
        for (int n=1; n < i_NBandC; ++n) X[n] = JF * X[n-1];
        for (int k=0; k < i_NBandC; ++k)
        {
            b[k] += sig * i_E[i] * X[k];
            for (int j=0; j < i_NBandC; ++j) a[k][j] += sig * X[j] * X[k];
        }
    }
    LUdcmp LU(a);
    LU.solve(b, x);
    double lFQS = calcEQSFQS(a, b, x);
    LU.mprove(b, x);
    double FQS = calcEQSFQS(a, b, x);
    while ((lFQS - FQS) / FQS > 1e-7)
    {
        lFQS = FQS;
        LU.mprove(b, x);
        FQS = calcEQSFQS(a, b, x);
    }
    LU.inverse(C);
    for (int n=0; n < i_NBandC; ++n)
    {
        o_BandConstants[n] = x[n];
        o_err[n] = sqrt(C[n][n]);
    }
}

void calcFFT(const double *const realInputData, const int N, const double delta, double **const realOutput, double **const imaginaryOutput)
{
	int n, m, DL = 2*N;
	Doub* data = new Doub[DL];
	for (n=0, m=-1; n<N; ++n)
	{
		data[++m] = realInputData[n];
		data[++m] = 0.0;
	}
	four1(data, N, 1);
	double f, fStep = 1.0 / (N * delta);
	realOutput[0][0] = imaginaryOutput[0][0] = 0.0;
	realOutput[0][1] = data[0];
	imaginaryOutput[0][1] = data[1];
	for (n=2, m=1, f = fStep; n < DL; n+=2, ++m, f += fStep)
	{
		realOutput[m][0] = imaginaryOutput[m][0] = f;
		realOutput[m][1] = data[n] + data[DL-n];
		imaginaryOutput[m][1] = data[n+1] - data[DL-n+1];
	}
	realOutput[N][0] = imaginaryOutput[N][0] = f;
	realOutput[N][1] = data[N];
	imaginaryOutput[N][1] = data[N+1];
}
