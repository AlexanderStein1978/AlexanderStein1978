//
// C++ Interface: PotWorker
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2007 - 2019
//
// Copyright: See README file that comes with this source code
//
//


#ifndef POTWORKER_H
#define POTWORKER_H


#include <QObject>

#include "fitresult.h"
#include "wantedcoeffvalue.h"
#include "constants.h"
#include "termenergy.h"
#include "fit.h"


class QMutex;

class PotFit;
class AnaPot;
class IsoTab;

struct SplinePoint;
struct TermEnergy;
struct TableLine;
struct WantedCoeffValue;
struct PotentialData;


class PotWorker : public QObject
{
	Q_OBJECT
	
public:
	PotWorker(PotFit *Fit, PotentialType PotType = NoPotential);
	PotWorker(const PotWorker &C);
	~PotWorker();
	void MonteCarloSimIteration(double *MCSTab, double UncFact = 1.0);
    double FitAnaPot(int NumWFPoints, bool robustWeighting = false, bool UseSVD = true, double threshhold = 1e-10, bool useLevenbergMarquardt = true,
					 int MaxIt = -1, bool adjustTEWeightFact = false, double hCi_cP = 0.0, bool FitRestarted = false);
    void FitSplinePot(int NumWFPoints, int MaxIt = 1, double threshFakt = 0.9999);
	void FitSplinePotToFunction(double *E, int NPoints, double MinR, double MaxR);
	double FitTangToenniesPot();
    double FitIsoMass(int Iso, double InitStepSize, double &RIsoMass, int NumWFPoints);
	void SetSplinePar(double *Par, double ****E, double *****K, int *mJ, int **mv, int Mv, int N, 
				      int NL, int **QN, double *Sig, double **EQS, double *EO, double *DE, 
                      double &minR, double &maxR, int NumWFPoints);
	AnaPot *getAnaPotForReading();
	AnaPot *getAnaPotForWriting();
	bool setAdCorr(int NAdCorr, double *adCorr, int TAdCorr, int PAdCorr, double RIso1, double RIso2,
				   double Rm = 0.0, double b = 1e-30, bool *AdCorrFree = 0);
	void getAdCorrForReading(int &NAdCorr, double *&adCorr, int &TAdCorr, int &PAdCorr, double &RIso1, double &RIso2, double &Rm, 
							 double &b, bool **AdCorrFree = 0, bool alreadyLocked = false);
	void getAdCorrForWriting(int &NAdCorr, double *&adCorr, int &TAdCorr, int &PAdCorr, double &RIso1, double &RIso2, double &Rm, 
							 double &b);
	bool setSplinePot(int NSplinePoints, SplinePoint *points, int NLRC, int *pLRC, double *LRC, double iA, double iO, double iExp,
					  double LongRangeF = 0.0);
	void getSplinePotForReading(int &NSplinePoints, SplinePoint *&points, int &NLRC, int *&pLRC, double *&LRC, double &iA, double &iO, 
								double &Exp);
	void getSplinePotForWriting(int &NSplinePoints, SplinePoint *&points, int &NLRC, int *&pLRC, double *&LRC, double &iA, double &iO, 
								double &Exp);
	void getLRCoeffForReading(int &N, int *&pLRC, double *&LRC, double *Ra = 0, bool alreadyLocked = false);
	void getLRCoeffForWriting(int &N, int *&pLRC, double *&LRC, bool **LRCFree = 0, double *Ra = 0);
	SplinePoint *getSplinePoints(int &numPoints);
	void setFitData(int *mJ, int **mv, int Mv, int nEs, int *NEL, TermEnergy **EL, int nGS, int *nGP, int *nGL, int **PS, 
					TableLine **GL, bool*ES, int MaxJ, int Maxv);
	int getRefIso();
	void getBadList(QString **&Data, int &NR, int &NC);
	void getDeviations(int &N, int *&Row, double *&dev, double *&RWErr);
	void getFitResult(double &FQS, double &WeightSum, int &NWeightF, int *&WeightFRow, double *&WeightF, int &NumFreePar);
	void getSplineMinMaxR(double &Min, double &Max);
    void showBadList(QString **&Data, int &NR, int &NC, int NumWFPoints);
	void badListChanged(QString *Items, int N, TableLine &L);
	void getSpinRGamma(int &rNSpinRGamma, double *&rSpinRGamma, double &rSpinRR1, double &rSpinRR2, double &rSpinRb, double &rSpinRRm);
	void setSpinRGamma(int NSpinRGamma, double *SpinRGamma, double SpinRR1, double SpinRR2, double SpinRb, double SpinRRm);
	void estMinErr(double ****&MinErr, int NFC, int *mJ, int **mv, int nEs, int *NEL, TermEnergy **EL, 
				   int nGS, int *nGL, TableLine **GL);
	void CreateEK(double ****&E, double *****&K, int *mJ, int **mv, int N);
	void CreateFitEQS(int &N, double **&EQS, double *&EO, double *&DE, double *&Sig, double *&Err, int **&QN, int NP, int nEs, 
					  int *NEL, TermEnergy **EL, int nGS, int *nGP, int *nGL, int **PS, TableLine **GL, int AdRows = 0);
    void calcE(double ****E, int *mJ, int **mv, int Mv, int NumWFPoints, double ****MinErr = 0, bool ****SFQSU = 0, double SFQSRad = 0.0);
	double *getPoints(double Rmin, double Rmax, int numPoints, int FC = 0);
    double *get_dVdR(const double Rmin, const double Rmax, const int numPoints) const;
	double getPoint(double R, int FC = 0);
    void getDiagFuncs(double *&WFS, double *&WWFS, int NumWFPoints);
	bool setMolData(IsoTab *IsoT, double *MU, double *IsoF);
	void takeMolData(IsoTab *&IsoT, double *&MU, double *&IsoF);
	void fastWaveFunc(double *U, double E, double M, double h, int NP, double *WF, double &Ri, 
					  double &Ra, double &ISq, double &F, int &i, int &c, int &a);
	double calcFQS(double ****E, int nES, int *NEL, TermEnergy **EL, int nGS, int *nGP, int *nGL, 
                   int **PS, TableLine **GL, int NumWFPoints, bool calcDiag = false, int MJ = 0, int Mv = 0, double *****K = 0,
	   			   int *rv = 0, double *Sig = 0, double TermEnergieF = 1.0, double hCi_cP = 0.0, bool ****SFQSU = 0);
    double getFQS(int NumWFPoints, bool calcDiagFuncs = false, double *StdDev = 0, bool Thread = false, double SFQSRad = 0.0,
				  double *sigma = 0, double *FQS_PAL = 0);
	double *getMCSE(int &N, int *&Rows);
	void unlock();
	void unlockBadList();
	virtual void cdConnectLR(int firstC);
	virtual void cdConnectLR1C();
	virtual void cdConnectSR();
	void getS(double &Rmin, double &Rmax, int &NPoints, int &NC, double **&S);
	double getMinR(double E, double Prec);
	void getMaximum(double &R, double &U);
	virtual void getMinimum(double &R, double &U);
	bool addPoint(double x, double y);
	bool movePoint(int &n, double newx, double newy);
	bool removePoint(int n);
	bool isFitStopped();
	PotWorker* getCopy();
	bool setOmega(int nOmega, double S);
	bool setUinf(double nUinf);
	virtual PotWorker* scalePotential(double newRe, double newDe);
	virtual void shift(double E);
	void getSFQSU(bool ****&SFQSU, int *&mJ, int **&mv);
	void CalcAdCorrPot(double* InPot, double* OutPot, double Rmin, double Step, int NPoints, double MIso1, double MIso2);
	virtual PotentialData *getPotentialData();
	void calcyss();
	virtual void setLRCoeff(double R, int numCoefficients, int *Exponents, double *Coefficients, bool *LRCFree = 0);
    virtual double getInnerConnectionRadius() const;
    virtual double getOuterConnectionRadius() const;
	
	inline void setSFQSU(bool ****nSFQSU)
	{
		SFQSU = nSFQSU;
	}
	
	inline FitResult getFitResult()
	{
		return Result;
	}
	
	inline bool isAdCorrA()
	{
		if (NAdCorr > 0) return true;
		return false;
	}
	
	inline bool isLambdaDoublingA()
	{
		return q_e - q_f != 0.0;
	}
	
	inline void getLambdaDoubling(double &qe, double &qf)
	{
		qe = q_e;
		qf = q_f;
	}
	
	inline void setLambdaDoubling(double qe, double qf)
	{
		q_e = qe;
		q_f = qf;
	}
	
	inline int getNSpinRGamma()
	{
		return NSpinRGamma;
	}
	
	inline int getNumFreePar()
	{
		return NumFreePar;
	}
	
	inline int getNIso()
	{
		return NIso;
	}
	
	inline double* getIsoF_ForReading()
	{
		return IsoF;
	}
	
	virtual inline double getUinf()
	{
		return Uinf;
	}
	
	inline IsoTab *getIsoT()
	{
		return IsoT;
	}
	
	inline bool isFitResultAv()
	{
		return dev != 0;
	}
	
	inline bool isFitDataAv()
	{
		return NEL != 0 || nGL != 0;
	}
	
	inline bool isAnaPotAv()
	{
		return Type == analyticalPotential;
	}
	
	inline bool isSplinePot()
	{
		return Type == SplinePotential;
	}
	
	inline bool isLPot()
	{
		return Type == MorseLongRangePotential;
	}
	
	inline bool isMTTPot()
	{
		return Type == ModifiedTangToenniesPotential;
	}
	
	inline void setCalcDiagFuncs(bool calc)
	{
		calcDiagFuncs = calc;
	}
	
	inline bool getCalcDiagFuncs()
	{
		return calcDiagFuncs;
	}
	
	inline void getmvmJ(int *&rmJ, int **&rmv)
	{
		rmJ = mJ;
		rmv = mv;
	}
	
	inline void getSplinePotPointandCoeffNumbers(int &NPoints, int &nLRC)
	{
		NPoints = numSplinePoints;
		nLRC = NLRC;
	}
	
	virtual inline double getIExp()
	{
		return iExp;
	}
	
	inline int *getSplinePLRC()
	{
		return PLRC;
	}
	
	inline bool isTTPot()
	{
		return Type == TangToenniesPotential;
	}
	
	virtual inline bool isExcIntAvailable()
	{
		return false;
	}
	
	inline bool getShortRVariable()
	{
		return shortRVariable;
	}
	
	inline void setLRCFree(bool *nLRCFree)
	{
		if (LRCfree != 0) delete[] LRCfree;
		LRCfree = nLRCFree;
	}
	
	inline void DestroySMap()
	{
		if (SMap != 0) delete[] SMap;
		SMap = 0;
	}
	
	inline void recalcSCLC()
	{
		reCalcLC = true;
		reCalcSC = true;
	}
	
	inline bool isLevenbergMarquardtRunning()
	{
		return LevenbergMarquardtRunning;
	}
	
	inline void setFit(PotFit *nFit)
	{
		Fit = nFit;
	}
	
	virtual inline bool *getLRCFree()
	{
		return LRCfree;
	}
	
	virtual inline bool isFSA()
	{
		return NSpinRGamma > 0;
	}
	
	inline PotentialType getType()
	{
		return Type;
	}
	
	inline double getSFQS()
	{
		return SFQS;
	}
	
	inline void setWantedValues(int numWantedValues, WantedCoeffValue *Values)
	{
		if (NumWantedValues > 0) delete[] WantedValues;
		WantedValues = Values;
		NumWantedValues = numWantedValues;
	}
	
	virtual inline void setWantedLRValues(int numValues, WantedCoeffValue *Values)
	{
		int n, p = numSplinePoints + 1;
		for (n=0; n < numValues; n++) Values[n].CoeffNum += p;
		setWantedValues(numValues, Values);
	}
	
	inline WantedCoeffValue *getWantedValues()
	{
		return WantedValues;
	}
	
	inline void CRotPot(int J, int Iso, double *U, double *R, double minR, double h, int n, int c=0, int FC = 0)
	{
		cRotPot(J, Iso, U, R, minR, h, n, Omega, q_e, q_f, IsoF, SpinRRm, SpinRb, NSpinRGamma, SpinRGamma, 
				SpinRR1, SpinRR2, c, FC);
	}
	
	inline void setR_o(double R_o)
	{
		Ro = R_o;
	}
	
public slots:
	
	void stopFit(bool Stop);
	
signals:
	void potentialImproved(PotentialData *Data);
	void potentialChanged();
	void badListChanged(double aFQS, double aLambda, double mLSFQS);
	void fitStopped();
	void fitRunning();
	void firstFCFCalculated(FitResult fitResult, double FQS);

protected:
	virtual void createSMap();
	virtual void updatePotential(double *C);
	virtual void restorePotential(double *bC);
	virtual void saveCoefficients(double *&bC);
	virtual int cAdRows();
	virtual void getLRCoeff(double &R, int &numCoefficients, int *&Exponents, double *&Coefficients);
	virtual void calc_hCi_cP_EQS(int NEQ, double **EQS, double *Sig, double *EDiff, double hCi_cP);
	virtual void calcS(int numPoints, double RMin, double RMax);
	virtual double Point(double R, int FC = 0);
	virtual void setCurValuesForWantedValues();
	virtual double *Points(double Rmin, double Rmax, int numPoints, int FC = 0);
    virtual double *dVdR(const double Rmin, const double Rmax, const int numPoints) const;
	
	PotentialType Type;
	SplinePoint *points;
	int numSplinePoints; 
	QMutex *Lock, *BLLock;
	int *SMap;
	int TAdCorr, NAdCorr, PAdCorr, NumFreePar, NumFreeLRC, NumFreeAdCorr, NSpinRGamma, NFitCoeff;
	double *adCorr, RIso1AdCorr, RIso2AdCorr, AdCorr_b, AdCorr_Rm, *AdCorrWF, *SpinRGamma, SpinRR1, SpinRR2, SpinRb, SpinRRm;
	PotFit *Fit;
	bool *adCorrFree, ShowBadList, shortRVariable, *LRCfree, calcDiagFuncs;
	int Omega, *mJ, **mv, Mv, nEs, *NEL, nGS, *nGP, *nGL, **PS, Maxv, MaxJ, NumS, NIso, BadListNR, BadListNC;
	QString **BadListData;
	TermEnergy **EL;
	TableLine **GL;
	double *MU, el_S, **S, SRMin, SRMax, *IsoF, q_e, q_f, *WeightF, FQS, WeightSum, minLSFQS;
	int NLRC, *PLRC, Ndev, *Row, NWeightF, *WeightFRow;
	double **L, *LRC, *dev, *RWErr;
	double Uinf;
	double iO, iA, iExp, Ro;
	IsoTab *IsoT;
	double *WFS, *WWFS, lambda;
	bool reCalcSC, reCalcLC, *ES, StopFit, LevenbergMarquardtRunning, FitStopped;
	FitResult Result;
	
private:
	void destroyFitData();
	void DestroyFitData(double ****E, double *****K, int *mJ, int **mv, int NL,
	   					double **EQS, double *EO, double *DE, double *Sig, double *Err, int **QN);
	void DestroyPotentials();
	void DestroyLS();
	void calcLMatrix();
	void calcLMEQS(int NEQ, double **EQS, double *Sig, double *Diff, double **LMEQS, double hCi_cP);
	int splineFit(bool *RPoints, double &FQSv, double &FQSn, int *&mJ, int **&mv, int &Mv, int &nEs,
				  int *&NEL, TermEnergy **&EL, int &nGS, int *&nGP, int *&nGL, int **&PS, 
                  TableLine **&GL, double &minR, double &maxR, int NumWFPoints, double ****MaxErr = 0, int MaxIt = 1, double threshFakt = 0.9999);
	double TTFit(double ****E, double *****K, int NEQ, int **QN, double *Sig, double **EQS, double *EObs,
				 double *EDiff, double ****MaxErr, double *C, double &FQSv, double b = 0.0);
    double IsoMFit(double ****E, double ****MaxErr, double &FQSv, int Iso, double UFakt, double Mass, int NumWFPoints);
	void fillWantedValueEQS(int NEQ, double *Sig, double **EQS, double *DE);
	bool isStopFit();
	void calcFitResult();
	
	WantedCoeffValue *WantedValues;
	int NumWantedValues;
	double SFQS;
	bool ****SFQSU;
};

#endif
