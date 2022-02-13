//
// C++ Interface: potential
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2007 - 2021
//
// Copyright: See README file that comes with this source code
//
//


#ifndef POTENTIAL_H
#define POTENTIAL_H


#include "PotFit.h"
#include "tools.h"
#include "tablewindow.h"
#include "montecarlosim.h"
#include "fit.h"
#include "potworker.h"


#include <QDialog>
#include <QThread>

class ElState;
class IsoTab;
class NaturalSpline;
class AnaPot;
class MTTPot;
class MLRPot;
class TangToenniesPot;
class SplinePot;
class CouplingFuncs;
class CoupledSineWaveFunc;
class PotFitLogWindow;
class TermTable;

struct SplinePoint;
struct WantedCoeffValue;
struct PotentialData;


class QListWidget;
class QLineEdit;
class QLabel;
class QPushButton;
class QRadioButton;
class QCheckBox;
class QFile;
class QTextStream;
class QMutex;
class QPlainTextEdit;
class QCloseEvent;
class QProgressBar;
class QTextStream;


class Potential : public TableWindow
{
	Q_OBJECT
    
    enum PreliminaryPotentialType{likelyASplinePotential = -2, unknown, anaPot, splinePot, morseLongRange, modTangToen, tangToennies};
    
public:
    Potential(MainWindow *MW = 0, Molecule *Mol = 0, int ThreadNum = -1);
	Potential(const Potential &C);
    ~Potential();
    bool init(QTextStream& inStream);
	Potential *scalePotential(double newRe, double newDe);
    void getEV(double ****&Ev, int &Nc, int &NI, int &Nv, int &NJ, int NumWFPoints);
    void exportAsymptoticLevels(QString FileName, int numv, int maxJ, int NumWFPoints);
	bool readPotData(QString Filename = "");
	bool readCoupledPotData(QString Filename = "");
	bool writePotData(QString Filename = "");
	bool readPointPot(QString Filename = "");
    void calcTermEnergies(TermTable *&TT, int NumWFPoints, int Mv = -1, int MJ = -1, bool show = true);
	void getWaveFunc(int I, int J, int &Nv, double *&E, double **&WF, int *&ri, int *&ra, int *&P, 
                     double *&F, double *&N, int NumWFPoints);
	void getFCF(int &Nv1, int &Nv2, double *E1, double *E2, double **WF1, double **WF2, int *ri1,
				int *ri2, int *ra1, int *ra2, int *P1, int *P2, double *F1, double *F2,
                double *N1, double *N2, double ***&FCF, double Min1 = 0.0, double Min2 = 0.0);
    void getFCF(Potential *uPot, int Iso, int Js, int Jss, int &Nvs, int &Nvss, double ***&FCF, int NumWFPoints);
	void getFastWaveFunc(int I, int J, int C, double E, double M, double h, int NP, double *WF, double &Ri, 
						 double &Ra, double &ISq, double &F, int &i, int &c, int &a);
	void getFastFCF(int I, Potential *uPot, int Js, int uC, double uE, int N, double *lE, int *lJss, 
                    double *FCF, int NumWFPoints, double **uWF = 0, double *TStr = 0, int NWFC = 0);
	double calcScatWaveFunc(bool write = true, double **nodePos = 0, double maxR = 5000.0);
	double **fitNodePos();
	bool readData(QString Filename = "");
	bool writeData(QString Filename, bool writeFitData);
	void writePoints(QString FileName, double RMin, double RMax, int NPoints, int NDigits, bool FirstRow, bool AdCorrE, int Iso, 
					 QString MCSDir);
	void setPointData(int nPoints, double *R, double *E);
	void setMTTCoefficients(double B, double alpha1, int NA, double *A, double alpha, double beta, double gamma, 
							int NLRC, int *pLRC, double *LRC, double Uinf, int NC = 0,
							double *C = 0, double Calpha = 0.0, double Cbeta = 0.0,
							int NCLRC = 0, int *pCLRC = 0, double *CLRC = 0, 
							double CInf = 0.0);
	bool setElState(ElState *nState);
	ElState *getElState();
	void setMolecule(Molecule *Mol);
    void getTexTable(int NumWFPoints, double FQS = -1.0);
	void UpdatePot(PreliminaryPotentialType type = unknown);
	double getMinR();
	double getMaxR();
    double getD0(int NumWFPoints);
	int getNFC();
	void getReDe(double &Re, double &De);
    double guessBe();
	double getMinimumE();
	double getAsymptote();
	double *getAdCorrEnergy(double Rmin, double Step, int numPoints, double mIso1, double mIso2);
    double getFQS(int NumWFPoints, bool calcDiagFuncs = false, double *StdDev = 0, double *sigma = 0, double *FQS_PAL = 0);
	void FitSplinePot(bool showResult = true, int MaxIt = 1, double threshFakt = 0.9999);
	double CreateSplineFromAnalyticalPot();
	void FitAnaPot(bool showResult = true, bool robustWeighting = false, bool UseSVD = true, double threshhold = 1e-10,
				   bool transform = false, bool ask = true, bool useLevenbergMarquardt = true, int MaxIt = -1,
				   bool adjustTEWeightFact = false, double hCi_cP = 0.0, bool WritePotFitTraceTab = false);
	double FitTangToenniesPot(bool showResult = true);
    double FitIsoMass(int Iso, double InitStepSize, double &RIsoMass, int NumWFPoints);
	void MonteCarloSimIteration(double *MCSTab, double UncFact = 1.0);
	void VaryCoefficients(bool vary, int N=0, int *Rows = 0);
	void getFixedCoefficients(int &N, int *&Rows);
	void setMinimum(double Energy);
	void setAsymptote(double Energy);
	void shiftEnergyOffset(double Energy);
	void cdConnectSR();
	void cdConnectLR(int firstC);
    void showBadList(int NumWFPoints);
    void autoCalcScatLengthsPotentialSet(int NumWFPoints);
	void setCoefficients();
	void addPoint(double x, double y);
	void movePoint(int &n, double newx, double newy);
	void removePoint(int n);
	void setCalcDiagFuncs(bool calc);
	bool getCalcDiagFuncs();
    void getDiagFuncs(double *&WaveFSum, double *&WWFS, double &RStart, double &h, int &NP, int NumWFPoints);
	FitData *getFitData();
	void showFitData();
	void updateFitData();
	void createMLRPot();
	void setCouplingData(CoupledSineWaveFunc *WF, int NCoupled, Potential **CoupledPot, int *CoupledComp);
	void setCouplingData(CoupledSineWaveFunc *WF, int NCoupled, ElState **CoupledStates, int *CoupledComp);
	CoupledSineWaveFunc *getCoupledSineWaveFunc();
	void getCoupledPotentials(int &nCoupled, Potential **&potentials, int *&components);
    void exportWaveFunction(int NumWFPoints);
	void getLRCoeffForReading(int &N, int *&pLRC, double *&LRC);
	void getLRCoeffForWriting(int &N, int *&pLRC, double *&LRC, bool **LRCFree = 0);
    void getSplinePotForWriting(int &NSplinePoints, SplinePoint *&points, int &NLRC, int *&pLRC, double *&LRC, double &iA, double &iO, 
								double &Exp);
	void getExchangeInt(double &A, double &alpha, double &gamma);
	void setExchangeInt(double A, double alpha, double gamma);
	bool setAdCorr(int NAdCorr, double *adCorr, int TAdCorr, int PAdCorr, double RIso1, double RIso2,
				   double Rm = 0.0, double b = 1e-30);
	void getAdCorrForWriting(int &NAdCorr, double *&adCorr, int &TAdCorr, int &PAdCorr, double &RIso1, double &RIso2, double &Rm, 
							 double &b);
	FitResult getFitResult();
	PotentialType getPotType();
	void cutRows(int &numRows, int &numColumns, QString **&Data);
	void insertRows(int numRows, int numColumns, QString **Data);
	void DeleteRows();
	void setWorker(PotWorker *NWorker);
	void calcFQS_SFQS(bool ****SFQSU = 0, double SFQSRad = 0.0);
    void calcyss(const bool movingPoints);
    void serialize(QTextStream& outStream) const;
	
	inline void getSFQSU(bool ****&SFQSU, int *&mJ, int **&mv)
	{
		Worker->getSFQSU(SFQSU, mJ, mv);
	}
	
	inline double *getPoints(double Rmin, double Rmax, int numPoints, int FC = 0, std::function<double(const double)> mapping = [](const double x){return x;})
	{
		return Worker->getPoints(Rmin, Rmax, numPoints, FC, mapping);
	}

    inline double *get_dVdR(const double Rmin, const double Rmax, const int numPoints, std::function<double(const double)> mapping = [](const double x){return x;}) const
    {
        return Worker->get_dVdR(Rmin, Rmax, numPoints, mapping);
    }
	
	inline double getPoint(double R, int FC = 0)
	{
		return Worker->getPoint(R, FC);
	}
	
	inline bool isAdCorrA()
	{
		return Worker->isAdCorrA();
	}
	
	inline int getNCoupledPotentials()
	{
		return NCoupled;
	}
	
	inline void setData(QString **Data, int NRows, int NCols)
	{
		TableWindow::setData(Data, NRows, NCols);
		updatePot();
	}
	
	inline bool writeData(QString Filename = "")
	{
		return writeData(Filename, true);
	}
	
	inline IsoTab *getIsoT()
	{
		return Worker->getIsoT();
	}
	
	inline void getAdCorrForReading(int &NAdCorr, double *&adCorr, int &TAdCorr, int &PAdCorr, double &RIso1, double &RIso2, double &Rm, 
							 double &b)
	{
		Worker->getAdCorrForReading(NAdCorr, adCorr, TAdCorr, PAdCorr, RIso1, RIso2, Rm, b);
	}
	
	inline void SetSplinePar(double *Par, double ****E, double *****K, int *mJ, int **mv, int Mv, int N, 
							 int NL, int **QN, double *Sig, double **EQS, double *EO, double *DE, 
                             double &minR, double &maxR, int NumWFPoints)
	{
        Worker->SetSplinePar(Par, E, K, mJ, mv, Mv, N, NL, QN, Sig, EQS, EO, DE, minR, maxR, NumWFPoints);
	}
	
	inline SplinePoint *getSplinePoints(int &numPoints)
	{
		return Worker->getSplinePoints(numPoints);
	}
	
	inline void getS(double &Rmin, double &Rmax, int &NPoints, int &NC, double **&S)
	{
		Worker->getS(Rmin, Rmax, NPoints, NC, S);
	}
	
	inline double getMinR(double E, double Prec)
	{
		return Worker->getMinR(E, Prec);
	}
	
	inline void getMaximum(double &R, double &U)
	{
		Worker->getMaximum(R, U);
	}
	
	virtual inline void getMinimum(double &R, double &U)
	{
		Worker->getMinimum(R, U);
	}
	
	inline void CRotPot(int J, int Iso, double *U, double *R, double minR, double h, int n, int c=0, int FC = 0)
	{
		Worker->CRotPot(J, Iso, U, R, minR, h, n, c, FC);
	}
	
	inline void CalcFQS()
	{
		if (molecule == 0 || Fit->isRunning()) return;
		FQSv = -1.0;
		bFQS = 1e300;
		setFitData();
		Fit->CalcFQS();
	}
	
	inline void setThreadNum(int N)
	{
		threadNum = N;
	}
	
	inline bool isExcIntAvailable()
	{
		return Worker->isExcIntAvailable();
	}
	
	inline int getLRCTabOffset()
	{
		return LRCTabOffs;
	}
	
	inline int getadCorrTabOffset()
	{
		return adCorrTabOffs;
	}
	
	inline int getExcIntATabPos()
	{
		return ExcIntATabPos;
	}
	
	inline bool isFitRunning()
	{
		return Fit->isRunning();
	}
	
	inline bool isFitStopped()
	{
		bool R = Worker->isFitStopped();
		Worker->unlock();
		return R && Fit->isFinished();
	}
	
	inline void cdConnectLR1C()
	{
		Worker->cdConnectLR1C();
	}
	
	inline void FitSplinePotToFunction(double *E, int NPoints, double minR, double maxR)
	{
		Worker->FitSplinePotToFunction(E, NPoints, minR, maxR);
	}
	
	inline double getUinf()
	{
		return Worker->getUinf();
	}
	
	inline void setWantedValues(int numWantedValues, WantedCoeffValue* Values)
	{
		Worker->setWantedValues(numWantedValues, Values);
	}
	
	inline void setWantedLRValues(int numWantedValues, WantedCoeffValue *Values)
	{
		Worker->setWantedLRValues(numWantedValues, Values);
	}
	
	inline WantedCoeffValue *getWantedValues()
	{
		return Worker->getWantedValues();
	}

    inline void SetCouplingFuncs(CouplingFuncs* couplingFuncs)
    {
        m_couplingFuncs = couplingFuncs;
    }

    inline bool Saving() const
    {
        return m_saving;
    }

    inline double getInnerConnectionRadius() const
    {
        return Worker->getInnerConnectionRadius();
    }

    inline double getOuterConnectionRadius() const
    {
        return Worker->getOuterConnectionRadius();
    }

    inline double getSplineSlope(const int p, const double A, const double B)
    {
        return Worker->getSplineSlope(p, A, B);
    }
	
signals:
	void FitFinished(int ThreadNum, double FQS);
	void FitTerminated(int ThreadNum);
	void SFQScalculated(int ThreadNum, double SFQS, double FQS);
	
private slots:
	void updatePot();
	void badListChanged(QString *Items, int N);
	void badListChanged(double aFQS, double aLambda, double mLSFQS);
	void potentialImproved(PotentialData *Data);
	void fitFinished();
	void fitTerminated();
	void UpdateTab(PotentialData *Data = 0);
	void CreateFirstPotFitTraceTabRows(FitResult Result, double FQS);
	void restartFit();
	void startFit();
	void cancelFit(bool wait);
	
	inline void stopFit()
	{
		Worker->stopFit(true);
	}
	
private:
	bool getFastWaveFuncs(int C, int I, int J, int Nv, double M, double h, int NP,
						  double **WF, double Ri, double Ra, double *WFsq, double *F,
						  int *i, int *c, int *a);
	void NCNF(double M, double E, double *Y, double *WF, int d);
	double PFakt(double M, double E);
	void calcDerivatives(TableLine *Lines, int NL, double *E, double **dyda);
	void setFitData(int maxv = 0, int maxJ = 0, int Iso = -1, bool OnlyTermEnergies = false);
	void initLogWindow();
	bool testIfSplinePot();
	double *calcSigFuncs(int numPoints, double rMin, double rMax);
						
	int fitN_C, fitN_I, fitN_N, fitN_MI, LRCTabOffs, adCorrTabOffs, ExcIntATabPos;
	double fitN_V, fitN_A, fitN_maxr, FQSv, bFQS;
	ElState *State;
	int maxSplinePoints, NFC, threadNum, NFitIt;
    bool m_saving, showFitResult, getTexTableSplineFitRunning, writePotFitTraceTab, mWasMoving;
	TableWindow *BadList;
	QPixmap *FixPix;
	FitData *fitData;
	CoupledSineWaveFunc *WaveFuncs;
	Potential **CoupledPot;
	int *CoupledComp, NCoupled;
	PotFitLogWindow *LogWindow;
	PotWorker *Worker;
	PotFit *Fit;
	QFile *DebugFile;
	QTextStream *DebugStream;
	QStringList PotFitTraceTab;
    CouplingFuncs *m_couplingFuncs;
};

#endif
