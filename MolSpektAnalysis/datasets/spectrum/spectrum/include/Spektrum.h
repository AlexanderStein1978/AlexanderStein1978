//
// C++ Interface: Spektrum
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2006 - 2019
//
// Copyright: See README file that comes with this source code
//
//

#ifndef SPEKTRUM_H
#define SPEKTRUM_H

#include <qvariant.h>
#include <qvalidator.h>
#include <qpixmap.h>
#include <QLabel>
#include <QKeyEvent>
#include <QDialog>
#include "TermPlot.h"
#include "tools.h"

#include "DiagWindow.h"
#include "tablewindow.h"
#include "simulationstateproperty.h"
#include "marker.h"
#include "progressions.h"
#include "pointwiselineprofile.h"

#include <qvarlengtharray.h>


struct AProg;
struct Band;

class ElState;
class LineTable;
class Molecule;
class Transition;
class Gaussian;
class IntProg;

class QSpacerItem;
class QAction;
class Picture;
class QScrollBar;
class QLabel;
class QLineEdit;

enum SpektrumTyp {LaserInducedFluorescence, Absorption, NormalizedAbsorption, ThermalEmission};


class Spektrum : public DiagWindow
{
	Q_OBJECT
			
public:
    Spektrum(MainWindow *MW = 0);
    ~Spektrum();
	
    IntProg* SearchLongestProg(Marker *M1);
    void editFind();
    bool readData(QString DateiName);
	bool writeData(QString FileName = "");
	void setFileName(QString FileName);
    void DisplayMarkers();
	void GetMarker(int &AnzahlMarker, Marker *&marker);
    int getNPoints();
	void ShowMarker();
    void FitLine( double LP );
    void PictureClicked( QPoint * P );
    void SearchLines();
    void FoundB();
    void FoundF();
    void ShowFound();
	void FocusSelected();
	void AcceptAssignments();
    void CSLSettings();
    void ClearMarked();
    void TTransition(int NumWFPoints);
	void FindProgressions();
	void FindLinesFromTable(LineTable *LTab);
	void SearchProgression();
    void FindSat();
	void ContinueBand();
	void ContinueProgressions();
    void ShowTermPlot();
    void ShowLineTable();
    void BRauschen();
	void SetLaserFrequency(double Frequency);
	double GetLaserFrequency();
	Marker *GetLaserLine();
	void SingleSLP();
	void MultiSLP();
	void AutoSLP();
	void MSLP(IntProg **&P, int &N);
	void FindEmissionLines(Molecule *mol, int Iso, ElState *EState, double MinH, double Tol, int MLPB);
	void IsLI();
	void addMarker(bool autoaccept = false);
	void MarkProgression(IntProg &P);
	void WriteIntensDist();
	int *getIntensDist(int numPoints, double minI, double maxI);
	void getMinMaxIntensities(double &MinI, double &MaxI);
	void getRMPH(double &Rauschen, double &MinPeakHeight);
	void AssignBands();
	void assignBand();
	void assignBandByDoubletPartners();
	Marker *getMarker(double SearchEnergy, double Tol = 0.0);
	void AverageWeighted(int NSpectra, QString *SpektFiles, int NLines, int **Lines, ElState *LS, 
						 ElState *US, int Iso, int Comp, int vs, int Js);
	void cut(double *Start, double *End, int NAreas);
	void cutAssignedLines();
	void cutStrongLines();
	double getMaxR();
	void getProgProfile(ElState *UState, ElState *LState, int Iso, int vs, int Js, int Comp, 
					    int *&vss, int *&Jss, double *&I, int &N);
	bool addByBand(ElState *UState, ElState *LState, int Iso, int vs, int vss, int Comp, double T, 
				   double SearchRadius, double Resolution, double **Result);
	bool addByProgression(ElState *LState, int Iso, int Js, int Comp, Spektrum *RSpektrum, 
						  double Ri, double Ra, double Resolution, int RIso, int RJs, int RComp,
						  double **Result);
	bool addByProgFCF(ElState *UState, ElState *LState, int Iso, int vs, int Js, int Comp, double T, double SearchRadius, 
                      double Resolution, double **Result, int NumWFPoints, double uE = 0.0, double **uWF = 0, double* TStr = 0,
                      int NC = 0, double *const o_intensF = 0, const bool i_simulate = false);
	void autoAddByProgFCF(ElState *UState, ElState *LState, int Iso, int Comp, double T, 
                          double SearchRadius, double Resolution, int NumWFPoints);
	bool add(double *OffSet, double *ScaleFakt, int N, double Ri, double Ra, double Resolution, 
			 double **Result);
	void optimizeWeighting(double **OffSet, double *TE, int NOff, int NTE, double *WF);
	void normalize(Spektrum *refSpect);
	void SetShowAssignmentsOnTop(bool V);
	bool GetShowAssignmentsOnTop();
	void LineTableSaved(ElState *UState, ElState *LState);
	int getType();
	void setType(int newType);
	void setData(double **Data, int numRows);
    double FitGaussianLineProfile(int &lineIndex, const bool considerSaturation = false, const int MaxIterations = 100, const double MinImprovements = 0.01, const double MinFreq = -1.0,
                                  const double MaxFreq = -1.0);
    void RemoveFittedLine(const int i_index);
    void SubtractFittedLine(const int i_index, const bool subtract);
    void SimulateAbsorption(const std::vector<SimulationStateProperty>& i_parametersVector);
    int GetLineFitData(double *&x, double *&y, double *&Sig, const double MinFreq, const double MaxFreq);
	    
	inline double getMinH()
	{
		return MinPeakHeight;
	}

    inline int GetNumFittedLines() const
    {
        return m_fittedLineVector.size();
    }

    Gaussian* GetFittedLine(const int i_index)
    {
        return m_fittedLineVector[i_index];
    }

signals:
	void SpectrumChanged(Spektrum *This);
    void NumberOfFittedLinesChanged();
    void SelectedRangeChanged(const double newMinE, const double newMaxE);

protected:
    virtual void PSpektrum(QPainter &P, const QRect &A, bool PrintFN );
	
private slots:
    void markRegion(QRect* i_regionToMark);

private:
    int SML(IntProg &P, int &I, int &J, double &T, double &E, int veu, double ***ELU);
	bool isLI(const int &I, int &v, const int &J, double T, double &d, double ***ELU, int veu);
	bool cbSat(int &I, int &v, int &J, Marker &L, Marker **&M, int &N, Marker *&ML);
	void GetMGMH();
	void GetMA(Marker **&A, int &N);
	int Assignvs(IntProg &P);
	bool excitedLevel(IntProg &P);
	void cFQS(IntProg &P);
	void chkDoubletts(IntProg &P);
	void eNGL(IntProg &P, Marker **&AL, int NAL, double ***ELU);
	void getMarkedLines(Marker **&M, int &N);
	void getUnassignedLines(Marker **&M, int &N);
	void FokusProg(Marker **M, int &N);
	void UpdateProg(IntProg &P);
	bool chkPrep();
	double getAUT(IntProg &P);
	void rSortP(IntProg &P);
	void ASLP(QString Filename);
	void continueBand();
	void Paint();
    int PQC();
    bool addSimulation(const double i_TUpLevel, const double *const i_OffSet, const double *const i_ScaleFakt, const int i_N, const double i_Ri,
                       const double i_Ra, const double i_Resolution, double **const io_Result) const;
	
	Marker *LaserLinie, MBeginn, MStop, **GMarker;
	int LPI, DC, numLineTablesAT, Type;
	int AssignmentStatus;
	bool useIntensities;
    double ST, GMH, MinPeakHeight, Rauschen, LaserFrequency, LaserFTol;
	Progressions Prog;
	AProg *Assigned;
	Band *band;
	LineTable **LineTablesAT;
    std::vector<Gaussian*> m_fittedLineVector;
    PointwiseLineProfile m_simulationProfile;
};

#endif // SPEKTRUM_H
