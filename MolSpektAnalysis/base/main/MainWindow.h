//
// C++ Interface: MainWindow
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2006 - 2019
//
// Copyright: See README file that comes with this source code
//
//

#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include "constants.h"

#include <QMainWindow>
#include <QAction>
#include <QPrinter>

#include "atom.h"

class QMenu;
class QMdiArea;
class QMdiSubWindow;
class QTableWidget;
class QSignalMapper;

class MDIChild;
class Molecule;
class TableWindow;
class TermTable;
class MeasuredTermEnergies;
class DunTable;
class Potential;
class LineTable;
class Spektrum;
class SpectList;
class SpectSimulation;
class ElState;
class Transition;
class DiagWindow;
class FCFTab;
class FitData;
class IsoTab;

 
class MainWindow : public QMainWindow
{
	Q_OBJECT

	public:
		MainWindow();
		~MainWindow();
		int getNumAtoms();
		Atom *getAtom(int Index);
		Atom *getAtom(QString Filename);
		Atom *getAtomN(QString Name);
		int getNumMolecules();
		Molecule *getMolecule(int Index);
		Molecule *getMolecule(QString Name);
		Molecule *getMoleculewFM(QString FileName);
		int getNumTermTables();
		TermTable *getTermTable(int Index);
		TermTable *getTermTable(QString Filename, Molecule *molecule, ElState *State);
		int getNumMeasuredTermTables();
		MeasuredTermEnergies *getMeasuredTermTable(int Index);
		MeasuredTermEnergies *getMeasuredTermTable(QString Filename, Molecule *molecule);
		int getNumDunTables();
		DunTable *getDunTable(int Index);
		DunTable *getDunTable(QString Filename, Molecule *molecule);
		int getNumPotentials();
		Potential *getPotential(int Index);
		Potential *getPotential(QString Filename, Molecule *molecule);
		Potential *getPotential(QString Name);
		int getNumFitDataSets();
		FitData *getFitData(int Index);
		FitData *getFitData(QString Filename, Molecule *molecule);
		int getNumLineTables();
		LineTable *getLineTable(int Index);
		LineTable *getLineTable(QString Filename, Molecule *molecule);
		int getNumFCFTables();
		FCFTab *getFCFTable(int Index);
		FCFTab *getFCFTable(QString Filename, Molecule *molecule);
		int getNumSpectra();
		Spektrum *getSpectrum(int Index);
		Spektrum *getSpectrum(QString Filename);
		Atom *CreateAtom();
		Molecule *CreateMolecule();
		TermTable *CreateTermTable();
		MeasuredTermEnergies *CreateMeasuredTermTable();
		DunTable *CreateDunTable();
		FitData *CreateFitData();
		Potential *CreatePotential(int ThreadNum = -1);
		LineTable *CreateLineTable();
		LineTable *getUnassignedLineTable();
		void loadLineTables();
		void LineTableSaved(ElState *UState, ElState *LState);
		FCFTab *CreateFCFTable(Molecule *Mol = 0);
		Spektrum *CreateSpectrum();
		SpectList *CreateSpectList();
		SpectSimulation *CreateSpectSimulation();
		QTableWidget *createTable();
		TableWindow *createTableWindow();
		ElState *CreateElState();
		Transition *CreateTransition();
		DiagWindow *CreateDiagWindow(MDIChild::Type = MDIChild::SimpleDiagWindow);
		DiagWindow *getActiveDiagWindow();
		
		void setResumeComputationLogfile(QString File);
		void setDir(QString Dir, MDIChild::Type = MDIChild::SimpleMDIChild);
		QString getDir(MDIChild::Type type = MDIChild::SimpleMDIChild);
		bool checkNotOpened(QString File, QString Filter = "");
		int getNDirRep();
		void addDirRep(QString DirRep, int NC);
		void getDirRep(int &N, QString &DirRep, int &NC);
		int getNumMDIChilds();
		MDIChild *getMDIChild(int i);
		void setActive(MDIChild *Window);
		QString getChildName(QString Suggestion, MDIChild *M);
		void showFCFTable(Molecule *Mol, QString lowerState, QString upperState);
		void showWaveFuncPlot(Potential *Pot, int Iso, int J, int v);
		void showMDIChild(QWidget *C);
		void showStatusText(QString Text = "");
		void ResumeCalculation(QString Logfile);
		void setAskForQuit(MDIChild *Window);
        void MdiChildClosed(QWidget* i_closingWindow);
        void MdiChildShown(QWidget* i_shownWindow);

        inline QMap<QString, LineTable*>& getLineTableTranslations()
        {
            return LineTableTranslations;
        }
		
	protected:
		void closeEvent(QCloseEvent *event);
		void timerEvent(QTimerEvent *event);

	private slots:
		void newAtom();
		void newMolecule();
		void newTermTable();
		void newMeasuredTermTable();
		void newDunTable();
		void newPotential();
		void newFitData();
		void newLineTable();
		void newFCFTable();
		void newSpectrum();
		void newSpectList();

		void open();
		void reload();
		void writePotData();
		void writeTFGS();
		void writeExpPotFitInput();
		void write2AtInput();
		
		void exportTFDun();
		void exportPotPoints();
		void exportAsymptoticLevels();
        void exportPicture();
		void exportTableData();
		void exportWaveFunction();
        void exportLineProfile();
		void exportObservedLevelLists();
		
		void importTermTable();
		void importDunTable();
		void importPotential();
		void importFitData();
		void importLineTable();
		void importSpectrum();
		void importSpectList();
		
		void importAbInitioPotSet();
		void importMTTPotentials();
		void importDunhamAsen();
		void read2AtOutput();
		void importCoupledPotfitOutput();
		void importCoupledTermTable();
		void importCoupledWaveFunctions();
		
		void quit();
		void save();
		void saveAs();
		void print();
		
		void showIsotopologues();
		void showLevelNumbers();
		void calcFCFTable();
        void checkAllConnections();
        void shrinkAllSpectRefs();
		
		void addRow();
		void deleteRows();
		void copyRows();
		void cutRows();
		void insertRows();
		void search();
		void shiftValues();
        void setTextInCells();
		void update();
		void findWrongData();
		void addCalculatedLevels();
		void removeDataFSource();
		void removeSingleLines();
		void fixJOffsets();
		void showNumLevels();
		void showUncertaintyStats();
		void compareLevelEnergies();
		void changeSourceOffset();
		void extractDataWithComponent();
        void extractChangedData();
        void extractNewData();
		void selectDataFSource();
		void sortTabByDeviation();
		void sortTabByDevRatio();
        void sortTabByElState();
        void sortTabByElStateAndProgNr();

		void addDunTableLine();
		void calcTermEnergies();
		void UpdateDunham();
		void ImproveDunham();
		void fitBorderLine();
		void removeFineStructure();
		void removeLambdaDoubling();
		void applyAdCorr();
		void testKratzer();
		void calcRe();
		void calcY00Te();
		void showBandConstants();
		
		void calcTermEnergiesPot();
		void calcScatWaveFunc();
		void showScatNodePos();
		void fitScatNodePos();
		void autoCalcScatLengthsPotentialSet();
		void fitControl();
		void monteCarloSim();
		void createAnaPotsFromMCSSplineSeries();
		void summarizePotInfo();
		void calcSFQSMCSSeries();
		void improvePotSeries();
		void plotPotential();
		void showFQS();
		void showTeRe();
        void calcRoughBe();
		void showBadList();
		void showFitData();
		void updateFitData();
		void fitSplinePot();
		void fitAnaPot();
		void fitTangToenniesPot();
		void fitWRobustWeighting();
		void fitIsoMass();
		void showSFuncs();
		void setEnergyOffset();
		void selPotCoeff();
		void scalePotential();
		void cdConnectSR();
		void cdConnectLR();
		void calcyss();
		void fixCoefficients();
		void varyCoefficients();
		void createMLRPot();
		
		void plotShowMouseCross(bool show);
		void plotAddPotential();
		void plotClearHistory();
		void plotShowDiagFuncs(bool checked);
		void plotShowHistory(bool checked);
		void plotShowPoints(bool checked);
		void plotPotSnapShot();
		void plotShowMarkerLabels(bool checked);

		void enableShowMenu(bool enable);
		void showTermView();
		void showTexTable();
		void showTermPlot();
		void showWaveFuncPlot();
		void showSpectSimulation();
		void showFCFTable();
		void showFCFJDependency();
		void showDataPlot();
		void showResidualPlot();
		
		void Assignvs();
		void AssignFC();
		void SplitLineTable();
    	void MarkSelected();
    	void TakeOnChanges();
		void TestProgressions();
    	void ShowUpTerm();
		void ShowGSDeviations();
    	void FindBigDiff();
		void ShowWeakProgressions();
		void FindSimilarProgressions();
    	void SelectFound();
    	void ShowUpTermTable();
		void ShowCalcRelInt();
		void RemoveDoubled();
		void sortUpTermIvJ();
		void SortProg();
		void SortIJvP();
		void SortIvPJ();
		void SortFPInt();
		void SortSpectrum();
		void SortfRemDoubled();
        void SortByProgNr();
		void SetError();
		void Shiftvup();
		void Shiftvdown();
		void ShiftJup();
		void ShiftJdown();
		void Shiftvsup();
		void Shiftvsdown();
		void ShiftIso();
		void SetvssAscending();
		void Delete();
		void setvs();
		void setFC();
		void sortbyvs();
		void setPN();
		/*void newFile();
		void cut();
		void copy();
		void paste();
		void about();
		void updateMenus();
		void updateWindowMenu();
		void switchLayoutDirection();*/
		
		void sortIvJfF();
		void sortvsIvJ();
		
		void SpectrumShowAssignmentsOnTop();
        void SpectrumDisplay_Marked();
    	void SpectrumClear_Marked();
    	void SpectrumSingleSLP();
		void SpectrumMultiSLP();
		void SpectrumAutoSLP();
    	void SpectrumTest_Transition();
    	void SpectrumFind_Satellites();
		void SpectrumFindEmissionLines();
		void SpectrumAssignBands();
		void SpectrumAssignBandsByDP();
		void SpectrumFindLinesFromTable();
		void SpectrumContinueProgressions();
		void SpectrumFindProgressions();
		void SpectrumSetLaserFrequency();
    	void SpectrumChange_Settings();
    	void SpectrumShow_Found();
		void SpectrumAcceptAssignments();
    	void SpectrumPrev_Progression();
		void SpectrumNext_Progression();
    	void SpectrumFindPeaks();
        void SpectrumFitGaussianLineProfile();
		void SpectrumCut();
		void SpectrumCutAssigned();
		void SpectrumCutStrong();
		void SpectrumAdd();
		void SpectrumSetType();
		void SpectrumNormalize();
		void SpectrumShowIntensityHistogram();

        void ShowAboutWindow();
		
		void WindowActivated(QMdiSubWindow *w);

	signals:
		void quitApp();
		void MoleculesChanged();
		void SpectrumChanged(Spektrum *Spectrum);
		void newSpectrum(Spektrum *Spectrum);
		void LineTableChanged();
		void LineTablesChanged();
		void TermTableChanged();
		void TermTablesChanged();
		void FitDataChanged();
		void DunTablesChanged();
		void PotentialsChanged();
		void FCFTablesChanged();
		void SpektrumDeleted(int Index);
		
	private:
		bool CreateAnaPotSeriesFromMCSSplinePotSeries(QString MolFN, QString StateN, QString FDDir, 
					QString SPDir, QString APDir, bool improveAnaPots, bool UseSvd, bool UseLeveMarq, int MaxIt, 
					double Prec, int nStart);
		void OpenProject(QString Filename);
		void enablePotMenu(bool enable, PotentialType PotType);
		void fillSummarizePotInfoRow(int &row, TableWindow* ResultTab, Potential *CPot, QString FQS = "", QString FQS_Bad = "",
									 QString N_Bad = "", QString N_Bad_PAL = "", QString sigma = "");
        Molecule *getCurrentMolecule();
        FitData* showSelectFitDataDialog(ElState *const i_elStateToSelectFitDataFrom, const QString& i_windowTitle, const QString& i_text);
        void disableMenues();

		MDIChild *activeMDIChild();
		//MdiChild *findMdiChild(const QString &fileName);
		LineTable *activeLineTable();
		Spektrum *activeSpectrum();
		bool checkSaved();
		bool checkSaved(MDIChild *Child);
        void calcFCFDiag(int NumWFPoints, bool update);
		IsoTab *selectIso();
						
		QMdiArea *workspace;
		QSignalMapper *windowMapper;

		QMenu *fileMenu;
		QMenu *fileNewMenu;
		QAction *newAtomAct;
		QAction *newMoleculeAct;
		QAction *newTermTableAct;
		QAction *newMeasuredTermTableAct;
		QAction *newDunTableAct;
		QAction *newPotentialAct;
		QAction *newFitDataAct;
		QAction *newLineTableAct;
		QAction *newFCFTableAct;
		QAction *newSpectrumAct;
		QAction *newSpectListAct;
		QAction *openAct;
		QAction *reloadAct;
		QAction *saveAct;
		QAction *saveAsAct;
		
		QMenu *exportMenu;
		QAction *writePotDataAct;
		QAction *writeTFGSAct;
		QAction *writeExpPotFitInputAct;
		QAction *write2AtInputAct;
		QAction *exportTFDunAct;
		QAction *exportPotPointsAct;
		QAction *exportAsymptoticLevelsAct;
		QAction *exportPictureAct;
		QAction *exportTableDataAct;
		QAction *exportWaveFunctionAct;
        QAction *exportLineProfileAct;
		QAction *exportObservedLevelListsAct;
		
		QMenu *importMenu;
		QAction *importTermTableAct;
		QAction *importDunTableAct;
		QAction *importPotentialAct;
		QAction *importFitDataAct;
		QAction *importLineTableAct;
		QAction *importSpectrumAct;
		QAction *importSpectListAct;
		
		QAction *read2AtOutputAct;
		QAction *importAbInitioPotSetAct;
		QAction *importMTTPotentialsAct;
		QAction *importDunhamAsenAct;
		QAction *importCoupledPotfitOutputAct;
		QAction *importCoupledTermTableAct;
		QAction *importCoupledWaveFunctionsAct;
		
		QAction *printAct;
		QAction *exitAct;
		
		QMenu *tableMenu;
		QAction *addRowAct;
		QAction *deleteRowsAct;
		QAction *copyRowsAct;
		QAction *cutRowsAct;
		QAction *insertRowsAct;
		QAction *searchAct;
		QAction *shiftValuesAct;
        QAction *setTextInCellsAct;
		QAction *updateAct;
		QAction *findWrongDataAct;
		QAction *addCalculatedLevelsAct;
		QAction *showNumLevelsAct;
		QAction *showUncertaintyStatsAct;
		QAction *compareLevelEnergiesAct;
		QAction *changeSourceOffsetAct;
		QAction *fixJOffsetsAct;
		QAction *removeSingleLinesAct;
		QAction *removeDataFSourceAct;
		QAction *selectDataFSourceAct;
		QAction *extractDataWithComponentAct;
        QAction *extractChangedDataAct;
        QAction *extractNewDataAct;
		QAction *sortTabByDeviationAct;
		QAction *sortTabByDevRatioAct;
        QAction *sortTabByElStateAct;
        QAction *sortTabByElStateAndProgNrAct;
		
		QMenu *MoleculeMenu;
		QAction *showIsotopologuesAct;
		QAction *showLevelNumbersAct;
		QAction *calcFCFTableAct;
        QAction *checkAllConnectionsAct;
        QAction *shrinkAllSpektRefsAct;
		
		QMenu *SpectrumMenu;
		QAction *SpectrumAutoSLPAct;
		QAction *SpectrumChangeSettingsAct;
		QAction *SpectrumClearMarkedAct;
		QAction *SpectrumDisplayMarkedAct;
		QAction *SpectrumFindPeaksAct;
		QAction *SpectrumFind_SatellitesAct;
		QAction *SpectrumAssignBandsAct;
		QAction *SpectrumAssignBandsByDPAct;
		QAction *SpectrumFindLinesFromTableAct;
		QAction *SpectrumFindEmissionLinesAct;
		QAction *SpectrumContinueProgressionsAct;
		QAction *SpectrumFindProgressionsAct;
		QAction *SpectrumMultiSLPAct;
		QAction *SpectrumPrev_ProgressionAct;
		QAction *SpectrumNext_ProgressionAct;
		QAction *SpectrumSetLaserFrequencyAct;
		QAction *SpectrumShowAssignmentsOnTopAct;
		QAction *SpectrumShowFoundAct;
		QAction *SpectrumAcceptAssignmentsAct;
		QAction *SpectrumSingleSLPAct;
		QAction *SpectrumTestProgressionAct;
        QAction *SpectrumFitGaussianLineProfileAct;
		QAction *SpectrumCutAct;
		QAction *SpectrumAddAct;
		QAction *SpectrumCutAssignedAct;
		QAction *SpectrumCutStrongAct;
		QAction *SpectrumSetTypeAct;
		QAction *SpectrumNormalizeAct;
		QAction *SpectrumShowIntensityHistogramAct;
		
		QMenu *DunhamMenu;
		QAction *addDunTableLineAct;
		QAction *calcTermEnergiesAct;
		QAction *updateDunhamAct;
		QAction *improveDunhamAct;
		QAction *fitBorderLineAct;
		QAction *removeFineStructureAct;
		QAction *removeLambdaDoublingAct;
		QAction *applyAdCorrAct;
		QAction *testKratzerAct;
		QAction *calcReAct;
		QAction *calcY00TeAct;
		QAction *showBandConstantsAct;
		
		QMenu *PotentialMenu;
		QAction *calcTermEnergiesPotAct;
		QAction *calcScatWaveFuncAct;
		QAction *showScatNodePosAct;
		QAction *fitScatNodePosAct;
		QAction *autoCalcScatLengthsPotentialSetAct;
		QAction *fitControlAct;
		QAction *monteCarloSimAct;
		QAction *createAnaPotsFromMCSSplineSeriesAct;
		QAction *summarizePotInfoAct;
		QAction *calcSFQSMCSSeriesAct;
		QAction *improvePotSeriesAct;
		QAction *plotPotentialAct;
		QAction *showFQSAct;
		QAction *showTeReAct;
        QAction *calcRoughBeAct;
		QAction *showBadListAct;
		QAction *showFitDataAct;
		QAction *updateFitDataAct;
		QAction *fitSplinePotAct;
		QAction *fitAnaPotAct;
		QAction *fitTangToenniesPotAct;
		QAction *fitWRobustWeightingAct;
		QAction *fitIsoMassAct;
		QAction *showSFuncsAct;
		QAction *setEnergyOffsetAct;
		QAction *selPotCoeffAct;
		QAction *scalePotentialAct;
		QAction *cdConnectSRAct;
		QAction *cdConnectLRAct;
		QAction *calcyssAct;
		QAction *fixCoefficientsAct;
		QAction *varyCoefficientsAct;
		QAction *createMLRPotAct;
		
		QMenu *PlotMenu;
		QAction *plotShowMouseCrossAct;
		QAction *plotAddPotentialAct;
		QAction *plotClearHistoryAct;
		QAction *plotShowDiagFuncsAct;
		QAction *plotShowHistoryAct;
		QAction *plotShowPointsAct;
		QAction *plotPotSnapShotAct;
		QAction *plotShowMarkerLabelsAct;
		
		QMenu *lineTableMenu;
		QMenu *setMenu;
		QAction *setvsAct;
		QAction *setFCAct;
		QAction *SetvssAscendingAct;
		QAction *SetErrorAct;
		QAction *SetPNAct;
		QMenu *shiftMenu;
		QAction *ShiftvupAct;
		QAction *ShiftvdownAct;
		QAction *ShiftJupAct;
		QAction *ShiftJdownAct;
		QAction *ShiftvsupAct;
		QAction *ShiftvsdownAct;
		QAction *ShiftIsoAct;
		QMenu *sortMenu;
		QAction *sortUpTermIvJAct;
		QAction *SortProgAct;
		QAction *sortbyvsAct;
		QAction *SortIJvPAct;
		QAction *SortIvPJAct;
		QAction *SortFPIntAct;
		QAction *SortSpectrumAct;
		QAction *SortfRemDoubledAct;
        QAction *SortByProgNrAct;
		
		QAction *AssignvsAct;
		QAction *AssignFCAct;
		QAction *SplitLineTableAct;
    	QAction *MarkSelectedAct;
    	QAction *TakeOnChangesAct;
		QAction *TestProgressionsAct;

		QAction *DeleteAct;
		QAction *RemoveDoubledAct;

 	  	QAction *ShowUpTermTableAct;
		QAction *ShowUpTermAct;
		QAction *ShowCalcRelIntAct;
    	
		QAction *FindBigDiffAct;
		QAction *ShowGSDeviationsAct;
		QAction *ShowWeakProgressionsAct;
    	QAction *SelectFoundAct;
		QAction *FindSimilarProgressionAct;
		
		QAction *sortIvJfFAct;
		QAction *sortvsIvJAct;
		
		QString Project;
			
		QMenu *ShowMenu;
		
		Atom *atoms[MaxAtoms];
		QMenu *ShowAtoms;
		ShowAction *ShowAtom[MaxAtoms];
		QString NAAtoms[MaxAtoms];
		
		Molecule *molecules[MaxMolecules];
		QMenu *ShowMolecules;
		ShowAction *ShowMolecule[MaxMolecules];
		QString NAMolecules[MaxMolecules];
		
		TermTable *termTables[MaxTermTables];
		QMenu *ShowTermTables;
		ShowAction *ShowTermTable[MaxTermTables];
		QString NATermTables[MaxTermTables];
		
		MeasuredTermEnergies *measuredTermTables[MaxTermTables];
		QMenu *ShowMeasuredTermTables;
		ShowAction *ShowMeasuredTermTable[MaxTermTables];
		
		DunTable *dunTables[MaxDunTables];
		QMenu *ShowDunTables;
		ShowAction *ShowDunTable[MaxDunTables];
		QString NADunTables[MaxDunTables];
		
		Potential *potentials[MaxPotentials];
		QMenu *ShowPotentials;
		ShowAction *ShowPotential[MaxPotentials];
		QString NAPotentials[MaxPotentials];
		
		FitData *fitDataSets[MaxFitDataSets];
		QMenu *ShowFitDataSets;
		ShowAction *ShowFitData[MaxFitDataSets];
		QString NAFitDataSets[MaxFitDataSets];
		
		LineTable *lineTables[MaxLineTables];
		QMenu *ShowLineTables;
		ShowAction *ShowLineTable[MaxLineTables];
		QString NALineTables[MaxLineTables];
		
		FCFTab *fcfTables[MaxFCFTables];
		QMenu *ShowFCFTables;
		ShowAction *ShowFCFTable[MaxFCFTables];
		QString NAFCFTables[MaxFCFTables];
		
		Spektrum *spectra[MaxSpectra];
		QMenu *ShowSpectra;
		ShowAction *ShowSpectrum[MaxSpectra];
		QString NASpectra[MaxSpectra];
		
		SpectList *spectLists[MaxSpectra];
		QMenu *ShowSpectLists;
		ShowAction *ShowSpectList[MaxSpectra];
		
		SpectSimulation *spectSimulations[MaxSpectSimulations];
		
		QAction *showTermViewAct;
		QAction *showTermPlotAct;
		QAction *showTexTableAct;
		QAction *showFCFTableAct;
		QAction *showFCFJDependencyAct;
		QAction *showWaveFuncPlotAct;
		QAction *showSpectSimulationAct;
		QAction *showDataPlotAct;
		QAction *showResidualPlotAct;

        QAction *aboutAct;
		
		int numAtoms, numMolecules, numTermTables, numDunTables, numPotentials; 
		int numFitDataSets, numLineTables, numSpectra, numSpectLists, numSpectSimulations;
		int numMeasuredTermTables, numFCFTables, NumWindowsToAskForQuit;
		int numNAAtoms, numNAMolecules, numNATermTables, numNADunTables, numNAPotentials, numNAFitDataSets;
		int numNALineTables, numNAFCFTables, numNASpectra, MaxDirRep, numDirRep, *DirRepNC, nCR, nCC;
		QString **cuttedRows;
		QString AtomDir, MolDir, DunDir, PotDir, TermDir, FitDataDir, LineDir, FCFDir, SpectDir, Dir; 
		QString *DirRep, SpectSimDir, PrintFile, PictureDir, TableDataDir, ResumeComputationLogfile;
		MDIChild **windowsToAskForQuit;

        QMap<QString, LineTable*> LineTableTranslations;
		//QPrinter Printer;
};

#endif
