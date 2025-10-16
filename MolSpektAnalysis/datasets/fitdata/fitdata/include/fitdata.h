//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//


#ifndef FITDATA_H
#define FITDATA_H


#include "tablewindow.h"
#include "levelcomb.h"
#include "fitdatasortfunctions.h"
#include "fitdatacore.h"

#include <QDialog>
#include <QCheckBox>


struct TLRef;
struct Assign_vs_CompTrans;

class TableLine;
class LineTable;
class TermEnergy;
class ResidualFit;
class Spektrum;

class QListWidget;


enum {eLevel, fLevel};


class FitData : public TableWindow
{
	Q_OBJECT
	
	public:
		FitData(ElState *State = 0, MainWindow* MW = 0, Molecule* M = 0);
		~FitData();
		
        void getData(TableLine *&Lines, int &NLines, int JD = -1, int F = -2, int v = -2, int mJ = 0, int Iso = -1, ElState* state = 0);
		void getData(TableLine *&Lines, int &NLines, int *&RowN,
                     bool sortFuncs(const FitDataCore *const Tab, const int n, const int m), int *mv = 0, int mJ = 0);
		void setData(TableLine *Lines, int NLines);
		void addData(TableLine *Lines, int NLines);
		ElState *getElState() override;
		void setElState(ElState *State);
        void setMolecule(Molecule *Mol) override;
		int getMaxJ();
		int getMaxv();
		void updateData();
		void updateRow(TableLine *Line);
		void updateEnergy(int N, int *r, double *Energy);
		void updateLevels(int N, int *r, double *Energy, int *vs, int *vss);
		void setDev(double *dev, double *DevR);
		void setDev(double *dev, int *RowN, int N);
        void setDev(TableLine** TL, TLRef *SortArray, int NE, double tol);
		void setRWErr(double *RWErr);
		void setWeightFactors(int N, int *Rows, double *Factors);
		void setMCSEnergies(double *MCSE, int *Rows, int N);
		void setUncertainty(int *Rows, double *Uncertaintys, int N);
		void findWrongData();
		void deleteRows(int *rows, int N);
		void AddRow() override;
		void DeleteRows() override;
		void removeDataFSource();
		void removeSingleLines();
		bool isDataAvailable();
		bool readData(QString Filename) override;
		bool writeData(QString Filename = "") override;
		bool writeExPotFitInput(QString fileName);
		bool writeTFGS(QString Filename);
        void extractNewData(const FitData* const i_fitDataToCompare, FitData* const i_FitDataToAddNewData) const;
        void extractChangedData(const FitData* const i_fitDataToCompare, FitData* const i_FitDataToAddChangedData) const;
		int addMarkedLevel(TermEnergy &TE, Spektrum *Source);
		void removeMarkedLevel(TermEnergy &TE, Spektrum *Source);
		void clearMarkedLevels();
		void RemoveDoubled() override;
        void Assign_v(double ****TE, int NC, int NI, int NJ, int *Nv, int *IsoTrans, Assign_vs_CompTrans *CompTrans, double Tolerancy,
					  double DoublAssTol);
        void getNumLevels(int **&LNum, int &NumIso, int &NumComp, int type = 0, int maxv = -1, bool ef = true, bool disComp = true,
                          bool disElSt = true, QList<ElState *> *statesList = 0);
        void getUncertaintyStats(QList<double> &Uncert, QList<int> &Numbers);
		double getUncertaintyOfvibLevel(int v, int I, Spektrum *Source);
		int getNumProgressions(int *mv = 0, int mJ = 0);
		int getNumLines(int *mv, int mJ);
		int getNumLines(int JD, int F=-2, int v = -2, int mJ = -1, int Iso = -1); 
		void getav(int &nv, int *&v);
		QList<int> getaFC();
		void getSourceOffset(QStringList &Names, double *&Offsets);
		void setSourceOffset(QStringList &Names, double *Offsets);
		QStringList getSources();
		void selectDataFSource(QString Source);
		void setUncertainty(double Uncertainty, bool Min);
		void setFC(int FC);
        int getNumLines() const override;
		void getavIso(int *&Iso, int &NIso);
		QList<LevelComb> getAvLevelComb(int iso, int comp, int vmin, int vmax);
        bool checkSourceConnections();
        bool containsState(ElState *State);
        void sortByLTabAndProg();
        void setResidualFit(ResidualFit *i_residualFit);
        ResidualFit* getResidualFit(ElState* const i_state, const int Iso, const int v, const int comp);

        inline bool checkAllConnections()
        {
            bool result = true;
            if (!checkSourceConnections()) result = false;
            if (!TableWindow::checkAllConnections(FitDataCore::fdcFile)) result = false;
            return result;
        }

        void shrinkAllSpectRefs()
        {
            TableWindow::shrinkAllSpectRefs(FitDataCore::fdcFile);
        }
		
		inline void sortIvJF()
		{
			sortTab(heapSort(sortIvJFreqF));
		}
		
		inline void sortvsIvJ()
		{
			sortTab(heapSort(SortvsIvJ));
		}
		
		inline void sortProg()
		{
			sortTab(heapSort(sortByProg));
		}
		
		inline void sortByDev()
		{
			sortTab(heapSort(sortbyDeviation));
		}
		
		inline void sortByDevRatio()
		{
			sortTab(heapSort(sortbyDevR));
		}

        inline void sortByLineElState()
        {
            if (Tab->columnCount() > FitDataCore::fdcLineElState)
            {
                sortTab(heapSort(sortByElState));
            }
        }

        inline bool containsDataForMoreThanOneState()
        {
            return (Tab->columnCount() > FitDataCore::fdcLineElState && LineElStates != 0);
        }

    protected:

        virtual bool ReadSpecialPart(QTextStream& i_stream, const QString& i_startString) override;
		bool readData(QTextStream& S) override;
		void writeData(QTextStream& S) override;

	signals:
		void AssignmentsAccepted(FitData*);
		
	private:
        typedef enum sortForExtractNewOrChangedO{SFENOCisSmaller, SFENOCenergyIsSmaller, SFENOCisEqual, SFENOCenergyIsBigger, SFENOCisBigger} sortForExtractNewOrChangedOrder;

        void getData(TableLine *Lines, int *SA, int i_SAL, int& NLines, int *RowN = 0, int mv = -1, int *Mv = 0,
                     int mJ = 0, int JD = -1, int F = -2, int Iso = -1, bool OnlyAssignedVss = false, ElState* state = 0);
		void sortTab(int *SArray) override;
        void copyDataFromTable(const int i_numLines, int* const i_Lines, const FitData* const i_fitDataToCopyFrom);
        void prepareForExtractNewOrChanged(const FitData *const i_fitDataOld, const int i_NRnew, const int i_NROld, const bool i_withSources, const bool i_subtractSourceOffsets, int *const io_Onew, int *const io_Oold) const;
        sortForExtractNewOrChangedOrder getForExtractNewOrChangedOrder(const FitData *const i_fitDataOld, const int i_RowNew, const int i_RowOld, const bool i_withSources, const bool i_subtractSourceOffset) const;
        bool AreSourcesAvailable() const;
		int *heapSort(bool sortFuncs(const FitDataCore *const, const int, const int)) const;
		
        ElState *State, **LineElStates;
		FitDataCore* fitDataCore;
		LineTable **Sources;
		int NMarkedLevels, *FC, lRow, NSourceOffset;
		QString *SourceOffsetNames;
		double *SourceOffset;
        QList<ResidualFit*> residualFits;
};

#endif // FITDATA_H
