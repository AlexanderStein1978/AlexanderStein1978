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


enum {fdcIso, fdcv, fdcJ, fdcvs, fdcJs, fdcSource, fdcProg, fdcFile, fdcEnergy,
        fdcUncert, fdcObsCalc, fdcDevR, fdcLineElState};

enum {eLevel, fLevel};


class FitData : public TableWindow
{
	Q_OBJECT
	
	public:
		FitData(ElState *State = 0, MainWindow* MW = 0, Molecule* M = 0);
		~FitData();
		
        void getData(TableLine *&Lines, int &NLines, int JD = -1, int F = -2, int v = -2, int mJ = 0, int Iso = -1, ElState* state = 0);
		void getData(TableLine *&Lines, int &NLines, int *&RowN,
                     bool sortFuncs(const QTableWidget *const Tab, const int n, const int m), int *mv = 0, int mJ = 0);
		void setData(TableLine *Lines, int NLines);
		void addData(TableLine *Lines, int NLines);
		ElState *getElState();
		void setElState(ElState *State);
        void setMolecule(Molecule *Mol);
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
		void AddRow();
		void DeleteRows();
		void removeDataFSource();
		void removeSingleLines();
		bool isDataAvailable();
		bool readData(QString Filename);
		bool writeData(QString Filename = "");
		bool writeExPotFitInput(QString fileName);
		bool writeTFGS(QString Filename);
        void extractNewData(const FitData* const i_fitDataToCompare, FitData* const i_FitDataToAddNewData) const;
        void extractChangedData(const FitData* const i_fitDataToCompare, FitData* const i_FitDataToAddChangedData) const;
		int addMarkedLevel(TermEnergy &TE, Spektrum *Source);
		void removeMarkedLevel(TermEnergy &TE, Spektrum *Source);
		void clearMarkedLevels();
		void RemoveDoubled();
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
        int getNumLines() const;
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
            if (!TableWindow::checkAllConnections(fdcFile)) result = false;
            return result;
        }

        void shrinkAllSpectRefs()
        {
            TableWindow::shrinkAllSpectRefs(fdcFile);
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
            if (Tab->columnCount() > fdcLineElState)
            {
                sortTab(heapSort(sortByElState));
            }
        }

        inline bool containsDataForMoreThanOneState()
        {
            return (Tab->columnCount() > fdcLineElState && LineElStates != 0);
        }

        enum TableColumns {Iso, vss, Jss, vs, Js, Source, PN, File, WN, err, dev, devR, WeightF};

    protected:

        virtual bool ReadSpecialPart(QTextStream& i_stream, const QString& i_startString);

        virtual inline QRegExp GetStartSpecialPartRegExp() const
        {
            return QRegExp("SourceOffsets:|Begin ResidualFit");
        }
		
	signals:
		void AssignmentsAccepted(FitData*);
		
	private:
        typedef enum sortForExtractNewOrChangedO{SFENOCisSmaller, SFENOCenergyIsSmaller, SFENOCisEqual, SFENOCenergyIsBigger, SFENOCisBigger} sortForExtractNewOrChangedOrder;

        void getData(TableLine *Lines, int *SA, int i_SAL, int& NLines, int *RowN = 0, int mv = -1, int *Mv = 0,
                     int mJ = 0, int JD = -1, int F = -2, int Iso = -1, bool OnlyAssignedVss = false, ElState* state = 0);
		void sortTab(int *SArray);
        void copyDataFromTable(const int i_numLines, int* const i_Lines, const FitData* const i_fitDataToCopyFrom);
        void prepareForExtractNewOrChanged(const FitData *const i_fitDataOld, const int i_NRnew, const int i_NROld, const bool i_withSources, const bool i_subtractSourceOffsets, int *const io_Onew, int *const io_Oold) const;
        sortForExtractNewOrChangedOrder getForExtractNewOrChangedOrder(const FitData *const i_fitDataOld, const int i_RowNew, const int i_RowOld, const bool i_withSources, const bool i_subtractSourceOffset) const;
        bool AreSourcesAvailable() const;
		
        ElState *State, **LineElStates;
		LineTable **Sources;
		int NMarkedLevels, NSources, *FC, lRow, NSourceOffset;
		QString *SourceOffsetNames;
		double *SourceOffset;
		QPixmap *NewPix;
        QList<ResidualFit*> residualFits;
};

#endif // FITDATA_H
