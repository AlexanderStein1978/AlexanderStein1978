//
// C++ Interface: elstate
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2006 - 2020
//
// Copyright: See README file that comes with this source code
//
//

#ifndef ELSTATE_H
#define ELSTATE_H 

#include <QString>
#include <QComboBox>
#include <QObject>
#include <QTableWidget>
#include <QLineEdit>

#include "mdichild.h"
#include "Combobox.h"

class TermTable;
class DunTable;
class Potential;
class FCFTable;
class MainWindow;
class Molecule;
class LineTable;
class FCFTab;
class FitData;
class Table;

class QKeyEvent;
class QLineEdit;
class QPoint;
class QAction;
class QMenu;
class QPushButton;
class QTabWidget;


class ElState : public MDIChild
{
	Q_OBJECT
			
	public:
		ElState();
		~ElState();
		void initialize(MainWindow *mw, Molecule *Mol, int S = 0, int lambda = 0, int Omega = 0, 
						int Parity = 0, int Symmetry = 0, QString Letter = "");
		void setMW(MainWindow *mw, Molecule *Mol);
		void setTermBox(ComboBox *Box);
		void setDunBox(ComboBox *Box);
		void setPotBox(ComboBox *Box);
		void setFitDataBox(ComboBox *Box);
		void setLambda(int nLambda);
		void setS(float nS);
		void setOmega(int O);
		void setParity(int P);
		void setSymmetry(int S);
		void setStateNum(int N);
		int getStateNum();
		ComboBox *getTermBox();
		ComboBox *getDunBox();
		ComboBox *getPotBox();
		ComboBox *getFitDataBox();
        Molecule *getMolecule() const;
		MainWindow *getMainWindow();
		int getLambda();
		float getS();
		int getOmega();
		int getParity();
		int getSymmetry();
        int getJStart(int Iso, int Comp) const;
		QString getName();
		void addTermTable(TermTable *table);
		void addTermTable(QString Name, QString FileName, QString Source);
		void setMainTermTable(TermTable *table);
		bool setMainTermTable(QString FileName);
		TermTable *getTermTable(bool calc = true);
		TermTable *getTermTable(int i);
		QString getTermTableName(int i = -1);
		QString getTermTableFileName(int i = -1);
		QString getTermTableSource(int i = -1);
		int getNumTermTables();
		void removeTermTable(int i);
		bool isTermTableLoaded(int i = -1);
		void addDunTable(DunTable *table);
		void addDunTable(QString Name, QString FileName, QString Source);
		void setMainDunTable(DunTable *table);
		bool setMainDunTable(QString FileName);
		DunTable *getDunTable();
		DunTable *getDunTable(int i);
		QString getDunTableName(int i = -1);
		QString getDunTableFileName(int i = -1);
		QString getDunTableSource(int i = -1);
		int getNumDunTables();
		void removeDunTable(int i);
		bool isDunTableLoaded(int i = -1);
		void addPotential(Potential *pot);
		void addPotential(QString Name, QString FileName, QString Source);
		void setMainPotential(Potential *pot);
		bool setMainPotential(QString FileName);
		Potential *getPotential();
		Potential *getPotential(int i);
		QString getPotentialName(int i = -1);
		QString getPotentialFileName(int i = -1);
		QString getPotentialSource(int i = -1);
		int getNumPotentials();
		void removePotential(int i);
		void removePotential(Potential *Pot);
		bool isPotentialLoaded(int i = -1);
		void addFitData(FitData *set);
		void addFitData(QString Name, QString FileName, QString Source);
		void setMainFitData(FitData *set);
		bool setMainFitData(QString FileName);
		FitData *getFitData();
		FitData *getFitData(int i);
		QString getFitDataName(int i=-1);
		QString getFitDataFileName(int i=-1);
		QString getFitDataSource(int i=-1);
		int getNumFitDataSets() const;
		void removeFitData(int i);
		bool isFitDataLoaded(int i = -1);
        double getBe();

		inline int getMainFitDataNum()
		{
			return mFD;
		}

        inline void setBe(const double i_Be)
        {
            m_Be = i_Be;
        }

        inline void setReading(const bool i_reading)
        {
            m_reading = i_reading;
        }
		
	public slots:
		void refreshDunBox();
		void refreshPotBox();
		void refreshTermBox();
		void refreshFitDataBox();
		void setMainPotential(int BoxIndex);
		void setMainDunTable(int BoxIndex);
		void setMainTermTable(int BoxIndex);
		void setMainFitData(int BoxIndex);
		void showShowMenu(QPoint P);
		void setName(QString NewName);
		bool writeData(QString FileName = "");
	
	signals:
		void PotentialChanged();
		void DunhamChanged();
		void TermTableChanged();
		void FitDataChanged();
		void PropertiesChanged(int StateNum);
		void Added(int r, int c);
		
	private slots:
		void updateTermTable(QString Name);
		void updateDunTable(QString Name);
		void updatePotential(QString Name);
		void updateFitData(QString Name);
		void deleteTermTables(bool *delR);
		void deleteDunTables(bool *delR);
		void deletePotentials(bool *delR);
		void deleteFitDataSets(bool *delR);
		void potTableChanged(int row, int column);
		void dunTableChanged(int row, int column);
		void termTableChanged(int row, int column);
		void fitDataTableChanged(int row, int column);
		void potTableClicked(int row, int column);
		void dunTableClicked(int row, int column);
		void termTableClicked(int row, int column);
		void fitDataTableClicked(int row, int column);
		void checkTermTables();
		void checkDunTables();
		void checkPotentials();
		void checkFitDataSets();
		void loadDunTable(int Number);
		void loadPotential(int Number);
		void loadTermTable(int Number);
		void loadFitData(int Number);
		void showTermTable();
		void showDunTable();
		void showPotential();
		void showFitData();
		void showDunMenu(QPoint P);
		void showTermMenu(QPoint P);
		void showPotMenu(QPoint P);
		void showFitDataMenu(QPoint P);
		void TypeChanged();
		void TermNameChanged();
		void DunNameChanged();
		void PotNameChanged();
		void FitDataNameChanged();
		
	private:
		void addDunTable(DunTable *table, int Index);
		void addTermTable(TermTable *table, int Index);
		void addPotential(Potential *pot, int Index);
		void addFitData(FitData *set, int Index);
		
		QMenu *ShowTermMenu;
		QAction *ShowTermAction;
		QAction *ShowTMAction;
		QMenu *ShowPotMenu;
		QAction *ShowPotAction;
		QAction *ShowPMAction;
		QMenu *ShowDunMenu;
		QAction *ShowDunAction;
		QAction *ShowDMAction;
		QMenu *ShowFitDataMenu;
		QAction *ShowFitDataAction;
		QAction *ShowFDMAction;
		QMenu *ShowMenu;
		QAction *ShowAction;
		QString Name;
		TermTable **termTable;
		DunTable **DunK;
		Potential **Pot;
		FitData **FitD;
		MainWindow *MW;
		Molecule *Mol;
		ComboBox *TermB, *DunB, *PotB, *FDB, **TermTB, **DunTB, **PotTB, **FDTB;
		QComboBox *lambdaB, *SB, *OmegaB, *ParB, *SymB;
		QLineEdit *StateName;
		QTabWidget *Register;
		Table *tTable, *dTable, *pTable, *fDTable;
		int lambda, mTerm, mDun, mPot, mFD, Par, Sym, Omega, numTerm, numDun, numPot, numFD, stateNum;
		int sPot, sDun, sTerm, sFD, nPotR, nDunR, nTermR, nFDR;
        bool TBlock, HS, m_reading;
		float S;
        double m_Be;
};

#endif
