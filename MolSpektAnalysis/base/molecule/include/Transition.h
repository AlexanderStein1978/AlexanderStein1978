//
// C++ Interface: Transition
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2006 - 2020
//
// Copyright: See README file that comes with this source code
//
//

#ifndef TRANSITION_H
#define TRANSITION_H


#include "mdichild.h"
#include "Table.h"

#include <QLineEdit>

class ComboBox;
class Molecule;
class ElState;
class LineTable;
class FCFTab;

class QPushButton;


class Transition : public MDIChild
{
	Q_OBJECT
			
	public:
		Transition();
		~Transition();
		void setMW(MainWindow *mw, Molecule *Mol);
		void setLineBox(ComboBox *Box);
		void setFCFBox(ComboBox *Box);
		void setLSB(ComboBox *Box);
		void setUSB(ComboBox *Box);
		ComboBox *getLineBox();
		ComboBox *getFCFBox();
		ComboBox *getLSB();
		ComboBox *getUSB();
		MainWindow *getMainWindow();
		Molecule *getMolecule();
		ElState *getUpperState();
		ElState *getLowerState();
		void setUpperState(ElState *uS);
		void setLowerState(ElState *lS);
		void addLineTable(LineTable *table);
		void addLineTable(QString Name, QString FileName, QString Source);
		void addFCFTable(FCFTab *table);
		void addFCFTable(QString Name, QString FileName, QString Source);
		void setMainLineTable(LineTable *table);
		void setMainLineTable(int BoxIndex);
		bool setMainLineTable(QString FileName);
		void setMainFCFTable(FCFTab *table);
		void setMainFCFTable(int BoxIndex);
		bool setMainFCFTable(QString FileName);
		LineTable *getLineTable();
		LineTable *getLineTable(int i);
		FCFTab *getFCFTable();
		FCFTab *getFCFTable(int i);
		QString getLineTableName(int i = -1);
		QString getFCFTableName(int i = -1);
		QString getLineTableFileName(int i = -1);
		QString getFCFTableFileName(int i = -1);
		QString getLineTableSource(int i = -1);
		QString getFCFTableSource(int i = -1);
		int getNumLineTables();
		int getNumFCFTables();
		void removeLineTable(int i);
		void removeFCFTable(int i);
		void setTransNum(int N);
		void refreshStateBox();
		
		inline bool isLineTableLoaded(int i = -1)
		{
			return (i >= 0 && i < nLines ? Lines[i] != 0 
					: (nLines > 0 ? Lines[mLine] != 0 : false));
		}
		
		inline bool isFCFTableLoaded(int i = -1)
		{
			return (i >= 0 && i < nFCF ? fcf[i] != 0 
					: (nFCF > 0 ? fcf[mFCF] != 0 : false));
		}
		
		inline double getTransitionStrength()
		{
			return transStrength->text().toDouble();
		}
		
		inline void setTransitionStrengthText(QString strength)
		{
			transStrength->setText(strength);
		}
		
		inline QString getTransitionStrengthText()
		{
			return transStrength->text();
		}
		
	signals:
		void Added(int TransNum, int Col);
		void LineTableChanged();
		void FCFTableChanged();
		
	public slots:
		void refreshLineBox();
		void refreshFCFBox();
		void showMenu(QPoint P);
		bool writeData(QString FileName = "");
		
	private slots:
		void deleteLineTables(bool *delR);
		void deleteFCFTables(bool *delR);
		void updateLineTable(QString Name);
		//void updateLowerState(QString Name);
		//void updateUpperState(QString Name);
		void updateFCFTable(QString Name);
		void lineTableChanged(int row, int column);
		void fcfTableChanged(int row, int column);
		void lineTableClicked(int row, int column);
		void fcfTableClicked(int row, int column);
		void checkLineTables();
		void checkFCFTables();
		void loadLineTable(int Number);
		void loadFCFTable(int Number);
		void USBChanged(int Index);
		void uSBChanged(int Index);
		void LSBChanged(int Index);
		void lSBChanged(int Index);
		void showLineTable();
		void showFCFTable();
		void showUSMenu(QPoint P);
		void showLSMenu(QPoint P);
		void showLineMenu(QPoint P);
		void showFCFMenu(QPoint P);
		void showUState();
		void showLState();
		void LineNameChanged();
		void FCFNameChanged();
		void mergeLineTables();
		
	private:
		void addLineTable(LineTable *table, int Index);
		void addFCFTable(FCFTab *table, int Index);
        void updateName();
		
		MainWindow *MW;
		Molecule *Mol;
		ComboBox *LineB, *FCFB, *USB, *LSB, **LineTB, **FCFTB, *uSB, *lSB;
		QLineEdit *transStrength;
		QPushButton *MergeTables;
		Table *lTable, *FCFT;
		ElState *upperState, *lowerState;
		LineTable **Lines;
		FCFTab **fcf;
		int nLines, mLine, transNum, sLine, mFCF, nFCF, sFCF, nLineR, nFCFR;
		bool TBlock;
		QMenu *ShowMenu, *ShowUSMenu, *ShowLSMenu, *ShowLineMenu, *ShowFCFMenu;
		QAction *ShowAction, *ShowUSTAction, *ShowUSSAction, *ShowLSTAction, *ShowLSSAction;
		QAction *ShowLineAction, *ShowLineTAction, *ShowFCFTAction, *ShowFCFAction;
};

#endif
