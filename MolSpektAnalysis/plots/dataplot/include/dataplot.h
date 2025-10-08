//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//

#ifndef DATAPLOT_H
#define DATAPLOT_H

#include <DiagWindow.h>

class Molecule;
class ElState;
class LineTable;
class TermTable;
class MainWindow;

class QRect;
class QColor;
class QPainter;

class DataPlot : public DiagWindow
{
	Q_OBJECT

	public:
    	DataPlot(MainWindow *MW);
		~DataPlot();
		void PSpektrum(QPainter &P, const QRect & A, bool PrintFN);
		
	public slots:
		void moleculesChanged();
		
	private slots:
		void molBoxChanged();
		void stateBoxChanged(int Index);
		void sourceBoxChanged();
		
	private:
		void drawSymbol(QPainter *P, int x, int y, int SymbolNR);
		int heapSort(int sortFuncs(int *L1, int *L2), int **&LineData, int NLines);
		
		Molecule *mol;
		ElState *state;
		LineTable *lineT;
		TermTable *termT;
		int **Data, Nv, NJ, NPoints, SymbolSize, SymbolCount;
		QColor *SymbolColors;
		bool Block;
};

#endif
