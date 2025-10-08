//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//

#ifndef MEASUREDTERMENERGIES_H
#define MEASUREDTERMENERGIES_H


#include "tablewindow.h"

class MainWindow;
class Molecule;


class MeasuredTermEnergies : public TableWindow
{
	public:
		MeasuredTermEnergies(MainWindow *MW = 0);
		bool readData(QString FileName);
		bool writeData(QString Filename = "");
		double ****getData();
		void setMolecule(Molecule *mol);
		void addData(int Iso, int Comp, int nTerm, double **nData, int **Ass);
		
	private:
		double ****Data;
		int NumComp, NumIso, Maxv, numVRows, MaxJ;
};

#endif
