//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//

#ifndef SELECTLINETABLEDIALOG_H
#define SELECTLINETABLEDIALOG_H


#include <QDialog>

class QComboBox;
class QPushButton;

class MainWindow;
class Molecule;
class Transition;
class LineTable;


class SelectLineTableDialog : public QDialog
{
	Q_OBJECT
	
public:
	SelectLineTableDialog(MainWindow *MW = 0);
	LineTable *getLineTable();
	
private slots:
	void MoleculeChanged(int n);
	void TransitionChanged(int n);
	
private:
	Molecule *Mol;
	Transition *Tr;
	MainWindow *MW;
	
	QComboBox *MolB, *TranB, *LTB;
	QPushButton *OK, *Cancel;
};

#endif
