//
// C++ Interface: vAssignDialog
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2008 - 2019
//
// Copyright: See README file that comes with this source code
//
//


#ifndef VASSIGNDIALOG_H
#define VASSIGNDIALOG_H


#include <QDialog>

class QListWidget;
class QComboBox;
class QPushButton;
class QLineEdit;

class ElState;
class TermTable;
class MainWindow;


class vAssignDialog : public QDialog
{
	Q_OBJECT
	
public:
	vAssignDialog(MainWindow *parent, ElState *S);
    bool getSelectedTables(TermTable **&T, int *&vMin, int *&vMax, int &N, double &Tol, double &DATol);
	
private slots:
	void addSource();
	void removeSource();
    void stateBoxChanged(int index);
	
private:
	QListWidget *List;
    QComboBox *Box, *StateBox;
    QPushButton *Add, *Remove, *OK, *Cancel;
    QLineEdit *Tolerance, *DATolerance, *vMin, *vMax;
	
	ElState *State;
};

#endif
