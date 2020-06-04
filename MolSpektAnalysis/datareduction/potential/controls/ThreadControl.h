//
// C++ Interface: ThreadControl
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2007 - 2019
//
// Copyright: See README file that comes with this source code
//
//


#ifndef THREADCONTROL_H
#define THREADCONTROL_H


#include <QWidget>


class ThreadControl : public QWidget
{
	Q_OBJECT
	
public:
    ThreadControl(MainWindow* MW, ElState *St, int NumIterations, int NumParFits);
    virtual ~ThreadControl();
	
public slots:
	void CalcTerminated(int ThreadNum);
	
protected:
	virtual void CalcFinished(int ThreadNum);
	void closeEvent(QCloseEvent *event);
	void StartCalcs();
	virtual void StartCalc(int p, int n);
	virtual void AllCalculationsFinished();
	
	int NumParF, MMParF, MaxParF, NumIt, NumTermF, *ItI, *PotI, ItFinished, NFC, *FCRows;
	int NumResultTabColumns, AdTabRows;
	Potential **Pots;
	ElState *St;
	QString **PotData, ResultDir, PotName, PotFileName;
	TableWindow *ResultTab;
	MainWindow *MW;
	FitData *MFD;
	
	QProgressBar *Progress;
	QLineEdit *NumParFE, *NumItE;
	QIntValidator *NumParFEValid, *NumItEValid;
	
	
protected slots:
	void NItChanged();
	void NParFChanged();
};

#endif
