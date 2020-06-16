//
// C++ Interface: PotFitLogWindow
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2007 - 2019
//
// Copyright: See README file that comes with this source code
//
//


#ifndef POTFITLOGWINDOW_H
#define POTFITLOGWINDOW_H


#include "mdichild.h"

class QPlainTextEdit;
class QPushButton;


class PotFitLogWindow : public MDIChild
{
	Q_OBJECT
	
public:
	PotFitLogWindow(QWidget *parent = 0);
	void addTextRow(QString Text);
	bool askForQuit();
	
signals:
	void CancelFit(bool wait);
	void StopFit();
	void StartFit();
	void RestartFit();
	
public slots:
	void FitStopped();
	void FitFinished();
	void FitTerminated();
	void FitRunning();
	
private slots:
	void cancelFit();
	void stopFit();

private:
	QPlainTextEdit *Log;
	QPushButton *Stop, *Cancel, *Close;
	bool fitFinished, fitStopped, fitRunning;
};

#endif
