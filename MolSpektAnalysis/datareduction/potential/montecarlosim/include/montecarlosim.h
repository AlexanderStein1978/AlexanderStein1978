//
// C++ Interface: MonteCarloSim
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2011 - 2019
//
// Copyright: See README file that comes with this source code
//
//

#ifndef MONTECARLOSIM_H
#define MONTECARLOSIM_H

#include <QPlainTextEdit>
#include <QDialog>
#include <QTableWidget>
#include <QProcess>

#include "tablewindow.h"
#include "includepot.h"
#include "processview.h"

struct FitResult;

class Potential;
class FitData;
class MainWindow;
class MCFSettingsDialog;
class FieldWindow;

class QLineEdit;
class QRadioButton;
class QPushButton;
class QFile;
class QTextStream;
class QMessageBox;
class QComboBox;
class QLabel;
class QTimer;
class QTcpServer;
class QTcpSocket;
class QCheckBox;
class QTime;


class MonteCarloSim : public TableWindow
{
	Q_OBJECT
	
public:
    MonteCarloSim(MainWindow *MW, MCFSettingsDialog *Dialog);
    ~MonteCarloSim();
	
protected:
	void closeEvent(QCloseEvent *Event);
	void timerEvent(QTimerEvent *Event);
	
private slots:
	void ProcessFinished(int N, int LP = -1, bool RecursiveCalls = true, FitResult *fitresult = 0);
	void FitTerminated(int N);
	void NParFitsChanged();
	void newConnection();
	void connected();
	void disconnected();
	void readyRead();
	
	inline void FitFinished(int Thread, double)
	{
		ProcessFinished(Thread, -1);
	}
	
private:
	
	void LoadPotential(Potential *Pot, int ResultNum, int ProcessNum);
	void LoadPotential(Potential *Pot, QString FileName, int ProcessNum);
	void CreatePotential(int N, Potential *Pot = 0);
	
	QFile *LogFile;
	QTextStream *LogStream;
	ProcessView **Processes;
	FitResult **BestResults, *Results;
	Potential **Pots;
	int NResults, Task, RAL, CMin, CMax, BorderReached, MaxResult, MinResult, NIpot, NFitsRunning, MaxParFits;
	int NProcesses, NLocalProcesses, *ClientLocalN, NParFits, NFixedCoefficients, *FixedCoefficients;
	int Method, Parameter, Parameter2, Status, MaxBestResults, NPos1, NPos2, *Pos1, *Pos2, ***Field;
	double Prec, RStart, RStart2, RStop, Step, Step2, Value, **Values, *IPotValue1, *IPotValue2;
	bool NoFit, ***Tried, isClient, *Started, internalFitRoutine, Improvement;
	QMessageBox *MBox;
	QString *InitPots, *ClientIP, NameB;
	FieldWindow *fieldWindow;
	QTcpServer *Server;
	QTcpSocket *Socket;
	QStringList MessagesToSend;
	QTime *Time;
	Potential *StartPot;
	QList<IncludePot> StartPotToInclude;
	QList<IncludePot>::const_iterator StartPotIt;
};

#endif // MONTECARLOSIM_H
