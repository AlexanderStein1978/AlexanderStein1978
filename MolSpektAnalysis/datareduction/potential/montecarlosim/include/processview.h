//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//


#ifndef PROCESSVIEW_H
#define PROCESSVIEW_H


#include <QPlainTextEdit>
#include <QStringList>
#include <QProcess>

class QTimer;

class Potential;


class ProcessView : public QPlainTextEdit
{
	Q_OBJECT
	
public:
	ProcessView(QWidget *parent, QString Directory, int ProcessNumber);
	~ProcessView();
	void start(int FitNum);
	void stop();
	void setData(QString FitProg, QString FitData, QString FitInFile1, QString FitInFile2, 
				 QString FitOutInFile, QString FitInFileName, QString FitOutFileName, 
			     QString BadListName, QString PotTempFileName, Potential *Pot);
	void setPotential(Potential *Pot);
	Potential *getPotential(int N = -1);
	double getSigma(int N = -1);
	QString getPotFileName(int N);
	bool isPotAvailable(int N);
	QStringList getBadList(int N = -1);
	bool isRunning();
	void saveResults(int ResultNum);
	
signals:
	void finished(int ProcessNumber, int LP);
	
private slots:
	void read();
	void ProcessFinished(int exitCode, QProcess::ExitStatus ExitStatus);
	void testWorking();
		
private:
	void finish();
	
	int Status, ProcessNum, FitNum;
	double lSigma, oSigma;
	bool working, killed;
	QProcess *Process;
	QString FitProg, FitData, FitInFile1, FitInFile2, FitOutInFile, FitInFileName;
	QString FitOutFileName, Directory, BadListName, PotTempFileName, PotFileName; 
	Potential *Pot;
	//QTextStream *LogStream;
	QTimer *Timer;
};

#endif
