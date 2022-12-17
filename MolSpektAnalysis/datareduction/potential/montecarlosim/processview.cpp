//
// C++ Implementation: ProcessView
//
// Description: 
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2006 - 2020
//
// Copyright: See COPYING file that comes with this distribution
//
//


#include "processview.h"
#include "potential.h"

#include <QTimer>
#include <QTextStream>


ProcessView::ProcessView(QWidget* parent, QString Dir, int PN) : QPlainTextEdit(parent)
{
	Status = 0;
	Directory = Dir;
	ProcessNum = PN;
	Process = new QProcess;
	Pot = 0;
	setReadOnly(true);
	setMaximumBlockCount(1000);
	//LogStream = 0;
	working = true;
	killed = false;
	Timer = new QTimer(this);
	connect(Timer, SIGNAL(timeout()), this, SLOT(testWorking()));
	connect(Process, SIGNAL(readyReadStandardError()), this, SLOT(read()));
	connect(Process, SIGNAL(readyReadStandardOutput()), this, SLOT(read()));
	connect(Process, SIGNAL(finished(int, QProcess::ExitStatus)), 
			this, SLOT(ProcessFinished(int,QProcess::ExitStatus)));
}

ProcessView::~ProcessView()
{
	delete Process;
	if (Pot != 0) delete Pot;
}

QStringList ProcessView::getBadList(int N)
{
	int n;
	QString BLN = (N == -1 ? BadListName : (n = BadListName.lastIndexOf('.')) == -1 ? 
				BadListName + QString::number(N) : BadListName.left(n) 
				+ QString::number(N) + BadListName.right(BadListName.length() - n));
	QFile F(BLN);
	QStringList R;
	if (F.exists())
	{
		F.open(QIODevice::ReadOnly);
		QTextStream S(&F);
		while (!S.atEnd()) R << S.readLine();
	}
	return R;
}

bool ProcessView::isRunning()
{
	if (Process->state() != QProcess::NotRunning) return true;
	return false;
}

QString ProcessView::getPotFileName(int N)
{
	int n = PotFileName.lastIndexOf('.');
	return (n == -1 ? PotFileName + QString::number(N) : PotFileName.left(n) 
					 + QString::number(N) + PotFileName.right(PotFileName.length() - n));
}

bool ProcessView::isPotAvailable(int N)
{
	QFile F(getPotFileName(N));
	return F.exists();
}

void ProcessView::finish()
{
	Timer->stop();
	emit finished(ProcessNum, -1);
}

void ProcessView::saveResults(int ResultNum)
{
	int n = PotFileName.lastIndexOf('.');
	QFile F(PotFileName);
	QString FitN = QString::number(ResultNum);
	F.copy(n == -1 ? PotFileName + FitN : PotFileName.left(n) + FitN 
										 + PotFileName.right(PotFileName.length() - n));
	F.setFileName(FitOutFileName);
	F.copy((n = FitOutFileName.lastIndexOf('.')) == -1 ? FitOutFileName + FitN
			: FitOutFileName.left(n) + FitN 
			  + FitOutFileName.right(FitOutFileName.length() - n));
	F.setFileName(BadListName);
	F.copy((n = BadListName.lastIndexOf('.')) == -1 ? BadListName + FitN
			: BadListName.left(n) + FitN + BadListName.right(BadListName.length() - n));
}

Potential* ProcessView::getPotential(int N)
{
	Pot->readData(N < 0 ? PotFileName : getPotFileName(N));
	return Pot;
}

double ProcessView::getSigma(int N)
{
	int n = FitOutFileName.lastIndexOf('.');
	QFile F(N==-1 ? FitOutFileName : n == -1 ? 
			FitOutFileName + QString::number(N) : FitOutFileName.left(n) 
			+ QString::number(N) + FitOutFileName.right(FitOutFileName.length() - n));
	double R = -1.0;
	if (F.exists())
	{
		F.open(QIODevice::ReadOnly);
		QTextStream S(&F);
		QString B;
		while (!S.atEnd() && (n = B.indexOf("Sigma=")) == -1) B = S.readLine();
		if (n >= 0)
            R = (B.right(B.length() - n - 6).split(' ', Qt::SkipEmptyParts))[0].toDouble();
	}
	return R;
}

void ProcessView::read()
{
	QString Buffer = Process->readAllStandardError();
	if (!Buffer.isEmpty() && Buffer.indexOf("argument too large in Airy!") == -1 && Buffer.indexOf("left limit not sufficient for v=") == -1)
	{
		appendPlainText(Buffer);
		//*LogStream << "Error: " << Buffer << '\n';
	}
	Buffer = Process->readAllStandardOutput();
	if (!Buffer.isEmpty() && Buffer.indexOf("argument too large in Airy!") == -1 && Buffer.indexOf("left limit not sufficient for v=") == -1)
	{
		appendPlainText(Buffer);
		//*LogStream << Buffer << '\n';
	}
	//LogStream->flush();
	working = true;
}

void ProcessView::setData(QString FProg, QString FData, QString FInFile1, 
						  QString FInFile2, QString FOutInFile, QString FInFileName,
						  QString FOutFileName, QString BListName, 
						  QString PTempFileName, Potential* P)
{
	if (Directory.right(1) != "/") Directory += '/';
	QFile F(FProg);
	FitProg = Directory + FProg.right(FProg.length() - FProg.lastIndexOf('/') - 1);
	QFile F2(FitProg);
	if (F2.exists()) F2.remove();
	F.copy(FitProg);
	F.setFileName(FData);
	FitData = Directory + FData.right(FData.length() - FData.lastIndexOf('/') - 1);
	F2.setFileName(FitData);
	if (F2.exists()) F2.remove();
	F.copy(FitData);
	FitInFile1 = FInFile1;
	FitInFile2 = FInFile2;
	FitOutInFile = FOutInFile;
	FitInFileName = Directory + FInFileName;
	FitOutFileName = Directory + FOutFileName;
	BadListName = Directory + BListName;
	PotTempFileName = Directory + PTempFileName;
	if (Pot != 0) delete Pot;
	Pot = P;
	Pot->writePotData(PotFileName = Directory + Pot->getFName());
	/*QFile *LogFile = new QFile(Directory + "ProcessView.log");
	LogFile->open(QIODevice::WriteOnly);
	LogStream = new QTextStream(LogFile);
	if (Pot->writePotData(PotFileName = Directory + Pot->getFName()))
		*LogStream << "Writing potential " << PotFileName << '\n';
	else *LogStream << "Error writing potential " << PotFileName << '\n';
	LogStream->flush();*/
	Process->setWorkingDirectory(Directory);
	//printf("Potential file: %s\n", PotFileName.toAscii().data());
}

void ProcessView::setPotential(Potential* P)
{
	Pot = P;
	Pot->writePotData(PotFileName = Directory + Pot->getFName());
	/*if (Pot->writePotData(PotFileName = Directory + Pot->getFName()))
		*LogStream << "Writing potential: " << PotFileName << '\n';
	else *LogStream << "Error writing potential " << PotFileName << '\n';
	LogStream->flush();*/
}

void ProcessView::start(int FN)
{
	QFile F(FitInFile1), F2(FitInFileName);
	if (F2.exists()) F2.remove();
	F.copy(FitInFileName);
    Process->start('\"' + FitProg + '\"', QStringList());
	//*LogStream << "Process started: " << FitProg << '\n';
	//LogStream->flush();
	Timer->start(600000);
	Status = 1;
	oSigma = lSigma = 1e10;
	FitNum = FN;
	working = true;
	killed = false;
	clear();
	if (!isVisible()) show();
}

void ProcessView::stop()
{
	QFile F(FitInFileName);
	Process->kill();
	Timer->stop();
	Status = 3;
}

void ProcessView::testWorking()
{
	if (!working && Process->state() != QProcess::NotRunning)
	{
		Process->kill();
		killed = true;
		//*LogStream << "Process killed!\n";
		//LogStream->flush();
	}
	working = false;
}

void ProcessView::ProcessFinished(int, QProcess::ExitStatus ExitStatus)
{
	if (Status == 3)
	{
		//*LogStream << "Process finished, Status = 3\n";
		//LogStream->flush();
		finish();
		return;
	}
	if (ExitStatus == QProcess::CrashExit || killed)
	{
		//*LogStream << "Process crashed, Status =" << QString::number(Status) << '\n';
		//LogStream->flush();
		QFile F(PotTempFileName);
		if (F.exists())
		{
			QFile F2(PotFileName);
			F2.remove();
			F.copy(PotFileName);
			F.remove();
            Process->start('\"' + FitProg + '\"', QStringList());
		}
		else
		{
			F.setFileName(FitInFile2);
			if (Status == 1 && F.exists())
			{
				QFile F2(FitInFileName);
				F2.remove();
				F.copy(FitInFileName);
                Process->start('\"' + FitProg + '\"', QStringList());
				Status = 2;
			}
			else 
			{
				F.setFileName(FitInFileName);
				F.remove();
				F.setFileName(FitOutInFile);
				F.copy(FitInFileName);
                Process->start('\"' + FitProg + '\"', QStringList());
				Status = 3;
			}
		}
		killed = false;
	}
	else
	{
		//*LogStream << "Process finished, Status =" << QString::number(Status) << '\n';
		//LogStream->flush();
		double Sigma = getSigma();
		if (Sigma < lSigma) 
		{
            Process->start('\"' + FitProg + '\"', QStringList());
			lSigma = Sigma;
		}
		else if (Status == 1)
		{
			QFile F(FitInFile2);
			if (F.exists())
			{
				QFile F2(FitInFileName);
				F2.remove();
				F.copy(FitInFileName);
                Process->start('\"' + FitProg + '\"', QStringList());
				Status = 2;
			}
			else finish();
		}
		else if (Status == 2 && lSigma < oSigma)
		{
			Status = 1;
			QFile F(FitInFileName);
			F.remove();
			F.setFileName(FitInFile1);
			F.copy(FitInFileName);
            Process->start('\"' + FitProg + '\"', QStringList());
			oSigma = lSigma;
		}
		else finish();
	}
}
