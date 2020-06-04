//
// C++ Implementation: MonteCarloSim
//
// Description: 
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2006 - 2019
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "montecarlosim.h"
#include "potential.h"
#include "MainWindow.h"
#include "utils.h"
#include "elstate.h"

#include <math.h>

#include <QProcess>
#include <QFile>
#include <QTextStream>
#include <QPushButton>
#include <QComboBox>
#include <QLineEdit>
#include <QMessageBox>
#include <QRadioButton>
#include <QFileDialog>
#include <QGridLayout>
#include <QTimer>
#include <QPainter>
#include <QCheckBox>
//#include <QTcpServer>
//#include <QTcpSocket>
#include <QButtonGroup>
#include <QTimerEvent>


MonteCarloSim::MonteCarloSim(MainWindow *MW, MCFSettingsDialog *D) : TableWindow(FitSeriesResultTable, MW)
{
	QString WorkDir, WorkDir2, WorkDir3, ResultLogFile, FitProg, FitData, FitInFile1;
	QString FitInFile2, InitPotDir, PotentialName, FileN, ServerIP;
	QString FitOutInFile, FitInFileName, FitOutFileName, BadListName, PotTempFileName;
	QStringList OldLogFile, BL1, BL2;
	QFile Fi;
	Potential *Pot;
	double RStop2;
	int n, m, i, j, k, RP[100], *FixC;
	int NLRC, NLRCc, *pLRC, NAdCorr, NAdCorrC, TAdCorr, PAdCorr;
	double *LRC, *adCorr, RIso1AdCorr, RIso2AdCorr, AdCorr_Rm, AdCorr_b;
	double EIA, EIalpha, EIgamma, R;
	bool DA = true, Started[100];
	Status = -1;
	BorderReached = 0;
	MBox = 0;
	Value = 0.0;
	CMin = CMax = -1;
	BestResults = 0;
	NFitsRunning = 0;
	D->getResults(WorkDir, WorkDir2, WorkDir3, ResultLogFile, FitProg, FitData, FitInFile1, FitInFile2, FitOutInFile, FitInFileName, 
				  FitOutFileName, BadListName, InitPotDir, PotentialName, PotTempFileName, Method, Prec, RStart, RStart2, RStop, 
			      RStop2, Step, Step2, Parameter, Parameter2, ServerIP, internalFitRoutine, Improvement, NParFits, Pot);
	delete D;
	Results = new FitResult[RAL = 100];
	NResults = Task = 0;
	StartPot = Pot;
	if (Pot != 0) 
	{
		if (WorkDir.right(1) != "/" && WorkDir.right(1) != "\\") WorkDir += DIRSEP;
		NameB = WorkDir + Pot->getFName();
	}
	LogFile = new QFile(ResultLogFile);
	if (LogFile->exists())
	{
		LogFile->open(QIODevice::ReadOnly);
		QTextStream readStream(LogFile);
		while (!readStream.atEnd()) OldLogFile << readStream.readLine();
		LogFile->close();
	}
	LogFile->open(QIODevice::Append);
	LogStream = new QTextStream(LogFile);
	if (!internalFitRoutine) Pot = new Potential(MW);
	QDir Dir(InitPotDir);
	QStringList L = Dir.entryList(QDir::Files);
	if (!internalFitRoutine) for (n=0; (n < L.count() ? !Pot->readData(InitPotDir + DIRSEP + L[n]) : false); n++) ;
	else n=0;
	Pot->getLRCoeffForReading(NLRC, pLRC, LRC);
	delete[] pLRC;
	delete[] LRC;
	Pot->getAdCorrForReading(NAdCorr, adCorr, TAdCorr, PAdCorr, RIso1AdCorr, RIso2AdCorr, AdCorr_Rm, AdCorr_b);
	delete[] adCorr;
	if (internalFitRoutine)
	{
		Pot->getFixedCoefficients(k, FixC);
		FixedCoefficients = new int[k+2];
		for (m=0; m<k; m++) FixedCoefficients[m] = FixC[m];
		if (k>0) delete[] FixC;
		for (m=0, NFixedCoefficients = k; m < (Method == 2 ? 2 : 1); m++)
		{
			i = (m==0 ? Parameter : Parameter2);
			if (i < NLRC) FixedCoefficients[NFixedCoefficients++] = Pot->getLRCTabOffset() + i;
			else if (i < NLRC + NAdCorr) FixedCoefficients[NFixedCoefficients++] = Pot->getadCorrTabOffset() + i - NLRC;
			else if (i == NLRC + NAdCorr) FixedCoefficients[NFixedCoefficients++] = Pot->getExcIntATabPos();
		}
	}
	setMaxParFits(Method == 2 ? (internalFitRoutine ? 100 : 3) : 1);
	switch (Method)
	{
		case 1:
			k = L.count();
			if (k==0) k=1;
			IPotValue1 = new double[k];
			InitPots = new QString[k];
			for (NIpot = 0; n < k; n++)
			{
				if (NIpot > 0) if (!Pot->readData(InitPots[NIpot] = InitPotDir + DIRSEP + L[n])) continue;
				Pot->getLRCoeffForReading(NLRCc, pLRC, LRC);
				Pot->getAdCorrForReading(NAdCorrC, adCorr, TAdCorr, PAdCorr, RIso1AdCorr, RIso2AdCorr, AdCorr_Rm, AdCorr_b);
				if (NLRC != NLRCc || NAdCorr != NAdCorrC) 
				{
					delete[] pLRC;
					delete[] LRC;
					delete[] adCorr;
					continue;
				}
				if (Parameter < NLRC) IPotValue1[NIpot] = LRC[Parameter];
				else if (Parameter < NLRC + NAdCorr) 
					IPotValue1[NIpot] = adCorr[Parameter - NLRC];
				if (Parameter >= NLRC + NAdCorr)
				{
					Pot->getExchangeInt(EIA, EIalpha, EIgamma);
					switch (Parameter - NLRC - NAdCorr)
					{
						case 0:
							IPotValue1[NIpot] = EIA;
							break;
						case 1:
							IPotValue1[NIpot] = EIalpha;
							break;
						case 2:
							IPotValue1[NIpot] = EIgamma;
							break;
						case 3:
							IPotValue1[NIpot] = 0.5 * a0_Angstrom * EIalpha; 
							break;
					}
				}
				if (IPotValue1[0] < RStart || IPotValue1[0] > RStop)
				{
					EIA = fabs(IPotValue1[0] < RStart ? RStart - IPotValue1[0] : IPotValue1[0] - RStop);
					if (fabs(IPotValue1[NIpot] < RStart ? RStart - IPotValue1[NIpot] : IPotValue1[NIpot] - RStop) < EIA)
					{
						IPotValue1[0] = IPotValue1[NIpot];
						InitPots[0] = InitPots[NIpot];
					}
				}
				else if (IPotValue1[NIpot] >= RStart && IPotValue1[NIpot] <= RStop) NIpot++;
				delete[] pLRC;
				delete[] LRC;
				delete[] adCorr;
			}
			Pot->readData(InitPots[0]);
			break;
		case 0:
			k = L.count();
			if (k==0) k=1;
			IPotValue1 = new double[k];
			InitPots = new QString[k];
			for (NIpot = 0; n<k; n++)
			{
				if (NIpot > 0) if (!Pot->readData(InitPots[NIpot] = InitPotDir + DIRSEP + L[n])) continue;
				Pot->getLRCoeffForReading(NLRCc, pLRC, LRC);
				Pot->getAdCorrForReading(NAdCorrC, adCorr, TAdCorr, PAdCorr, RIso1AdCorr, RIso2AdCorr, AdCorr_Rm, AdCorr_b);
				if (NLRC != NLRCc || NAdCorr != NAdCorrC) 
				{
					delete[] pLRC;
					delete[] LRC;
					delete[] adCorr;
					continue;
				}
				if (Parameter < NLRC) IPotValue1[NIpot] = LRC[Parameter];
				else if (Parameter < NLRC + NAdCorr) 
					IPotValue1[NIpot] = adCorr[Parameter - NLRC];
				if (Parameter >= NLRC + NAdCorr)
				{
					Pot->getExchangeInt(EIA, EIalpha, EIgamma);
					switch (Parameter - NLRC - NAdCorr)
					{
						case 0:
							IPotValue1[NIpot] = EIA;
							break;
						case 1:
							IPotValue1[NIpot] = EIalpha;
							break;
						case 2:
							IPotValue1[NIpot] = EIgamma;
							break;
						case 3:
							IPotValue1[NIpot] = 0.5 * a0_Angstrom * EIalpha; 
							break;
					}
				}
				if (IPotValue1[0] > IPotValue1[1])
				{
					WorkDir2 = InitPots[0];
					InitPots[0] = InitPots[1];
					InitPots[1] = WorkDir2;
					EIA = IPotValue1[0];
					IPotValue1[0] = IPotValue1[1];
					IPotValue1[1] = EIA;
				}
				if (NIpot == 2)
				{
					if (IPotValue1[NIpot] < IPotValue1[0])
					{
						IPotValue1[0] = IPotValue1[NIpot];
						InitPots[0] = InitPots[NIpot];
					}
					else if (IPotValue1[NIpot] > IPotValue1[1])
					{
						IPotValue1[1] = IPotValue1[NIpot];
						InitPots[1] = InitPots[NIpot];
					}
				}
				else NIpot++;
				delete[] pLRC;
				delete[] LRC;
				delete[] adCorr;
			}
			Pot->readData(InitPots[0]);
			break;
		case 2:
			Values = Create(MaxParFits = NParFits, 2);
			NPos1 = int((RStop - RStart) / Step) + 1;
			NPos2 = int((RStop2 - RStart2) / Step2) + 1;
			Field = CreateInt(NPos1, NPos2, 5);
			Pos1 = new int[NParFits];
			Pos2 = new int[NParFits];
			for (m=0; m < NParFits; m++) Pos1[m] = -1;
			if (internalFitRoutine)
			{
				Pots = new Potential*[NParFits];
				Processes = 0;
			}
			else
			{
				Processes = new ProcessView*[NParFits];
				Pots = 0;
			}
			setNumParIt(NParFits);
			connect(NumParFits, SIGNAL(editingFinished()), this, SLOT(NParFitsChanged()));
			k = L.count();
			if (k==0) k=1;
			IPotValue1 = new double[k];
			IPotValue2 = new double[k];
			InitPots = new QString[k];
			for (NIpot = 0, --n; n<k; n++)
			{
				if (Pot == 0) Pot = new Potential();
				if (Pot != StartPot) 
				{
					if (!Pot->readData(InitPots[NIpot] = InitPotDir + DIRSEP + L[n])) continue;
				}
				else InitPots[NIpot] = Pot->getFileName();
				if (InitPots[NIpot] == StartPot->getFileName() && StartPot != Pot) continue;
				Pot->getLRCoeffForReading(NLRCc, pLRC, LRC);
				Pot->getAdCorrForReading(NAdCorrC, adCorr, TAdCorr, PAdCorr, RIso1AdCorr, RIso2AdCorr, AdCorr_Rm, AdCorr_b);
				if (NLRC != NLRCc || NAdCorr != NAdCorrC) 
				{
					delete[] pLRC;
					delete[] LRC;
					delete[] adCorr;
					continue;
				}
				if (Parameter < NLRC) IPotValue1[NIpot] = LRC[Parameter];
				else if (Parameter < NLRC + NAdCorr) 
					IPotValue1[NIpot] = adCorr[Parameter - NLRC];
				if (Parameter2 < NLRC) IPotValue2[NIpot] = LRC[Parameter2];
				else if (Parameter2 < NLRC + NAdCorr) 
					IPotValue2[NIpot] = adCorr[Parameter2 - NLRC];
				if (Parameter >= NLRC + NAdCorr || Parameter2 >= NLRC + NAdCorr)
				{
					Pot->getExchangeInt(EIA, EIalpha, EIgamma);
					if ((m = Parameter - NLRC - NAdCorr) >= 0)
					{
						switch (m)
						{
							case 0:
								IPotValue1[NIpot] = EIA;
								break;
							case 1:
								IPotValue1[NIpot] = EIalpha;
								break;
							case 2:
								IPotValue1[NIpot] = EIgamma;
								break;
							case 3:
								IPotValue1[NIpot] = 0.5 * a0_Angstrom * EIalpha; 
								break;
						}
					}
					if ((m = Parameter2 - NLRC - NAdCorr) >= 0)
					{
						switch (m)
						{
							case 0:
								IPotValue2[NIpot] = EIA;
								break;
							case 1:
								IPotValue2[NIpot] = EIalpha;
								break;
							case 2:
								IPotValue2[NIpot] = EIgamma;
								break;
							case 3:
								IPotValue2[NIpot] = 0.5 * a0_Angstrom * EIalpha; 
								break;
						}
					}
				}
				if (IPotValue1[NIpot] < RStart - Step || IPotValue1[NIpot] > RStop + Step 
					|| IPotValue2[NIpot] < RStart2 - Step2 || IPotValue2[NIpot] > RStop2 + Step2) continue; 
				IncludePot IP;
				for (R = RStart, IP.pos1 = 0; R < RStop && R < IPotValue1[NIpot]; R += Step, IP.pos1++) ;
				if (R != IPotValue1[NIpot]) IP.pos1 = -1;
				for (R = RStart2, IP.pos2 = 0; R < RStop2 && R < IPotValue2[NIpot]; R += Step2, IP.pos2++) ;
				if (R != IPotValue2[NIpot]) IP.pos2 = -1;
				if (IP.pos1 >= 0 && IP.pos2 >= 0)
				{
					IP.FileName = InitPots[NIpot];
					StartPotToInclude.append(IP);
					printf("StartPotToInclude[%d]=%s\n", StartPotToInclude.count() - 1, IP.FileName.toLatin1().data());
				}
				else NIpot++;
				delete[] pLRC;
				delete[] LRC;
				delete[] adCorr;
				if (Pot == StartPot) Pot = 0;
			}
			if (Pot != StartPot)
			{
				if (Pot != 0) delete Pot;
				Pot = StartPot;
			}
			Tried = CreateBool(NPos1, NPos2, NIpot);
			for (n=0; n < NPos1; n++) for (m=0; m < NPos2; m++)
			{
				for (i=0; i < NIpot; i++) Tried[n][m][i] = false;
				for (i=0; i<5; i++) Field[n][m][i] = -1;
			}
			fieldWindow = new FieldWindow(Field, NPos1, NPos2);
			MW->showMDIChild(fieldWindow);
			Pot->setFileName(PotentialName);
			if (internalFitRoutine) for (n=0; n < NParFits; n++) CreatePotential(n, Pot);
			else for (n=0; n<3; n++)
			{
				Processes[n] = new ProcessView(MW, (n==0 ? WorkDir 
													: (n==1 ? WorkDir2 : WorkDir3)), n);
				Processes[n]->setData(FitProg, FitData, FitInFile1, FitInFile2, 
									FitOutInFile, FitInFileName, FitOutFileName, 
									BadListName, PotTempFileName, new Potential(*Pot));
				connect(Processes[n], SIGNAL(finished(int, int)), 
						this, SLOT(ProcessFinished(int, int)));
				MW->showMDIChild(Processes[n]);
			}
			for (n=0; n<100; n++) 
			{
				Started[n] = false;
				RP[n] = -1;
			}
			for (n=0, m=-1; n < OldLogFile.count(); n++)
			{
				BL1 = OldLogFile[n].split(' ');
				while (BL1.count() >= 5 ? BL1[0] == "Process" : false)
				{
					i = BL1[1].left(BL1[1].length() - 1).toInt();
					if (i>=0 && i<100)
					{
						FileN = QStringList(BL1.mid(4)).join(" ");
						j = PotentialName.lastIndexOf('.');
						if (j!=-1)
						{
							j += FileN.indexOf(PotentialName.left(j));
							RP[i] = FileN.mid(j, j - FileN.lastIndexOf('.')).toInt();
						}
						else
						{
							j = FileN.indexOf(PotentialName) + PotentialName.length();
							RP[i] = FileN.right(FileN.length() - j).toInt();
						}
					}
					if (++n == OldLogFile.count()) break;
					BL1 = OldLogFile[n].split(' ');
				}
				while (BL1.count() >= 5 ? BL1[0] == "Started" : false)
				{
					i = BL1[2].toInt();
					if (i>=0 && i<100) 
					{
						Started[0] = true;
						Pos1[i] = int(round((BL1[3].toDouble() - RStart) / Step));
						Pos2[i] = int(round((BL1[4].toDouble() - RStart2) / Step2));
					}
					if (++n == OldLogFile.count()) break;
					BL1 = OldLogFile[n].split(' ');
				}
				while (BL1.count() >= 3 ? BL1[0] == "Finished" : false)
				{
					i = BL1[2].toInt();
					if (i>=0 && i<100) Started[i] = false;
					if (++n >= OldLogFile.count()) break;
					BL1 = OldLogFile[n].split(' ');
				}
				while (BL1.count() >= 5 ? BL1[0].left(8) != "NResults=" : true)
				{
					if (++n >= OldLogFile.count()) break;
					BL1 = OldLogFile[n].split(' ');
				}
				BL2.clear();
				while (BL2.count() >= 10 ? BL2[0].toDouble() == 0.0 : true)
				{
					if (BL2.count() >= 10 ? BL2[0].left(2) == "a1" : false)
					{
						Tab->setRowCount(RAL = 100);
						Tab->setColumnCount(BL2.count());
						Tab->setHorizontalHeaderLabels(BL2);
					}
					if (++n >= OldLogFile.count()) break;
					BL2 = OldLogFile[n].split('\t');
				}
				if (n < OldLogFile.count())
				{
					m = BL1[0].right(BL1[0].length() - 8).toInt();
					if (m == RAL)
					{
						FitResult *NRA = new FitResult[RAL + 100];
						for (i=0; i < RAL; i++) NRA[i] = Results[i];
						delete[] Results;
						Results = NRA;
						Tab->setRowCount(RAL += 100);
					}
					Results[m].Index = m;
					if (BL1.count() > 5) Results[m].ProcessNum = BL1[5].right(BL1[5].length() - 11).toInt();
					else if (internalFitRoutine) Results[m].ProcessNum = 0;
					else
					{
						for (i=0; (i<3 ? !Processes[n]->isPotAvailable(i) : false); i++) ;
						Results[m].ProcessNum = i;
					}
					Results[m].FQS_Bad = BL2[BL2.count() - 1].toDouble();
					Results[m].NBad = BL2[BL2.count() - 3].toInt();
					Results[m].NBadPAL = BL2[BL2.count() - 2].toInt();
					Results[m].Sigma = BL2[BL2.count() - 4].toDouble();
					for (i=0; i < BL2.count(); i++) Tab->setItem(m, i, new QTableWidgetItem(BL2[i]));
					j = int(round((BL1[3].right(BL1[3].length() - 7).toDouble() - RStart) / Step));
					k = int(round((BL1[4].right(BL1[4].length() - 7).toDouble() - RStart2) / Step2));
					if (j >= 0 && j < NPos1 && k >= 0 && k < NPos2)
					{ 
						if (Field[j][k][0] == -1 ? true : Results[m].better(Results + Field[j][k][0])) Field[j][k][0] = m;
						if ((i = RP[Results[m].ProcessNum]) == -1) 
							for (i=m-1; (i>=0 ? Results[m].ProcessNum != Results[i].ProcessNum : false); i--) ;
						if (i>=0)
						{
							if (j>0 ? Field[j-1][k][0] == i : false) 
							{
								Field[j][k][1] = i;
								*LogStream << "Field[" << QString::number(j) << "][" << QString::number(k) << "][1] = " 
										   << QString::number(i) << '\n';
							}
							else if(k>0 ? Field[j][k-1][0] == i : false)
							{
								Field[j][k][2] = i;
								*LogStream << "Field[" << QString::number(j) << "][" << QString::number(k) << "][2] = " 
										   << QString::number(i) << '\n';
							}
							else if(j < NPos1 - 1 ? Field[j+1][k][0] == i : false) 
							{
								Field[j][k][3] = i;
								*LogStream << "Field[" << QString::number(j) << "][" << QString::number(k) << "][3] = " 
										   << QString::number(i) << '\n';
							}
							else if(k < NPos2 - 1 ? Field[j][k+1][0] == i : false) 
							{
								Field[j][k][4] = i;
								*LogStream << "Field[" << QString::number(j) << "][" << QString::number(k) << "][4] = " 
										   << QString::number(i) << '\n';
							}
						}
					}
					RP[Results[m].ProcessNum] = -1;
				}
			}
			if ((NResults = m+1) > 0) fieldWindow->update(Results, NResults, Pos1, Pos2, MaxParFits);
			for (i=0; (i<100 ? !Started[i] : false); i++) ;
			if (i==100) for (m++; DA; m++)
			{
				if (internalFitRoutine) 
				{
					Fi.setFileName(NameB.left(NameB.indexOf(".")) + QString::number(m) + ".pot");
					DA = Fi.exists();
				}
				else for (n=0, DA = false; (n < 3 ? !(DA = Processes[n]->isPotAvailable(m)) : false); n++) ;
				if (DA)
				{
					NoFit = true;
					Status = 0;
					ProcessFinished(n, m, false);
				}
			}
			if (NResults == 0 && internalFitRoutine && !StartPotToInclude.isEmpty())
			{
				//FitResult startResult;
				StartPotIt = StartPotToInclude.constBegin();
			}
			else StartPotIt = StartPotToInclude.constEnd();
			for (n=0; n < NParFits; n++)
			{
				if (Started[n] && internalFitRoutine)
				{
					Fi.setFileName(NameB.left(NameB.indexOf(".")) + ".temp.pot");
					if (!Fi.exists()) Started[n] = false;
				}
				if (Started[n])
				{
					if (internalFitRoutine)
					{
						Pots[n]->readData(NameB.left(NameB.indexOf(".")) + ".temp.pot");
						Pots[n]->setFileName(NameB.left(NameB.indexOf(".")) + QString::number(NResults));
						if (NFixedCoefficients > 0) 
							Pots[n]->VaryCoefficients(false, NFixedCoefficients, FixedCoefficients);
						/*if (Pots[n]->getPotType() == SplinePotential) Pots[n]->FitSplinePot(false, 100);
						else*/
						Pots[n]->FitAnaPot(false, false, false, 1e-10, false, false);
					}
					else Processes[n]->start(NResults);
				}
				else
				{
					Pos1[n] = Pos2[n] = -1;
					NoFit = true;
					//if (m>0) Processes[n]->getPotential();
					ProcessFinished(n, -2, false);
				}
			}
			if (!internalFitRoutine) delete Pot;
			break;
	}
	if (Method != 2)
	{
		if (!internalFitRoutine)
		{
			Pot->setFileName(PotentialName);
			Processes = new ProcessView*[1];
			Processes[0] = new ProcessView(MW, WorkDir, 0);
			Processes[0]->setData(FitProg, FitData, FitInFile1, FitInFile2, FitOutInFile, 
						FitInFileName, FitOutFileName, BadListName, PotTempFileName, Pot);
			connect(Processes[0], SIGNAL(finished(int, int)), this, SLOT(ProcessFinished(int, int)));
			MW->showMDIChild(Processes[0]);
		}
		else
		{
			Pots = new Potential*[1];
			Pots[0] = Pot;
		}
		NoFit = true;
		ProcessFinished(0);
	}
}

MonteCarloSim::~MonteCarloSim()
{
	delete[] Results;
	if (MBox != 0) delete MBox;
	if (BestResults != 0) delete[] BestResults;
	if (internalFitRoutine) delete[] Pots;
	else delete[] Processes;
	delete[] IPotValue1;
	delete[] InitPots;
	if (Method == 2)
	{
		Destroy(Tried, NPos1, NPos2);
		Destroy(Field, NPos1, NPos2);
		delete[] Pos1;
		delete[] Pos2;
		Destroy(Values, MaxParFits);
		delete[] IPotValue2;
	}
	if (internalFitRoutine && NFixedCoefficients > 0) delete[] FixedCoefficients;
}

void MonteCarloSim::closeEvent(QCloseEvent* Event)
{
   if (Status >= 0 ? QMessageBox::Ok == QMessageBox::information(this, "MolSpektAnalysis",
	   "Closing this window will abort the calculations", 
	   QMessageBox::Ok | QMessageBox::Cancel, QMessageBox::Ok) : true) 
   {
	   Status = -2;
	   int n, N = (Method == 2 ? 3 : 1);
	   if (!internalFitRoutine) for (n=0; n<N; n++) Processes[n]->stop();
	   Event->accept();
   }
   else Event->ignore();
}

void MonteCarloSim::CreatePotential(int n, Potential *Pot)
{
	if (Pot == 0) Pot = Pots[n-1];
	ElState *St = Pot->getElState();
	Pots[n] = new Potential(*Pot);
	Pots[n]->setThreadNum(n);
	St->addPotential(Pots[n]);
	connect(Pots[n], SIGNAL(FitFinished(int,double)), this, SLOT(FitFinished(int,double)));
	connect(Pots[n], SIGNAL(FitTerminated(int)), this, SLOT(FitTerminated(int)));
	MW->showMDIChild(Pots[n]);
}

void MonteCarloSim::FitTerminated(int N)
{
	*LogStream << "Fit" << QString::number(N) << " got terminated!";
	if (Method != 2) QMessageBox::information(this, "MolSpektAnalysis", "The calculations ended because the fit got terminated!");
	else
	{
		if (Field[Pos1[N]][Pos2[N]][0] == -2) Field[Pos1[N]][Pos2[N]][0]++;
		Pos1[N] = Pos2[N] = -1;
		NFitsRunning--;
		setNumParIt(--NParFits);
	}
}

void MonteCarloSim::NParFitsChanged()
{
	int n = NParFits;
	if ((NParFits = getNumParIt()) > MaxParFits)
	{
		double **nValues = new double*[NParFits];
		int n, *nPos1 = new int[NParFits], *nPos2 = new int[NParFits];
		Potential **nPots = new Potential*[NParFits];
		for (n=0; n < MaxParFits; n++)
		{
			nValues[n] = Values[n];
			nPos1[n] = Pos1[n];
			nPos2[n] = Pos2[n];
			nPots[n] = Pots[n];
		}
		delete[] Pos1;
		delete[] Pos2;
		delete[] Values;
		delete[] Pots;
		Pos1 = nPos1;
		Pos2 = nPos2;
		Values = nValues;
		Pots = nPots;
		while (n < NParFits)
		{
			nPos1[n] = nPos2[n] = -1;
			CreatePotential(n);
			nValues[n++] = new double[2];
		}
		MaxParFits = NParFits;
	}
	if (n < NParFits)
	{
		for (n=0; (n < NParFits ? Pos1[n] != -1 : false); n++) ;
		if (n < NParFits)
		{
			NoFit = true;
			ProcessFinished(n, -2);
		}
	}
}

void MonteCarloSim::ProcessFinished(int N, int LP, bool RecursiveCalls, FitResult *fitResult)
{
	int n, m=0, c;
	double V;
	int NLRC, *pLRC, NAdCorr, TAdCorr, PAdCorr, NPoints = 5000;
	double *LRC, *AdCorr, RIso1AdCorr, RIso2AdCorr, AdCorr_Rm, AdCorr_b, *E = 0, minR = 0.0, maxR = 20.0;
	Potential *Pot = (internalFitRoutine ? Pots[N] : Processes[N]->getPotential(LP));
	QString Buffer;
	bool PreFit = true;
	if (LP == -1 && fitResult == 0) NFitsRunning--;
	if (LP >= 0)
	{
		Pot->getLRCoeffForReading(NLRC, pLRC, LRC);
		Pot->getAdCorrForReading(NAdCorr, AdCorr, TAdCorr, PAdCorr, RIso1AdCorr, RIso2AdCorr, AdCorr_Rm, AdCorr_b);
		if (Parameter < NLRC) Values[N][0] = LRC[Parameter];
		else if (Parameter < NLRC + NAdCorr) Values[N][0] = AdCorr[Parameter - NLRC];
		if (Parameter2 < NLRC) Values[N][1] = LRC[Parameter2];
		else if (Parameter2 < NLRC + NAdCorr) Values[N][1] = AdCorr[Parameter2 - NLRC];
		n = NLRC + NAdCorr;
		if (Parameter >= n || Parameter2 >= n)  
		{	
			double EIA, EIalpha, EIgamma;
			Pot->getExchangeInt(EIA, EIalpha, EIgamma);
			if (Parameter >= n)
			{
				switch (Parameter - n)
				{
					case 0:
						Values[N][0] = EIA;
						break;
					case 1:
						Values[N][0] = EIalpha;
						break;
					case 2:
						Values[N][0] = EIgamma;
						break;
					case 3:
						Values[N][0] = 0.5 * EIalpha * a0_Angstrom;
						break;
				}
			}
			if (Parameter2 >= n)
			{
				switch (Parameter2 - n)
				{
					case 0:
						Values[N][1] = EIA;
						break;
					case 1:
						Values[N][1] = EIalpha;
						break;
					case 2:
						Values[N][1] = EIgamma;
						break;
					case 3:
						Values[N][1] = 0.5 * EIalpha * a0_Angstrom;
						break;
				}
			}
		}
		Pos1[N] = int(round((Values[N][0] - RStart) / Step));
		Pos2[N] = int(round((Values[N][1] - RStart2) / Step2));
		delete[] pLRC;
		delete[] LRC;
		delete[] AdCorr;
	}
	if (NoFit && LP < 0) NoFit = false;
	else
	{
		if (LP == -1) 
		{
			if (internalFitRoutine) Pot->writeData(NameB.left(NameB.indexOf('.')) + QString::number(NResults) + ".pot");
			else Processes[N]->saveResults(NResults);
		}
		*LogStream << "NResults=" << QString::number(NResults) << " Method=" << QString::number(Method) 
				   << " Status=" << QString::number(Status);
		if (Method != 2) *LogStream << " Value=" << QString::number(Value) << '\n';
		else *LogStream << " Value1=" << QString::number(Values[N][0]) << " Value2=" 
					<< QString::number(Values[N][1]) << " ProcessNum=" << QString::number(N) << '\n';
		LogStream->flush();
		if (Status >= 0)
		{
			int NRows = 0, NCols = 2;
			QString **PotData = Pot->getData(NRows, NCols);
			if (NResults == 0) 
			{
				Tab->setColumnCount(NRows + 4);
				QStringList L;
				for (n=0; n < NRows; n++) 
				{
					L << PotData[n][0];
					*LogStream << PotData[n][0] << '\t';
				}
				L << "Sigma" << "NBad" << "NBad_PAL" << "FQS_Bad";
				*LogStream << "Sigma\tNBad\tNBad_PAL\tFQS_Bad\n";
				Tab->setHorizontalHeaderLabels(L);
				Tab->setRowCount(100);
			}
			else if (NResults == RAL)
			{
				FitResult *NRA = new FitResult[RAL + 100];
				for (n=0; n < RAL; n++) NRA[n] = Results[n];
				delete[] Results;
				Results = NRA;
				Tab->setRowCount(RAL += 100);
			}
			for (n=0; n < NRows; n++)
			{
				*LogStream << PotData[n][1] << '\t';
				Tab->setItem(NResults, n, new QTableWidgetItem(PotData[n][1]));
			}
			if (fitResult != 0) Results[NResults] = *fitResult;
			else if (internalFitRoutine) Results[NResults] = Pot->getFitResult();
			else Results[NResults].Sigma = Processes[N]->getSigma(LP);
			Buffer = QString::number(Results[NResults].Sigma, 'f', 5);
			*LogStream << Buffer << '\t';
			Tab->setItem(NResults, NRows, new QTableWidgetItem(Buffer));
			if (!internalFitRoutine)
			{
				QStringList BadL = Processes[N]->getBadList(LP);
				Results[NResults].FQS_Bad = 0.0;
				Results[NResults].NBad = BadL.count() - 1;
				Results[NResults].NBadPAL = 0;
				for (n=0; n < Results[NResults].NBad; n++)
				{
					if (BadL[n].mid(8, 5) == " 9999") Results[NResults].NBadPAL++;
					V = BadL[n].mid(48, 11).toDouble() / BadL[n].mid(39, 9).toDouble();
					Results[NResults].FQS_Bad += V*V;
				}
			}
			Buffer = QString::number(Results[NResults].NBad);
			Tab->setItem(NResults, NRows + 1, new QTableWidgetItem(Buffer));
			*LogStream << Buffer << '\t';
			Buffer = QString::number(Results[NResults].NBadPAL);
			Tab->setItem(NResults, NRows + 2, new QTableWidgetItem(Buffer));
			*LogStream << Buffer << '\t';
			Buffer = QString::number(Results[NResults].FQS_Bad, 'f', 5);
			Tab->setItem(NResults, NRows + 3, new QTableWidgetItem(Buffer));
			*LogStream << Buffer << '\n';
			LogStream->flush();
			Results[NResults].Index = NResults;
			Results[NResults].ProcessNum = N;
			if (LP != -2) Results[NResults++].ParValue = Value;
		}
	}
	if (LP >= 0)
	{
		if (Pos1[N] >= 0 && Pos1[N] < NPos1 && Pos2[N] >= 0 && Pos2[N] < NPos2 ? 
				(Field[Pos1[N]][Pos2[N]][0] != -1 ? Results[NResults - 1].better(Results + Field[Pos1[N]][Pos2[N]][0]) : true) : false)
			Field[Pos1[N]][Pos2[N]][0] = LP;
		return;
	}
	if (Method == 0)
	{
		switch (Status)
		{
			case -2:
				return;
			case 0:
				if (Results[NResults-1].NBad >= Results[0].NBad + 3
					|| Results[NResults-1].NBadPAL > Results[0].NBadPAL)
				{
					Status = 1;
					CMin = NResults - 2;
				}
				else
				{
					Value += Step;
					if (Value > RStop)
					{
						BorderReached = 1;
						Value = Results[0].ParValue - Step;
						Status = 2;
						Buffer = 
	"The upper maximum paramter value is reached and thus the estimated uncertainty may not be reliable!";
						*LogStream << Buffer << '\n';
						LogStream->flush();
						MBox = new QMessageBox(QMessageBox::Information, "MolSpektAnalysis", Buffer, 
											   QMessageBox::Ok);
						MBox->show();
					}
					break;
				}
			case 1:
				if (Results[NResults-1].NBad >= Results[0].NBad + 3
					|| Results[NResults-1].NBadPAL > Results[0].NBadPAL)
					CMax = NResults - 1;
				else CMin = NResults - 1;
				if (Results[CMax].ParValue - Results[CMin].ParValue <= Prec)
				{
					Status = 2;
					MaxResult = CMax;
					Value = Results[0].ParValue;
				}
				else
				{
					Value = (Results[CMin].ParValue + Results[CMax].ParValue) * 0.5;
					break;
				}
				/* Falls through. */
			case 2:
				if (Results[NResults-1].NBad >= Results[0].NBad + 3
					|| Results[NResults-1].NBadPAL > Results[0].NBadPAL)
				{
					Status = 3;
					CMax = NResults - 2;
				}
				else
				{
					Value -= Step;
					if (Value < RStart)
					{
						BorderReached += 2;
						Status = -2;
						if (BorderReached == 3)
						{
							Buffer = 
"Finished, both borders are reached, the uncertainty of the parameter is larger than the search interval.";
							*LogStream << Buffer << '\n';
							LogStream->flush();
							QMessageBox::information(this, "MolSpektAnalysis", Buffer);
						}
						else
						{
							Buffer = 
					"The lower maximum paramter value is reached and thus the estimated uncertainty of " 
							         + QString::number(Results[CMax].ParValue - Results[0].ParValue, 'g',
													   1 - int(log(Prec))) + " may not be reliable!";
							*LogStream << Buffer << '\n';
							LogStream->flush();
							QMessageBox::information(this, "MolSpektAnalysis", Buffer);
						}
						return;
					}
					break;
				}
				/* Falls through. */
			case 3:
				if (Results[NResults-1].NBad >= Results[0].NBad + 3
					|| Results[NResults-1].NBadPAL > Results[0].NBadPAL)
					CMin = NResults - 1;
				else CMax = NResults - 1;
				if (Results[CMax].ParValue - Results[CMin].ParValue <= Prec)
				{
					Status = -2;
					MinResult = CMin;
					Value = Results[0].ParValue;
					if (BorderReached == 0)
					{
						Buffer = "Calculatuions succesful finished. The estimated uncertainty is +/- " 
								+ QString::number(0.5 * (Results[CMax].ParValue - Results[CMin].ParValue),
								                  'g', 1 - int(log(Prec)))
								+ ". The result with maximum parameter value is fit number " 
								+ QString::number(CMax)
								+ " and the result with minimum parameter value fit number " 
								+ QString::number(CMin) + ".";
						*LogStream << Buffer << '\n';
						LogStream->flush();
						QMessageBox::information(this, "MolSpektAnalysis", Buffer);
					}
					else
					{
						Buffer = "Calculations finished. The estimated uncertainty of "
								+ QString::number(Results[0].ParValue - Results[CMin].ParValue, 'g',
												  1 - int(log(Prec)))
				+ " may not be reliable, because the upper border of the search interval was reached.";
						*LogStream << Buffer << '\n';
						LogStream->flush();
						QMessageBox::information(this, "MolSpektAnalysis", Buffer);
					}
					return;
				}
				else Value = (Results[CMin].ParValue + Results[CMax].ParValue) * 0.5;
				break;
		}
	}
	else if (Method == 1)
	{
		int BRIndex = 0;
		if (Status >= 0 && Status <= 3)
		{
			BRIndex = (fabs(double(c = int((Value - RStart) / Step)) * Step + RStart - Value) 
						< 0.01 * Step ? c : c+1); 
			if (BestResults[BRIndex] != 0 ? Results[NResults - 1].better(BestResults[BRIndex]) : true) 
				 BestResults[BRIndex] = Results + NResults - 1;
		}
		switch (Status)
		{
			case -2:
				return;
			case -1:
				if (BestResults != 0) delete[] BestResults;
				BestResults = new FitResult*[MaxBestResults = int((RStop - RStart) / Step) + 2];
				for (n=0; n < MaxBestResults; n++) BestResults[n] = 0;
				break;
			case 0:
				Value -= Step;
				if (Value < RStart)
				{
					if (fabs(Value + Step - RStart) < 0.01 * Step)
					{
						Value = Results[0].ParValue;
						Status = 1;
						/*FitResult FBB;
						for (n=0; 2 * n < NResults - 1; n++)
						{
							FBB = Results[n];
							Results[n] = Results[NResults - 1 - n];
							Results[NResults - 1 - n] = FBB;
						}*/
					}
					else
					{
						Value = RStart;
						break;
					}
				}
				else break;
                /* Falls through. */
			case 1:
				Value += Step;
				if (Value > RStop) 
				{
					if (fabs(Value - Step - RStop) < 0.01 * Step)
					{
						if (NResults < 3) 
						{
							Status = -2;
							Buffer = 
"Calculations finished. No minumum could be found, since the step size is too large compared to the search interval.";
							*LogStream << Buffer << '\n';
							LogStream->flush();
							QMessageBox::information(this, "MolSpektAnalysis", Buffer);
							return;
						}
						else
						{
							for (n=1, c=0; n < NResults - 1; n++) 
							{
								if (BestResults[n]->better(BestResults[n-1]) 
									&& BestResults[n]->better(BestResults[n+1])) c++;
								else if (BestResults[n]->better(BestResults[n-1]) 
									&& BestResults[n]->equal(BestResults[n+1]))
								{
									while(n < NResults - 1 ? BestResults[n]->equal(BestResults[n+1]) 
											: false) n++;
									if (n < NResults - 1 ? BestResults[n]->better(BestResults[n+1]) 
										: false) c++;
								}
							}
							if (c==1) 
							{
								Status = 4;
								NoFit = true;
								ProcessFinished(0);
								return;
							}
							else 
							{
								Status = 2;
								Value = RStart;
								if (internalFitRoutine) 
									Buffer = NameB.left(NameB.indexOf('.')) + QString::number(BestResults[0]->Index) + ".pot";
								else Buffer = Processes[N]->getPotFileName(BestResults[0]->Index);
								*LogStream << "Reading potential " << Buffer << '\n';
								LogStream->flush();
								QString FName = Pot->getFileName();
								Pot->readData(Buffer);
								if (!internalFitRoutine) Pot->setFileName(FName);
								else Pot->VaryCoefficients(false, NFixedCoefficients, FixedCoefficients);
							}
						}
					}
					else Value = RStop;
				}
				break;
			case 2:
				if (BestResults[BRIndex] == Results + NResults - 1 && RStop - Value > 0.01 * Step)
				{
					Value = BestResults[BRIndex + 1]->ParValue;
					break;
				}
				for (n = BRIndex; (n+1 < MaxBestResults ? (BestResults[n+1] != 0 ? BestResults[n+1] 
					- BestResults[n] == 1 : false)
						: false); n++) ;
				if (RStop - BestResults[n]->ParValue > 0.01 * Step) Value = BestResults[n+1]->ParValue;
				else
				{
					while (n>0 ? BestResults[n-1] - BestResults[n] == 1 : false) n--;
					if (n==0)
					{
						Status = 4;
						NoFit = true;
						ProcessFinished(0);
						return;
					}
					Value = BestResults[n-1]->ParValue;
					Status = 3;
				}
				if (BestResults[n]->Index != NResults - 1)
				{
					Buffer = (internalFitRoutine ? NameB.left(NameB.indexOf('.')) + QString::number(BestResults[n]->Index) + ".pot"
												 : Processes[N]->getPotFileName(BestResults[n]->Index));
					*LogStream << "Reading potential " << Buffer << '\n';
					LogStream->flush();
					QString FName = Pot->getFileName();
					Pot->readData(Buffer);
					if (internalFitRoutine) Pot->VaryCoefficients(false, NFixedCoefficients, FixedCoefficients);
					else Pot->setFileName(FName);
				}
				break;
			case 3:
				if (BestResults[BRIndex] == Results + NResults - 1 && Value - RStart > 0.01 * Step)
				{
					Value = BestResults[BRIndex - 1]->ParValue;
					break;
				}
				for (n = BRIndex; (n-1 >= 0 ? BestResults[n-1] - BestResults[n] == 1 : false); n--) ;
				if (n!=0) Value = BestResults[n-1]->ParValue;
				else
				{
					while (n+1 < MaxBestResults ? (BestResults[n+1] != 0 ? BestResults[n+1] 
							- BestResults[n] == 1 : false) : false) n++;
					if (n+1 < MaxBestResults ? BestResults[n+1] != 0 : false) NoFit = true;
					else
					{
						Value = BestResults[n+1]->ParValue;
						Status = 2;
					}
				}
				if (NoFit) 
				{
					NoFit = false;
					Status = 4;
				}
				else
				{
					if (BestResults[n]->Index != NResults - 1)
					{
						Buffer = (internalFitRoutine ? NameB.left(NameB.indexOf('.')) + QString::number(BestResults[n]->Index) + ".pot"
												     : Processes[N]->getPotFileName(BestResults[n]->Index));
						*LogStream << "Reading potential " << Buffer << '\n';
						LogStream->flush();
						QString FName = Pot->getFileName();
						Pot->readData(Buffer);
						if (internalFitRoutine) Pot->VaryCoefficients(false, NFixedCoefficients, FixedCoefficients);
						else Pot->setFileName(FName);
					}
					break;
				}
				/* Falls through. */
			case 4:
				if (CMin == -1)
				{
					FitResult *BestResult;
					for (BestResult = BestResults[n=c=0]; (n < MaxBestResults ? BestResults[n] != 0 
							: false); n++)
						if (BestResults[n]->better(BestResult))
					{
						BestResult = BestResults[n];
						c=n;
					}
					if (c==0 || (c+1 < MaxBestResults ? BestResults[c+1] == 0 : true))
					{
						Buffer = 
			"Calculations unsuccessfully finished, because the selected parameter range is too small!";
						*LogStream << Buffer << '\n';
						LogStream->flush();
						QMessageBox::information(this, "MolSpektAnalysis", Buffer);
						return;
					}
					for (n=1; (c+n < MaxBestResults ? (BestResults[c+n] != 0 ? 
							BestResults[c]->equal(BestResults[c+n]) : false) : false); n++) ;
					if (n>1)
					{
						c += n/2;
						if (2 * (n/2) != n)
						{
							Buffer = "Calculations succesfully finished. The estimated paramter value is " 
								+ QString::number(BestResults[c]->ParValue, 'g', 1 - int(log(Prec))) + '.';
							*LogStream << Buffer << '\n';
							LogStream->flush();
							QMessageBox::information(this, "MolSpektAnalysis", Buffer);
							return;
						}
						CMin = BestResults[c]->Index;
						CMax = BestResults[c+1]->Index;
						Value = 0.5 * (BestResults[c]->ParValue + BestResults[CMax]->ParValue);
						if (CMin < NResults - 1 && CMax < NResults - 1)
						{
							Buffer = (internalFitRoutine ? 
										NameB.left(NameB.indexOf('.')) + QString::number(BestResults[n]->Index) + ".pot"
										 : Processes[N]->getPotFileName(BestResults[n]->Index));
							*LogStream << "Reading potential " << Buffer << '\n';
							LogStream->flush();
							QString FName = Pot->getFileName();
							Pot->readData(Buffer);
							if (internalFitRoutine) Pot->VaryCoefficients(false, NFixedCoefficients, FixedCoefficients);
							else Pot->setFileName(FName);
						}
						break;
					}
					CMin = BestResults[c-1]->Index;
					CMax = BestResults[c+1]->Index;
				}
				else c = NResults - 1;
				if (Results[CMin].equal(Results + c) && Results[c].equal(Results + CMax))
				{
					Buffer = "Calculations succesfully finished. The estimated paramter value is " 
							+ QString::number(BestResults[c]->ParValue, 'g', 1 - int(log(Prec))) + '.';
					*LogStream << Buffer << '\n';
					LogStream->flush();
					QMessageBox::information(this, "MolSpektAnalysis", Buffer);
					return;
				}
				if (Results[c].better(Results + CMin) || Results[c].better(Results + CMax))
				{
					if (Results[CMin].NBadPAL == Results[c].NBadPAL 
						&& Results[c].NBadPAL == Results[CMax].NBadPAL
						&& Results[CMin].NBad == Results[c].NBad
						&& Results[c].NBad == Results[CMax].NBad)
					{
						if (Results[c].FQS_Bad < Results[CMin].FQS_Bad 
								|| Results[c].FQS_Bad < Results[CMax].FQS_Bad)
							ParabInterpol(Results[CMin].ParValue, Results[CMin].FQS_Bad,
											Results[c].ParValue, Results[c].FQS_Bad,
											Results[CMax].ParValue, Results[CMax].FQS_Bad,
											Value, V);
							
						else ParabInterpol(Results[CMin].ParValue, Results[CMin].Sigma,
											Results[c].ParValue, Results[c].Sigma,
											Results[CMax].ParValue, Results[CMax].Sigma, Value, V);
						if (fabs((Value - Results[c].ParValue) / Value) < Prec) 
						{
							Status = 5;
							break;
						}
						if (Value > Results[CMin].ParValue && Value < Results[CMax].ParValue)
						{
							if (Value < Results[c].ParValue) CMax = c;
							else CMin = c;
							break;
						}
					}
					else 
					{
						if (Results[CMin].NBadPAL < Results[CMax].NBadPAL
							|| (Results[CMin].NBadPAL == Results[CMax].NBadPAL
							&& (Results[CMin].NBad < Results[CMax].NBad
							|| (Results[CMin].NBad == Results[CMax].NBad
							&& Results[CMin].FQS_Bad < Results[CMax].FQS_Bad))))
						{
							CMax = c;
							Value = 0.5 * (Results[CMin].ParValue + Results[c].ParValue);
						}
						else
						{
							CMin = c;
							Value = 0.5 * (Results[c].ParValue + Results[CMax].ParValue);
						}
						break;
					}
				}
				if (Results[NResults - 1].ParValue != Results[NResults - 2].ParValue)
				{
					n = (CMin == NResults - 2 ? CMax : CMin);
					Buffer = (internalFitRoutine ? NameB.left(NameB.indexOf('.')) + QString::number(n) + ".pot"
												 : Processes[N]->getPotFileName(n));
					*LogStream << "Reading potential " << Buffer << '\n';
					LogStream->flush();
					QString FName = Pot->getFileName();
					Pot->readData(Buffer);
					if (internalFitRoutine) 
						Pot->VaryCoefficients(false, NFixedCoefficients, FixedCoefficients);
					else Pot->setFileName(FName);
					Value = Results[NResults - 1].ParValue;
				}
				else
				{
					FitResult *BestResult = 0;
					if (Results[CMin].better(Results + CMax))
					{
						for (n=0; n < NResults - 2; n++)
							if (Results[n].ParValue < Results[CMin].ParValue && (BestResult != 0 ? 
								(fabs(BestResult->ParValue - Results[n].ParValue) < Prec ? 
									Results[n].better(BestResult) 
										: Results[n].ParValue > BestResult->ParValue) : true)) 
								BestResult = Results + n;
						CMax = CMin;
						CMin = BestResult->Index;
						if (Value > BestResult->ParValue && Value < Results[CMax].ParValue)
						{
							if (Value - Results[CMin].ParValue 
									< fabs(Value - Results[NResults - 1].ParValue)
								|| Results[CMax].ParValue - Value 
									< fabs(Value - Results[NResults - 1].ParValue))
							{
								n = (Value - Results[CMin].ParValue < Results[CMax].ParValue - Value ? CMin : CMax);
								Buffer = (internalFitRoutine ? NameB.left(NameB.indexOf('.')) + QString::number(n) + ".pot"
												 : Processes[N]->getPotFileName(n));
								*LogStream << "Reading potential " << Buffer << '\n';
								LogStream->flush();
								QString FName = Pot->getFileName();
								if (internalFitRoutine) 
									Pot->VaryCoefficients(false, NFixedCoefficients, FixedCoefficients);
								else Pot->readData(Buffer);
								Pot->setFileName(FName);
							}
							break;
						}
						c = CMax;
						CMax = NResults - 1;
					}
					else
					{
						for (n=0; n < NResults - 2; n++)
							if (Results[n].ParValue > Results[CMax].ParValue && (BestResult != 0 ? 
								(fabs(BestResult->ParValue - Results[n].ParValue) < Prec ? 
									Results[n].better(BestResult)
									: Results[n].ParValue < BestResult->ParValue) : true)) 
								BestResult = Results + n;
						CMin = CMax;
						CMax = BestResult->Index;
						if (Value < BestResult->ParValue && Value > Results[CMin].ParValue)
						{
							if (Value - Results[CMin].ParValue 
									< fabs(Value - Results[NResults - 1].ParValue)
								|| Results[CMax].ParValue - Value 
									< fabs(Value - Results[NResults - 1].ParValue))
							{
								n = (Value - Results[CMin].ParValue < Results[CMax].ParValue - Value ? CMin : CMax);
								Buffer = (internalFitRoutine ? NameB.left(NameB.indexOf('.')) + QString::number(n) + ".pot"
												 : Processes[N]->getPotFileName(n));
								*LogStream << "Reading potential " << Buffer << '\n';
								LogStream->flush();
								QString FName = Pot->getFileName();
								if (internalFitRoutine) 
									Pot->VaryCoefficients(false, NFixedCoefficients, FixedCoefficients);
								else Pot->readData(Buffer);
								Pot->setFileName(FName);
							}
							break;
						}
						c = CMin;
						CMin = NResults - 1;
					}
					if (Results[CMin].NBadPAL == Results[c].NBadPAL 
						&& Results[c].NBadPAL == Results[CMax].NBadPAL
						&& Results[CMin].NBad == Results[c].NBad
						&& Results[c].NBad == Results[CMax].NBad)
					{
						if (Results[c].FQS_Bad < Results[CMin].FQS_Bad || Results[c].FQS_Bad 
								< Results[CMax].FQS_Bad)
							ParabInterpol(Results[CMin].ParValue, Results[CMin].FQS_Bad,
											Results[c].ParValue, Results[c].FQS_Bad,
											Results[CMax].ParValue, Results[CMax].FQS_Bad,
											Value, V);
							
						else ParabInterpol(Results[CMin].ParValue, Results[CMin].Sigma,
											Results[c].ParValue, Results[c].Sigma,
											Results[CMax].ParValue, Results[CMax].Sigma, Value, V);
					}
					else 
					{
						if (Results[CMin].NBadPAL < Results[CMax].NBadPAL
							|| (Results[CMin].NBadPAL == Results[CMax].NBadPAL
							&& (Results[CMin].NBad < Results[CMax].NBad
							|| (Results[CMin].NBad == Results[CMax].NBad
							&& Results[CMin].FQS_Bad < Results[CMax].FQS_Bad))))
						{
							CMax = c;
							Value = 0.5 * (Results[CMin].ParValue + Results[c].ParValue);
						}
						else
						{
							CMin = c;
							Value = 0.5 * (Results[c].ParValue + Results[CMax].ParValue);
						}
					}
				}
				break;
			case 5:
				Buffer = "Calculations succesfully finished. The estimated paramter value is " 
						+ QString::number(Value, 'g', 1 - int(log(Prec))) + '.';
				*LogStream << Buffer << '\n';
				LogStream->flush();
				QMessageBox::information(this, "MolSpektAnalysis", Buffer);
				return;
				break;
		}
	}
	else if (LP < 0)
	{
		if (NFitsRunning >= NParFits) return;
		FitResult *BResult;
		double V1, V2;
		if (Pos1[N] >= 0 ? (Field[Pos1[N]][Pos2[N]][0] >= 0 ? 
					Results[NResults - 1].better(Results + Field[Pos1[N]][Pos2[N]][0]) : true) : false)
		{
			Field[Pos1[N]][Pos2[N]][0] = NResults - 1;
			if (fitResult == 0 && StartPotIt == StartPotToInclude.constEnd() 
				&& (Improvement || Results[NResults - 1].Sigma <= 1) && Results[NResults - 1].LRelParDev < 10.0)
			{
				if (Pos1[N] > 0 ? (Field[Pos1[N] - 1][Pos2[N]][0] >= 0 ? 
					Improvement && Results[NResults - 1].better(Results + Field[Pos1[N] - 1][Pos2[N]][0]) 
						: Field[Pos1[N] - 1][Pos2[N]][0] != -2) : false)
				{
					Values[N][0] -= Step;
					Field[--Pos1[N]][Pos2[N]][1] = NResults - 1;
					if (Field[Pos1[N]][Pos2[N]][0] == -1) Field[Pos1[N]][Pos2[N]][0]--;
				}
				else if (Pos2[N] > 0 ? (Field[Pos1[N]][Pos2[N] - 1][0] >= 0 ?
					Improvement && Results[NResults - 1].better(Results + Field[Pos1[N]][Pos2[N] - 1][0])
						: Field[Pos1[N]][Pos2[N] - 1][0] != -2) : false)
				{
					Values[N][1] -= Step2;
					Field[Pos1[N]][--Pos2[N]][2] = NResults - 1;
					if (Field[Pos1[N]][Pos2[N]][0] == -1) Field[Pos1[N]][Pos2[N]][0]--;
				}
				else if (Pos1[N] < NPos1 - 1 ? (Field[Pos1[N] + 1][Pos2[N]][0] >= 0 ?
					Improvement && Results[NResults - 1].better(Results + Field[Pos1[N] + 1][Pos2[N]][0])
						: Field[Pos1[N] + 1][Pos2[N]][0] != -2) : false)
				{
					Values[N][0] += Step;
					Field[++Pos1[N]][Pos2[N]][3] = NResults - 1;
					if (Field[Pos1[N]][Pos2[N]][0] == -1) Field[Pos1[N]][Pos2[N]][0]--;
				}
				else if (Pos2[N] < NPos2 - 1 ? (Field[Pos1[N]][Pos2[N] + 1][0] >= 0 ?
					Improvement && Results[NResults - 1].better(Results + Field[Pos1[N]][Pos2[N] + 1][0])
						: Field[Pos1[N]][Pos2[N] + 1][0] != -2) : false)
				{
					Values[N][1] += Step2;
					Field[Pos1[N]][++Pos2[N]][4] = NResults - 1;
					if (Field[Pos1[N]][Pos2[N]][0] == -1) Field[Pos1[N]][Pos2[N]][0]--;
				}
				else Pos1[N] = -1;
			}
			else Pos1[N] = -1;
		}
		else Pos1[N] = -1;
		if (StartPotIt != StartPotToInclude.constEnd())
		{
			printf("LoadStartPot[%d], FilName=%s\n", StartPotIt - StartPotToInclude.constBegin(), StartPotIt->FileName.toLatin1().data());
			LoadPotential(Pot, StartPotIt->FileName, N);
			Pos1[N] = StartPotIt->pos1;
			Pos2[N] = StartPotIt->pos2;
			StartPotIt++;
			Values[N][0] = RStart + Pos1[N] * Step;
			Values[N][1] = RStart2 + Pos2[N] * Step2;
			Field[Pos1[N]][Pos2[N]][0] = -2;
			PreFit = false;
		}
		if (Pos1[N] == -1 && fitResult == 0)
		{
			for (V1 = RStart, n=0, BResult = 0; n < NPos1; n++, V1 += Step)
				for (V2 = RStart2, m=0; m < NPos2; m++, V2 += Step2) 
					if (Field[n][m][0] == -1)
			{
				if (n>0 ? Field[n-1][m][0] >= 0 : false)
					if ((Improvement || Results[Field[n-1][m][0]].Sigma <= 1) && Results[Field[n-1][m][0]].LRelParDev < 10.0
						&& (BResult != 0 ? Results[Field[n-1][m][0]].better(BResult) : true))
				{
					BResult = Results + Field[n-1][m][0];
					Pos1[N] = n;
					Pos2[N] = m;
					Values[N][0] = V1;
					Values[N][1] = V2;
					c=1;
				}
				if (m>0 ? Field[n][m-1][0] >= 0 : false)
					if ((Improvement || Results[Field[n][m-1][0]].Sigma <= 1) && Results[Field[n][m-1][0]].LRelParDev < 10.0
						&& (BResult != 0 ? Results[Field[n][m-1][0]].better(BResult) : true))
				{
					BResult = Results + Field[n][m-1][0];
					Pos1[N] = n;
					Pos2[N] = m;
					Values[N][0] = V1;
					Values[N][1] = V2;
					c=2;
				}
				if (n < NPos1 - 1 ? Field[n+1][m][0] >= 0 : false)
					if ((Improvement || Results[Field[n+1][m][0]].Sigma <= 1) && Results[Field[n+1][m][0]].LRelParDev < 10.0
						&& (BResult != 0 ? Results[Field[n+1][m][0]].better(BResult) : true))
				{
					BResult = Results + Field[n+1][m][0];
					Pos1[N] = n;
					Pos2[N] = m;
					Values[N][0] = V1;
					Values[N][1] = V2;
					c=3;
				}
				if (m < NPos2 - 1 ? Field[n][m+1][0] >= 0 : false)
					if ((Improvement || Results[Field[n][m+1][0]].Sigma <= 1) && Results[Field[n][m+1][0]].LRelParDev < 10.0
						&& (BResult != 0 ? Results[Field[n][m+1][0]].better(BResult) : true))
				{
					BResult = Results + Field[n][m+1][0];
					Pos1[N] = n;
					Pos2[N] = m;
					Values[N][0] = V1;
					Values[N][1] = V2;
					c=4;
				}
			}
			if (BResult != 0)
			{
				LoadPotential(Pot, BResult->Index, N);
				Field[Pos1[N]][Pos2[N]][c] = BResult->Index;
				Field[Pos1[N]][Pos2[N]][0] = -2;
			}
			else if (Improvement)
			{
				for (Values[N][n=0] = RStart; n < NPos1; n++, Values[N][0] += Step)
				{
					for (Values[N][1] = RStart2, m=0; m < NPos2; m++, Values[N][1] += Step2) 
						if (Field[n][m][0] >= 0 && Results[Field[n][m][0]].NBad > 0
							&& ((n>0 && Field[n-1][m][0] >= 0 ? Results[Field[n-1][m][0]].NBad <= 0 : false)
							|| (m>0 && Field[n][m-1][0] >= 0 ? Results[Field[n][m-1][0]].NBad <= 0 : false)
							|| (n < NPos1 - 1 && Field[n+1][m][0] >= 0 ? Results[Field[n+1][m][0]].NBad <= 0 : false)
							|| (m < NPos2 - 1 && Field[n][m+1][0] >= 0 ? Results[Field[n][m+1][0]].NBad <= 0 : false)))
					{
						if (n>0 ? Field[n][m][1] != Field[n-1][m][0] 
							&& Field[n-1][m][0] >= 0 : false)
						{
							LoadPotential(Pot, Field[n][m][1] = Field[n-1][m][0], N);
							break;
						}
						if (m>0 ? Field[n][m][2] != Field[n][m-1][0]
							&& Field[n][m-1][0] >= 0 : false)
						{
							LoadPotential(Pot, Field[n][m][2] = Field[n][m-1][0], N);
							break;
						}
						if (n < NPos1 - 1 ? Field[n][m][3] != Field[n+1][m][0]
							&& Field[n+1][m][0] >= 0 : false)
						{
							LoadPotential(Pot, Field[n][m][3] = Field[n+1][m][0], N);
							break;
						}
						if (m < NPos2 - 1 ? Field[n][m][4] != Field[n][m+1][0] 
							&& Field[n][m+1][0] >= 0 : false)
						{
							LoadPotential(Pot, Field[n][m][4] = Field[n][m+1][0], N);
							break;
						}
						for (c=0; (c < NIpot ? fabs(Values[N][0] - IPotValue1[c]) > Step 
							|| fabs(Values[N][1] - IPotValue2[c]) > Step2 || Tried[n][m][c] : false); c++) ;
						if (c < NIpot)
						{
							if (!InitPots[c].isEmpty()) LoadPotential(Pot, InitPots[c], N);
							else if (StartPot != 0)
							{
								delete Pot;
								Pot = new Potential(*Pot);
								Buffer = NameB.right(NameB.length() - NameB.lastIndexOf(QRegExp("[\\/]")) - 1);
								Pot->VaryCoefficients(false, NFixedCoefficients, FixedCoefficients);
								Pot->setFileName(NameB.left(NameB.lastIndexOf(".")) + QString("w%0.pot").arg(N));
								Pot->setName(Buffer.left(Buffer.lastIndexOf(".")) + "w" + QString::number(N));
							}
							else continue;
							Tried[n][m][c] = true;
							break;
						}
					}
					if (m < NPos2) break;
				}
				if (m == NPos2)
				{
					for (Values[N][n=0] = RStart; n < NPos1; n++, Values[N][0] += Step)
					{
						for (Values[N][1] = RStart2, m=0; m < NPos2; m++, Values[N][1] += Step2) 
							if ((n>0 && n < NPos1 - 1 && Field[n-1][m][0] >= 0 && Field[n+1][m][0] >= 0 ? 
								(Results[Field[n-1][m][0]].better(Results + Field[n][m][0])
								&& Results[Field[n+1][m][0]].betterEqual(Results + Field[n][m][0]))
								|| (Results[Field[n-1][m][0]].betterEqual(Results + Field[n][m][0])
								&& Results[Field[n+1][m][0]].better(Results + Field[n][m][0]))
										: false)
							|| (m>0 && m < NPos2 - 1 && Field[n][m-1][0] >= 0 && Field[n][m+1][0] >= 0 ? 
								(Results[Field[n][m-1][0]].better(Results + Field[n][m][0])
								&& Results[Field[n][m+1][0]].betterEqual(Results + Field[n][m][0]))
								|| (Results[Field[n][m-1][0]].betterEqual(Results + Field[n][m][0])
								&& Results[Field[n][m+1][0]].better(Results + Field[n][m][0]))
										: false))
						{
							if (n>0 ? Field[n][m][1] != Field[n-1][m][0] 
								&& Field[n-1][m][0] >= 0 : false)
							{
								LoadPotential(Pot, Field[n][m][1] = Field[n-1][m][0], N);
								break;
							}
							if (m>0 ? Field[n][m][2] != Field[n][m-1][0]
								&& Field[n][m-1][0] >= 0 : false)
							{
								LoadPotential(Pot, Field[n][m][2] = Field[n][m-1][0], N);
								break;
							}
							if (n < NPos1 - 1 ? Field[n][m][3] != Field[n+1][m][0]
								&& Field[n+1][m][0] >= 0 : false)
							{
								LoadPotential(Pot, Field[n][m][3] = Field[n+1][m][0], N);
								break;
							}
							if (m < NPos2 - 1 ? Field[n][m][4] != Field[n][m+1][0] 
								&& Field[n][m+1][0] >= 0 : false)
							{
								LoadPotential(Pot, Field[n][m][4] = Field[n][m+1][0], N);
								break;
							}
							for (c=0; (c < NIpot ? fabs(Values[N][0] - IPotValue1[c]) > Step 
								|| fabs(Values[N][1] - IPotValue2[c]) > Step2 || Tried[n][m][c] : false); c++) ;
							if (c < NIpot)
							{
								if (!InitPots[c].isEmpty()) LoadPotential(Pot, InitPots[c], N);
								else if (StartPot != 0)
								{
									delete Pot;
									Pot = new Potential(*StartPot);
									Buffer = NameB.right(NameB.length() - NameB.lastIndexOf(QRegExp("[\\/]")) - 1);
									Pot->VaryCoefficients(false, NFixedCoefficients, FixedCoefficients);
									Pot->setFileName(NameB.left(NameB.lastIndexOf(".")) + QString("w%0.pot").arg(N));
									Pot->setName(Buffer.left(Buffer.lastIndexOf(".")) + "w" + QString::number(N));
								}
								else continue;
								Tried[n][m][c] = true;
								break;
							}
						}
						if (m < NPos2) break;
					}
				}
				if (m == NPos2)
				{
					for (Values[N][n=0] = RStart; n < NPos1; n++, Values[N][0] += Step)
					{
						for (Values[N][1] = RStart2, m=0; m < NPos2; m++, Values[N][1] += Step2) 
							if (Field[n][m][0] >= 0 &&
								((n>0 && Field[n-1][m][0] >= 0 ? Results[Field[n-1][m][0]].better(Results + Field[n][m][0]) : false)
								|| (m>0 && Field[n][m-1][0] >= 0 ? Results[Field[n][m-1][0]].better(Results + Field[n][m][0]) : false)
								|| (n < NPos1 - 1 && Field[n+1][m][0] >= 0 ? Results[Field[n+1][m][0]].better(Results + Field[n][m][0])
										: false)
								|| (m < NPos2 - 1 && Field[n][m+1][0] >= 0 ? Results[Field[n][m+1][0]].better(Results + Field[n][m][0])
										: false)))
						{
							if (n>0 ? Field[n][m][1] != Field[n-1][m][0] 
								&& Field[n-1][m][0] >= 0 : false)
							{
								LoadPotential(Pot, Field[n][m][1] = Field[n-1][m][0], N);
								break;
							}
							if (m>0 ? Field[n][m][2] != Field[n][m-1][0]
								&& Field[n][m-1][0] >= 0 : false)
							{
								LoadPotential(Pot, Field[n][m][2] = Field[n][m-1][0], N);
								break;
							}
							if (n < NPos1 - 1 ? Field[n][m][3] != Field[n+1][m][0]
								&& Field[n+1][m][0] >= 0 : false)
							{
								LoadPotential(Pot, Field[n][m][3] = Field[n+1][m][0], N);
								break;
							}
							if (m < NPos2 - 1 ? Field[n][m][4] != Field[n][m+1][0] 
								&& Field[n][m+1][0] >= 0 : false)
							{
								LoadPotential(Pot, Field[n][m][4] = Field[n][m+1][0], N);
								break;
							}
							for (c=0; (c < NIpot ? fabs(Values[N][0] - IPotValue1[c]) > Step 
								|| fabs(Values[N][1] - IPotValue2[c]) > Step2 || Tried[n][m][c] : false); c++) ;
							if (c < NIpot)
							{
								if (!InitPots[c].isEmpty()) LoadPotential(Pot, InitPots[c], N);
								else if (StartPot != 0)
								{
									delete Pot;
									Pot = new Potential(*StartPot);
									Buffer = NameB.right(NameB.length() - NameB.lastIndexOf(QRegExp("[\\/]")) - 1);
									Pot->VaryCoefficients(false, NFixedCoefficients, FixedCoefficients);
									Pot->setFileName(NameB.left(NameB.lastIndexOf(".")) + QString("w%0.pot").arg(N));
									Pot->setName(Buffer.left(Buffer.lastIndexOf(".")) + "w" + QString::number(N));
								}
								else continue;
								Tried[n][m][c] = true;
								break;
							}
						}
						if (m < NPos2) break;
					}
				}
				Pos1[N] = n;
				Pos2[N] = m;
			}
			else Pos2[N] = NPos2;
			if (Pos2[N] == NPos2)
			{
				for (Values[N][n=0] = RStart; n < NPos1; n++, Values[N][0] += Step)
				{
					for (Values[N][1] = RStart2, m=0; m < NPos2; m++, Values[N][1] += Step2) if (Field[n][m][0] == -1)
					{
						for (c=0; (c < NIpot ? fabs(Values[N][0] - IPotValue1[c]) > Step 
								|| fabs(Values[N][1] - IPotValue2[c]) > Step2 || Tried[n][m][c] : false); c++) ;
						if (c < NIpot)
						{
							if (!InitPots[c].isEmpty()) LoadPotential(Pot, InitPots[c], N);
							else 
							{
								Buffer = NameB.right(NameB.length() - NameB.lastIndexOf(QRegExp("[\\/]")) - 1);
								Pot->VaryCoefficients(false, NFixedCoefficients, FixedCoefficients);
								Pot->setFileName(NameB.left(NameB.lastIndexOf(".")) + QString("w%0.pot").arg(N));
								Pot->setName(Buffer.left(Buffer.lastIndexOf(".")) + "w" + QString::number(N));
							}
							break;
						}
					}
					if (m < NPos2) break;
				}
				Pos1[N] = n;
				Pos2[N] = m;
				if (m < NPos2) Field[n][m][0] = -2;
			}
		}
		if (fitResult != 0) Pos1[N] = Pos2[N] = -1;
		fieldWindow->update(Results, NResults, Pos1, Pos2, MaxParFits);
		if (fitResult != 0) return;
		if (Pos2[N] == NPos2)
		{
			Pos1[N] = -1;
			if (Field[0][0][0] == -1) return;
			*LogStream << "Finished Process " << QString::number(N) << '\n';
			LogStream->flush();
			for (n=0; n < NParFits; n++) 
				if (internalFitRoutine ? Pots[n]->isFitRunning() && !Pots[n]->isFitStopped() : Processes[n]->isRunning()) return;
			*LogStream << "Calculations finished\n";
			LogStream->flush();
			QMessageBox::information(this, "MolSpektAnalysis", "Calculations finished");
			return;
		}
	}
	if (NResults > 0 && Method == 0)
	{
		double B;
		for (n=1, c=0, V = fabs(Value - Results[0].ParValue); n < NResults; n++)
			if ((B = fabs(Value - Results[n].ParValue)) < V)
		{
			V=B;
			c=n;
		}
		if (c != NResults - 1) 
		{
			Buffer = (internalFitRoutine ? NameB.left(NameB.indexOf('.')) + QString::number(Results[c].Index) + ".pot"
										 : Processes[N]->getPotFileName(Results[c].Index));
			*LogStream << "Reading potential " << Buffer << '\n';
			LogStream->flush();
			QString FName = Pot->getFileName();
			Pot->readData(Buffer);
			if (internalFitRoutine) Pot->VaryCoefficients(false, NFixedCoefficients, FixedCoefficients);
			else Pot->setFileName(FName);
		}
	}
	if (internalFitRoutine && PreFit ? Pot->getPotType() == SplinePotential : false)
	{
		minR = Pot->getMinR(Pot->getUinf(), 0.1);
		E = Pot->getPoints(minR, maxR, NPoints);
	}
	Pot->getLRCoeffForWriting(NLRC, pLRC, LRC);
	Pot->getAdCorrForWriting(NAdCorr, AdCorr, TAdCorr, PAdCorr, RIso1AdCorr, RIso2AdCorr, AdCorr_Rm, AdCorr_b);
	if (Status == -1)
	{
		Status = 0;
		if (Parameter < NLRC) Value = LRC[Parameter];
		else if (Parameter < NLRC + NAdCorr) Value = AdCorr[Parameter - NLRC];
		else 
		{
			double EIA, EIalpha, EIgamma;
			Pot->getExchangeInt(EIA, EIalpha, EIgamma);
			switch (Parameter - NLRC - NAdCorr)
			{
				case 0:
					Value = EIA;
					break;
				case 1:
					Value = EIalpha;
					break;
				case 2:
					Value = EIgamma;
					break;
				case 3:
					Value = 0.5 * EIalpha * a0_Angstrom;
					break;
			}
		}
	}
	if (Method != 2)
	{
		if (Parameter < NLRC) LRC[Parameter] = Value;
		else if (Parameter < NLRC + NAdCorr) AdCorr[Parameter - NLRC] = Value;
		else
		{
			double EIA, EIalpha, EIgamma;
			Pot->getExchangeInt(EIA, EIalpha, EIgamma);
			switch (Parameter - NLRC - NAdCorr)
			{
				case 0:
					EIA = Value;
					break;
				case 1:
					EIalpha = Value;
					break;
				case 2:
					EIgamma = Value;
					break;
				case 3:
					EIalpha = 2.0 * Value / a0_Angstrom;
					EIgamma = 3.5 / Value - 1.0;
					break;
			}
			Pot->setExchangeInt(EIA, EIalpha, EIgamma);
		}
	}
	/*else if (internalFitRoutine && Pots[N]->getPotType() == SplinePotential)
	{
		WantedCoeffValue *WW = new WantedCoeffValue[2];
		WW[0].CoeffNum = Parameter;
		WW[0].Value = Values[N][0];
		WW[0].Precision = 0.001 * Step;
		WW[1].CoeffNum = Parameter2;
		WW[1].Value = Values[N][1];
		WW[1].Precision = 0.001 * Step2;
		Pots[N]->setWantedLRValues(2, WW);
	}*/
	else
	{
		if (Parameter < NLRC) LRC[Parameter] = Values[N][0];
		else if (Parameter < NLRC + NAdCorr) AdCorr[Parameter - NLRC] = Values[N][0];
		if (Parameter2 < NLRC) LRC[Parameter2] = Values[N][1];
		else if (Parameter2 < NLRC + NAdCorr) AdCorr[Parameter2 - NLRC] = Values[N][1];
		if (Parameter >= NLRC + NAdCorr || Parameter2 >= NLRC + NAdCorr)
		{
			double EIA, EIalpha, EIgamma;
			Pot->getExchangeInt(EIA, EIalpha, EIgamma);
			for (n=0; n<2; n++) if ((m = (n==0 ? Parameter : Parameter2) - NLRC - NAdCorr) >= 0)
				switch (m)
			{
				case 0:
					EIA = Values[N][n];
					break;
				case 1:
					EIalpha = Values[N][n];
					break;
				case 2:
					EIgamma = Values[N][n];
					break;
				case 3:
					EIalpha = 2.0 * Values[N][n] / a0_Angstrom;
					EIgamma = 3.5 / Values[N][n] - 1.0;
					break;
			}
			Pot->setExchangeInt(EIA, EIalpha, EIgamma);
		}
		Pot->cdConnectLR1C();
	}
	if (Method == 2)
	{
		*LogStream << "Started Process " << QString::number(N) << ' ' << QString::number(Values[N][0])
				   << ' ' << QString::number(Values[N][1]) << '\n';
		LogStream->flush();
	}
	if (internalFitRoutine)
	{
		if (PreFit && Pots[N]->getPotType() == SplinePotential) Pots[N]->FitSplinePotToFunction(E, NPoints, minR, maxR);
		//	Pots[N]->FitSplinePot(false, 100);
		//}
		Pots[N]->FitAnaPot(false, false, false, 1e-10, false, false, true);
	}
	else
	{
		Processes[N]->setPotential(Pot);
		Processes[N]->start(NResults);
	}
	if (++NFitsRunning < NParFits && RecursiveCalls)
	{
		for (N=0; (N < NParFits ? Pos1[N] != -1 : false); N++) ;
		if (N < NParFits)
		{
			NoFit = true;
			ProcessFinished(N, -2);
		}
	}
	if (E != 0) delete[] E;
}

void MonteCarloSim::LoadPotential(Potential* Pot, int ResultNum, int PN)
{
	QString Buffer = (internalFitRoutine ? 
						NameB.left(NameB.indexOf('.')) + QString::number(Results[ResultNum].Index) + ".pot"
						: Processes[Results[ResultNum].ProcessNum]->getPotFileName(Results[ResultNum].Index));
	*LogStream << "Process " + QString::number(PN) + ": Reading potential " 
			   << Buffer << '\n';
	QString FName = Pot->getFileName();
	if (!Pot->readData(Buffer))
	{
		printf(("Error reading File \"" + Buffer + "\"").toLatin1().data());
		*LogStream << "Error reading File \"" << Buffer << "\"\n";
	}
	LogStream->flush();
	if (internalFitRoutine) 
	{
		Pot->VaryCoefficients(false, NFixedCoefficients, FixedCoefficients);
		Pot->setFileName(Buffer.left(Buffer.length() - 4) + "w" + QString::number(PN) + ".pot");
		Pot->setName(Pot->getName() + "w" + QString::number(PN));
	}
	else Pot->setFileName(FName);
}

void MonteCarloSim::LoadPotential(Potential* Pot, QString FileName, int PN)
{
	*LogStream << "Process " + QString::number(PN) + ": Reading potential " 
			   << FileName << '\n';
	QString FName = Pot->getFileName();
	if (!Pot->readData(FileName))
	{
		printf(("Error reading File \"" + FileName + "\"").toLatin1().data());
		*LogStream << "Error reading File \"" << FileName << "\"\n";
	}
	LogStream->flush();
	if (internalFitRoutine) 
	{
		Pot->VaryCoefficients(false, NFixedCoefficients, FixedCoefficients);
		Pot->setFileName(FileName.left(FileName.length() - 4) + "w" + QString::number(PN) + ".pot");
		Pot->setName(Pot->getName() + "w" + QString::number(PN));
	}
	else Pot->setFileName(FName);
}

void MonteCarloSim::connected()
{
	/*if (!isClient) return;
	int n, N = MessagesToSend.count();
	QTextStream S(Socket);
	for (n=0; n<N; n++) S << MessagesToSend[n] << '\n';*/
}

void MonteCarloSim::disconnected()
{

}

void MonteCarloSim::newConnection()
{

}

void MonteCarloSim::readyRead()
{

}

void MonteCarloSim::timerEvent(QTimerEvent* /*Event*/)
{
   
}
