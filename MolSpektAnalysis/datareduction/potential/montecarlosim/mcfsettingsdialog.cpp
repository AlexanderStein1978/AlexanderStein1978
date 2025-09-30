//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//


#include "mcfsettingsdialog.h"
#include "MainWindow.h"
#include "potential.h"

#include <QPushButton>
#include <QLineEdit>
#include <QLabel>
#include <QRadioButton>
#include <QComboBox>
#include <QCheckBox>
#include <QIntValidator>
#include <QGridLayout>
#include <QButtonGroup>
#include <QMessageBox>
#include <QDir>
#include <QFileDialog>


MCFSettingsDialog::MCFSettingsDialog(MainWindow* parent): QDialog(parent)
{
	int n, N = parent->getNumPotentials();
	MW = parent;
	setWindowTitle("Settings for fitting procedure");
	QGridLayout *L = new QGridLayout(this), *L1 = new QGridLayout();
	QButtonGroup *intExt = new QButtonGroup(this);
	intExt->addButton(internalFitRoutine = new QRadioButton("internal fit routine", this));
	L1->addWidget(internalFitRoutine, 0, 0);
	intExt->addButton(externalProgram = new QRadioButton("external fit program", this));
	L1->addWidget(externalProgram, 0, 1);
	internalFitRoutine->setChecked(true);
	L->addLayout(L1, 0, 0, 1, 3);
	L->addWidget(new QLabel("Dir of init potentials:", this), 0, 3);
	L->addWidget(InitPotDirE = new QLineEdit(this), 0, 4);
	L->addWidget(sInitPotDirB = new QPushButton("...", this), 0, 5);
	L->addWidget(new QLabel("Result logfile:", this), 1, 0);
	L->addWidget(ResultLogFileE = new QLineEdit(this), 1, 1);
	L->addWidget(sResultLogFileB = new QPushButton("...", this), 1, 2);
	L->addWidget(PotentialLabel = new QLabel("Potential:", this), 1, 3);
	L->addWidget(PotentialNameE = new QLineEdit(this), 1, 4);
	L->addWidget(sPotentialNameB = new QPushButton("...", this), 1, 5);
	L->addWidget(PotBox = new QComboBox(this), 1, 4);
	for (n=0; n<N; n++) PotBox->addItem(parent->getPotential(n)->getName());
	PotentialLabel->setVisible(false);
	sPotentialNameB->setEnabled(false);
	L->addWidget(FitProgLabel = new QLabel("Fit program:", this), 2, 0);
	L->addWidget(FitProgE = new QLineEdit(this), 2, 1);
	L->addWidget(sFitProgB = new QPushButton("...", this), 2, 2);
	FitProgLabel->setEnabled(false);
	FitProgE->setEnabled(false);
	sFitProgB->setEnabled(false);
	L->addWidget(withImprovement = new QCheckBox("with improvment", this), 2, 3);
	withImprovement->setChecked(false);
	L->addWidget(new QLabel("Num parallel fits:", this), 2, 4);
	L->addWidget(numParIt = new QLineEdit("1", this), 2, 5);
	numParIt->setValidator(NumParItValid = new QIntValidator(1, 1, numParIt));
	L->addWidget(DataFileLabel = new QLabel("Data file:", this), 3, 0);
	L->addWidget(FitDataE = new QLineEdit(this), 3, 1);
	L->addWidget(sFitDataB = new QPushButton("...", this), 3, 2);
	DataFileLabel->setEnabled(false);
	FitDataE->setEnabled(false);
	sFitDataB->setEnabled(false);
	L->addWidget(isClient = new QCheckBox("Client", this), 3, 3);
	L->addWidget(new QLabel("Server IP:", this), 3, 4);
	L->addWidget(ServerIP = new QLineEdit("", this), 3, 5);
	ServerIP->setEnabled(false);
	isClient->setCheckState(Qt::Unchecked);
	L->addWidget(inputFile1Label = new QLabel("Input file 1:", this), 4, 0);
	L->addWidget(FitInFile1E = new QLineEdit(this), 4, 1);
	L->addWidget(sFitInFile1B = new QPushButton("...", this), 4, 2);
	inputFile1Label->setEnabled(false);
	FitInFile1E->setEnabled(false);
	sFitInFile1B->setEnabled(false);
	L->addWidget(new QLabel("Procedure to apply:", this), 4, 3, 1, 3);
	L->addWidget(inputFile2Label = new QLabel("Input file 2:", this), 5, 0);
	L->addWidget(FitInFile2E = new QLineEdit(this), 5, 1);
	L->addWidget(sFitInFile2B = new QPushButton("...", this), 5, 2);
	inputFile2Label->setEnabled(false);
	FitInFile2E->setEnabled(false);
	sFitInFile2B->setEnabled(false);
	QButtonGroup *method = new QButtonGroup(this);
	method->addButton(searchMin = new QRadioButton("Search min:", this));
	L->addWidget(searchMin, 5, 3);
	searchMin->setChecked(true);
	method->addButton(estimateUncert = new QRadioButton("Esitmate uncert:", this));
	L->addWidget(estimateUncert, 5, 4);
	method->addButton(search2DPlane = new QRadioButton("Scan 2D plane:", this));
	L->addWidget(search2DPlane, 5, 5);
	L->addWidget(inffoutLabel = new QLabel("In.f.f. output:", this), 6, 0);
	L->addWidget(FitOutInFileE = new QLineEdit(this), 6, 1);
	L->addWidget(sFitOutInFileB = new QPushButton("...", this), 6, 2);
	inffoutLabel->setEnabled(false);
	FitOutInFileE->setEnabled(false);
	sFitOutInFileB->setEnabled(false);
	L->addWidget(new QLabel("Parameter:", this), 6, 3);
	L->addWidget(ParameterB = new QComboBox(this), 6, 4);
	L->addWidget(ParameterB2 = new QComboBox(this), 6, 5);
	ParameterB2->setEnabled(false);
	L->addWidget(inFileLabel = new QLabel("In file name:", this), 7, 0);
	L->addWidget(FitInFileNameE = new QLineEdit(this), 7, 1);
	L->addWidget(sFitInFileNameB = new QPushButton("...", this), 7, 2);
	inFileLabel->setEnabled(false);
	FitInFileNameE->setEnabled(false);
	sFitInFileNameB->setEnabled(false);
	L->addWidget(new QLabel("S.int. begin:", this), 7, 3);
	L->addWidget(RStartE = new QLineEdit(this), 7, 4);
	L->addWidget(RStartE2 = new QLineEdit(this), 7, 5);
	RStartE2->setEnabled(false);
	L->addWidget(FitOutFileLabel = new QLabel("Out file name:", this), 8, 0);
	L->addWidget(FitOutFileNameE = new QLineEdit(this), 8, 1);
	L->addWidget(sFitOutFileNameB = new QPushButton("...", this), 8, 2);
	FitOutFileLabel->setEnabled(false);
	FitOutFileNameE->setEnabled(false);
	sFitOutFileNameB->setEnabled(false);
	L->addWidget(new QLabel("S.int. end:", this), 8, 3);
	L->addWidget(RStopE = new QLineEdit(this), 8, 4);
	L->addWidget(RStopE2 = new QLineEdit(this), 8, 5);
	L->addWidget(BadListNameLabel = new QLabel("File n. badlist:", this), 9, 0);
	L->addWidget(BadListNameE = new QLineEdit(this), 9, 1);
	L->addWidget(sBadListNameB = new QPushButton("...", this), 9, 2);
	BadListNameLabel->setEnabled(false);
	BadListNameE->setEnabled(false);
	sBadListNameB->setEnabled(false);
	L->addWidget(new QLabel("Init step s.:", this), 9, 3);
	L->addWidget(StepE = new QLineEdit(this), 9, 4);
	L->addWidget(StepE2 = new QLineEdit(this), 9, 5);
	L->addWidget(PotTempFileNameLabel = new QLabel("F.n. temp pot:", this), 10, 0);
	L->addWidget(PotTempFileNameE = new QLineEdit(this), 10, 1);
	L->addWidget(sPotTempFileNameB = new QPushButton("...", this), 10, 2);
	PotTempFileNameLabel->setEnabled(false);
	PotTempFileNameE->setEnabled(false);
	sPotTempFileNameB->setEnabled(false);
	L->addWidget(new QLabel("Rel. precision:", this), 10, 3);
	L->addWidget(PrecE = new QLineEdit("0.01", this), 10, 4);
	L->addWidget(new QLabel("Work directory:", this), 11, 0);
	L->addWidget(WorkDirE = new QLineEdit(this), 11, 1);
	L->addWidget(sWorkDirB = new QPushButton("...", this), 11, 2);
	QGridLayout *L2 = new QGridLayout;
	L->addLayout(L2, 11, 3, 1, 3);
	L2->addWidget(WorkDirE2 = new QLineEdit(this), 0, 0);
	L2->addWidget(sWorkDirB2 = new QPushButton("...", this), 0, 1);
	L2->addWidget(WorkDirE3 = new QLineEdit(this), 0, 2);
	L2->addWidget(sWorkDirB3 = new QPushButton("...", this), 0, 3);
	WorkDirE2->setEnabled(false);
	WorkDirE3->setEnabled(false);
	sWorkDirB2->setEnabled(false);
	sWorkDirB3->setEnabled(false);
	L->setRowMinimumHeight(12, 20);
	L->addWidget(OK = new QPushButton("OK", this), 13, 0);
	L->addWidget(Cancel = new QPushButton("Cancel", this), 13, 5);
	ParameterB2->setEditable(false);
	ParameterB->setEditable(false);
	connect(isClient, SIGNAL(toggled(bool)), this, SLOT(statusChanged(bool)));
	connect(sResultLogFileB, SIGNAL(clicked()), this, SLOT(sResultLogFile()));
	connect(sFitProgB, SIGNAL(clicked()), this, SLOT(sFitProgD()));
	connect(sFitDataB, SIGNAL(clicked()), this, SLOT(sFitDataD()));
	connect(sFitInFile1B, SIGNAL(clicked()), this, SLOT(sFitInFile1D()));
	connect(sFitInFile2B, SIGNAL(clicked()), this, SLOT(sFitInFile2D()));
	connect(sFitOutInFileB, SIGNAL(clicked()), this, SLOT(sFitOutInFileD()));
	connect(sFitInFileNameB, SIGNAL(clicked()), this, SLOT(sFitInFileNameD()));
	connect(sFitOutFileNameB, SIGNAL(clicked()), this, SLOT(sFitOutFileNameD()));
	connect(sBadListNameB, SIGNAL(clicked()), this, SLOT(sBadListNameD()));
	connect(sPotTempFileNameB, SIGNAL(clicked()), this, SLOT(sPotTempFileNameD()));
	connect(sWorkDirB, SIGNAL(clicked()), this, SLOT(sWorkDir()));
	connect(sWorkDirB2, SIGNAL(clicked()), this, SLOT(sWorkDir2()));
	connect(sWorkDirB3, SIGNAL(clicked()), this, SLOT(sWorkDir3()));
	connect(InitPotDirE, SIGNAL(editingFinished()), this, SLOT(potentialChanged()));
	connect(PotBox, SIGNAL(currentIndexChanged(int)), this, SLOT(potentialChanged()));
	connect(sPotentialNameB, SIGNAL(clicked()), this, SLOT(sPotentialName()));
	connect(sInitPotDirB, SIGNAL(clicked()), this, SLOT(sInitPotDir()));
	connect(search2DPlane, SIGNAL(toggled(bool)), this, SLOT(MethodChanged(bool)));
	connect(externalProgram, SIGNAL(toggled(bool)), this, SLOT(intExtFitProgChanged(bool)));
	connect(numParIt, SIGNAL(editingFinished()), this, SLOT(numParItChanged()));
	connect(OK, SIGNAL(clicked()), this, SLOT(accept()));
	connect(Cancel, SIGNAL(clicked()), this, SLOT(reject()));
	potentialChanged();
}

void MCFSettingsDialog::accept()
{
	QFile F;
	if (ParameterB->count() == 0)
	{
		QMessageBox::information(this, "MolSpektAnalysis", 
			"You have to specify a valid directory with initial potentials!");
		return;
	}
	QString B;
	if (ResultLogFileE->text().isEmpty())
	{
		QMessageBox::information(this, "MolSpektAnalysis", "You have to specify a name for the logfile of the results!");
		return;
	}
	if (InitPotDirE->text().isEmpty())
	{
		QMessageBox::information(this, "MolSpektAnalysis", "You have to specify a directory to search for start potentials!");
		return;
	}
	if (externalProgram->isChecked())
	{
		F.setFileName(FitProgE->text());
		if (!F.exists())
		{
			QMessageBox::information(this, "MolSpektAnalysis", "The selected program for fitting does not exist!");
			return;
		}
		F.setFileName(FitDataE->text());
		if (!F.exists())
		{
			QMessageBox::information(this, "MolSpektAnalysis", "The selected data file does not exist!");
			return;
		}
		F.setFileName(FitInFile1E->text());
		if (!F.exists())
		{
			QMessageBox::information(this, "MolSpektAnalysis", "The selected first input file does not exist!");
			return;
		}
		if (!(B = FitInFile2E->text()).isEmpty())
		{
			F.setFileName(B);
			if (!F.exists())
			{
				QMessageBox::information(this, "MolSpektAnalysis", "The selected second input file does not exist!");
				return;
			}
		}
		F.setFileName(FitOutInFileE->text());
		if (!F.exists())
		{
			QMessageBox::information(this, "MolSpektAnalysis", 
				"The selected input file for obtaining results in the case of a program crash (In.f.f. output) does not exist!");
			return;
		}
		if (FitInFileNameE->text().isEmpty())
		{
			QMessageBox::information(this, "MolSpektAnalysis", "You have to specify the file name for the fit input files!");
			return;
		}
		if (FitOutFileNameE->text().isEmpty())
		{
			QMessageBox::information(this, "MolSpektAnalysis", "You have to give the ouput file name of the fitting program!");
			return;
		}
		if (BadListNameE->text().isEmpty())
		{
			QMessageBox::information(this, "MolSpektAnalysis", "You have to give the badlist file name of the fitting program!");
			return;
		}
		if (PotTempFileNameE->text().isEmpty())
		{
			QMessageBox::information(this, "MolSpektAnalysis", 
									 "You have to give the file name for temporary potentials of the fitting program!");
			return;
		}
		if (PotentialNameE->text().isEmpty())
		{
			QMessageBox::information(this, "MolSpektAnalysis",
									"You have to give the file name the potential has in the selected input files!");
			return;
		}
	}
	if (PrecE->text().toDouble() <= 0.0 && !isClient->isChecked())
	{
		QMessageBox::information(this, "MolSpektAnalysis", "You have to specify a valid precision >0 for paramter to opimize!");
		return;
	}
	if (RStartE->text().toDouble() - RStopE->text().toDouble() >= 0.0 && searchMin->isChecked())
	{
		QMessageBox::information(this, "MolSpektAnalysis", "You have to specify a valid interval for searching!");
		return;
	}
	if (StepE->text().toDouble() == 0.0 && !isClient->isChecked())
	{
		QMessageBox::information(this, "MolSpektAnalysis", "You have to specify a valid initial step size != 0!");
		return;
	}
	if (internalFitRoutine->isChecked() && PotBox->currentIndex() == -1)
	{
		QMessageBox::information(this, "MolSpektAnalysis", 
			"If you want to use the internal fitting routine you have to have loaded and selected an initial potential.");
		return;
	}
	int m, NLRC, *pLRC, NAdCorr, TAdCorr, PAdCorr, PT1 = ParameterB->currentIndex(), PT2 = ParameterB2->currentIndex(), *FixC, NFixC = 0;
	double *LRC, *adCorr, RIso1AdCorr, RIso2AdCorr, AdCorr_Rm, AdCorr_b;
	double Par1Start = RStartE->text().toDouble(), Par1Stop = RStopE->text().toDouble(), Par1;
	double Par2Start = RStartE2->text().toDouble(), Par2Stop = RStopE2->text().toDouble(), Par2;
	double Step1 = StepE->text().toDouble(), Step2 = StepE2->text().toDouble();
	bool external = externalProgram->isChecked(), nothingInside = true, M2 = search2DPlane->isChecked();
	QDir D(InitPotDirE->text());
	QStringList L = D.entryList(QDir::Files), L2;
	Potential *Pot = StartPot;
	if (Pot != 0) Pot->getFixedCoefficients(NFixC, FixC);
	for (m = (external ? 0 : -1); m < L.count() && nothingInside; m++)
	{
		if (m==0) Pot = (external ? new Potential() : new Potential(*Pot));
		Par1 = Par2 = 0.0;
		if (m>=0) if (!Pot->readData(InitPotDirE->text() + DIRSEP + L[m])) continue;
		Pot->getLRCoeffForReading(NLRC, pLRC, LRC);
		if (PT1 < NLRC) Par1 = LRC[PT1];
		if (PT2 < NLRC) Par2 = LRC[PT2];
		if (PT1 >= NLRC || PT2 >= NLRC)
		{
			Pot->getAdCorrForReading(NAdCorr, adCorr, TAdCorr, PAdCorr, RIso1AdCorr, RIso2AdCorr, AdCorr_Rm, AdCorr_b);
			if (PT1 >= NLRC && PT1 < NLRC + NAdCorr) Par1 = adCorr[PT1 - NLRC];
			if (PT2 >= NLRC && PT2 < NLRC + NAdCorr) Par2 = adCorr[PT2 - NLRC];
			if ((PT1 >= NLRC + NAdCorr || PT2 >= NLRC + NAdCorr) && Pot->isExcIntAvailable())
			{
				double A, alpha, gamma;
				Pot->getExchangeInt(A, alpha, gamma);
				if (PT1 >= NLRC + NAdCorr) switch(PT1 - NLRC - NAdCorr)
				{
					case 0:
						Par1 = A;
						break;
					case 1:
						Par1 = alpha;
						break;
					case 2:
						Par1 = gamma;
						break;
					case 3:
						Par1 = 0.5 * alpha * a0_Angstrom;
						break;
					default:
						Par1 = 0.0;
						break;
				}
				if (PT2 >= NLRC + NAdCorr) switch(PT2 - NLRC - NAdCorr)
				{
					case 0:
						Par2 = A;
						break;
					case 1:
						Par2 = alpha;
						break;
					case 2:
						Par2 = gamma;
						break;
					case 3:
						Par2 = 0.5 * alpha * a0_Angstrom;
						break;
					default:
						Par2 = 0.0;
						break;
				}
			}
		}
		if (Par1 >= Par1Start - Step1 && Par1 <= Par1Stop + Step1 && (!M2 || (Par2 >= Par2Start - Step2 && Par2 <= Par2Stop + Step2)))
				nothingInside = false;
		delete[] pLRC;
		delete[] LRC;
		if ((PT1 > NLRC || PT2 > NLRC) && NAdCorr > 0) delete[] adCorr;
	}
	if (m>0) 
	{
		StartPot = Pot;
		Pot->VaryCoefficients(false, NFixC, FixC);
	}
	else delete Pot;
	if (nothingInside)
	{
		QMessageBox::information(this, "MolSpektAnalysis", 
			"Error: No start potential found which lies inside the selected parameter range!");
		return;
	}
    QDialog::accept();
}

void MCFSettingsDialog::getResults(QString &WorkDir, QString &WorkDir2, QString &WorkDir3, 
								   QString &ResultLogFile, QString& FitProg, QString& FitData, 
								   QString& FitInFile1, QString& FitInFile2, 
								   QString& FitOutInFile, QString& FitInFileName, 
								   QString& FitOutFileName, QString& BadListName, 
								   QString &InitPotDir, QString &PotentialName, 
								   QString& PotTempFileName, int& Method, double& Prec, 
								   double& RStart, double &RStart2, double& RStop, 
								   double &RStop2, double& Step, double &Step2, int& Parameter, int &Parameter2, QString &serverIP, 
								   bool &internalFit, bool &Improvement, int &NumParIt, Potential *&Pot)
{
	WorkDir = WorkDirE->text();
	WorkDir2 = WorkDirE2->text();
	WorkDir3 = WorkDirE3->text();
	ResultLogFile = ResultLogFileE->text();
	FitProg = FitProgE->text();
	FitData = FitDataE->text();
	FitInFile1 = FitInFile1E->text();
	FitInFile2 = FitInFile2E->text();
	FitOutInFile = FitOutInFileE->text();
	FitInFileName = FitInFileNameE->text();
	FitOutFileName = FitOutFileNameE->text();
	BadListName = BadListNameE->text();
	InitPotDir = InitPotDirE->text();
	PotentialName = PotentialNameE->text();
	PotTempFileName = PotTempFileNameE->text();
	Method = (estimateUncert->isChecked() ? 0 : (searchMin->isChecked() ? 1 : 2));
	Prec = PrecE->text().toDouble();
	RStart = RStartE->text().toDouble();
	RStart2 = RStartE2->text().toDouble();
	RStop = RStopE->text().toDouble();
	RStop2 = RStopE2->text().toDouble();
	Step = StepE->text().toDouble();
	Step2 = StepE2->text().toDouble();
	Parameter = ParameterB->currentIndex();
	Parameter2 = ParameterB2->currentIndex();
	internalFit = internalFitRoutine->isChecked();
	Improvement = withImprovement->isChecked();
	NumParIt = numParIt->text().toInt();
	Pot = StartPot;
	if (isClient->isChecked()) serverIP = ServerIP->text();
}

void MCFSettingsDialog::intExtFitProgChanged(bool external)
{
	if (!external && PotBox->count() == 0)
	{
		QMessageBox::information(this, "MolSpektAnalysis", 
			"If you want to use the internal fitting routine you have to have loaded an initial potential.");
		externalProgram->setChecked(true);
		return;
	}
	int n;
	FitProgLabel->setEnabled(external);
	FitProgE->setEnabled(external);
	sFitProgB->setEnabled(external);
	DataFileLabel->setEnabled(external);
	FitDataE->setEnabled(external);
	sFitDataB->setEnabled(external);
	inputFile1Label->setEnabled(external);
	FitInFile1E->setEnabled(external);
	sFitInFile1B->setEnabled(external);
	inputFile2Label->setEnabled(external);
	FitInFile2E->setEnabled(external);
	sFitInFile2B->setEnabled(external);
	inffoutLabel->setEnabled(external);
	FitOutInFileE->setEnabled(external);
	sFitOutInFileB->setEnabled(external);
	inFileLabel->setEnabled(external);
	FitInFileNameE->setEnabled(external);
	sFitInFileNameB->setEnabled(external);
	FitOutFileLabel->setEnabled(external);
	FitOutFileNameE->setEnabled(external);
	sFitOutFileNameB->setEnabled(external);
	BadListNameLabel->setEnabled(external);
	BadListNameE->setEnabled(external);
	sBadListNameB->setEnabled(external);
	PotTempFileNameLabel->setEnabled(external);
	PotTempFileNameE->setEnabled(external);
	sPotTempFileNameB->setEnabled(external);
	PotBox->setVisible(!external);
	PotentialNameE->setVisible(external);
	sPotentialNameB->setEnabled(external);
	PotentialLabel->setText(external ? "Potential filename:" : "Potential:");
	if (search2DPlane->isChecked())
	{
		n = numParIt->text().toInt();
		if (n>=2)
		{
			WorkDirE2->setEnabled(external);
			sWorkDirB2->setEnabled(external);
			if (n==3)
			{
				sWorkDirB3->setEnabled(external);
				WorkDirE3->setEnabled(external);
			}
		}
		if (external && n>3) numParIt->setText("3");
		NumParItValid->setTop(external ? 3 : 100);
	}
	potentialChanged();
}

void MCFSettingsDialog::MethodChanged(bool s2DChecked)
{
	if (externalProgram->isChecked())
	{
		int n = numParIt->text().toInt();
		if (n>=2)
		{
			WorkDirE2->setEnabled(s2DChecked);
			sWorkDirB2->setEnabled(s2DChecked);
			if (n==3)
			{
				WorkDirE3->setEnabled(s2DChecked);
				sWorkDirB3->setEnabled(s2DChecked);
			}
		}
	}
	PrecE->setEnabled(!s2DChecked);
	RStartE2->setEnabled(s2DChecked);
	RStopE2->setEnabled(s2DChecked);
	StepE2->setEnabled(s2DChecked);
	ParameterB2->setEnabled(s2DChecked);
	NumParItValid->setTop(s2DChecked ? (externalProgram->isChecked() ? 3 : 100) : 1);
}

void MCFSettingsDialog::numParItChanged()
{
	if (internalFitRoutine->isChecked()) return;
	int n = numParIt->text().toInt();
	WorkDirE2->setEnabled(n>=2);
	sWorkDirB2->setEnabled(n>=2);
	WorkDirE3->setEnabled(n==3);
	sWorkDirB3->setEnabled(n==3);
}

void MCFSettingsDialog::potentialChanged()
{
	int n, m, NLRC, *pLRC, NAdCorr, TAdCorr, PAdCorr;
	double *LRC, *adCorr, RIso1AdCorr, RIso2AdCorr, AdCorr_Rm, AdCorr_b;
	bool external = externalProgram->isChecked();
	QDir D(InitPotDirE->text());
	if (!D.exists() && external) return;
	QStringList L = D.entryList(QDir::Files), L2;
	StartPot = (external ? 0 : MW->getPotential(PotBox->currentIndex()));
	Potential *Pot = (external ? new Potential() : new Potential(*StartPot));
	QString Buffer;
	ParameterB->clear();
	ParameterB2->clear();
	for (m = (external ? 0 : -1); m < L.count() && ParameterB->count() == 0; m++)
	{
		if (m>=0) if (!Pot->readData(InitPotDirE->text() + DIRSEP + L[m])) continue;
		Pot->getLRCoeffForReading(NLRC, pLRC, LRC);
		Pot->getAdCorrForReading(NAdCorr, adCorr, TAdCorr, PAdCorr, RIso1AdCorr, RIso2AdCorr, AdCorr_Rm, AdCorr_b);
		for (n=0; n < NLRC; n++) 
		{
			ParameterB->addItem(Buffer = "C" + QString::number(pLRC[n]));
			ParameterB2->addItem(Buffer);
		}
		for (n=0; n < NAdCorr; n++) 
		{
			ParameterB->addItem(Buffer = "AdCorr" + QString::number(n));
			ParameterB2->addItem(Buffer);
		}
		if (Pot->isExcIntAvailable())
		{
			ParameterB->addItems(L2 << "ExcInt A" << "ExcInt alpha" << "ExcInt gamma" << "ExcInt beta");
			ParameterB2->addItems(L2);
		}
		if (NLRC > 0)
		{
			delete[] pLRC;
			delete[] LRC;
		}
		if (NAdCorr > 0) delete[] adCorr;
	}
	delete Pot;
}

void MCFSettingsDialog::sBadListNameD()
{
	QString CurText = BadListNameE->text();
	QString Name = QFileDialog::getOpenFileName(this, "Please select file name of badlist", 
												(!CurText.isEmpty() ? CurText : (!Dir.isEmpty() ? Dir : MW->getDir())));
	if (!Name.isEmpty())
	{
		int n = Name.lastIndexOf('/');
		if (n > 0) MW->setDir(Dir = Name.left(n+1));
		BadListNameE->setText(Name.right(Name.length() - n - 1));
	}
}

void MCFSettingsDialog::sFitDataD()
{
	QString CurText = FitDataE->text();
	QString Name = QFileDialog::getOpenFileName(this, "Please select file with fit data", 
												(!CurText.isEmpty() ? CurText : (!Dir.isEmpty() ? Dir : MW->getDir())));
	if (!Name.isEmpty())
	{
		int n = Name.lastIndexOf('/');
		if (n > 0) MW->setDir(Dir = Name.left(n+1));
		FitDataE->setText(Name);
	}
}

void MCFSettingsDialog::sFitInFile1D()
{
	QString CurText = FitInFile1E->text();
	QString Name = QFileDialog::getOpenFileName(this, "Please select fit input file 1", 
												(!CurText.isEmpty() ? CurText : (!Dir.isEmpty() ? Dir : MW->getDir())));
	if (!Name.isEmpty())
	{
		int n = Name.lastIndexOf('/');
		if (n > 0) MW->setDir(Dir = Name.left(n+1));
		FitInFile1E->setText(Name);
	}
}

void MCFSettingsDialog::sFitInFile2D()
{
	QString CurText = FitInFile2E->text();
	QString Name = QFileDialog::getOpenFileName(this, "Please select fit input file 2", 
												(!CurText.isEmpty() ? CurText : (!Dir.isEmpty() ? Dir : MW->getDir())));
	if (!Name.isEmpty())
	{
		int n = Name.lastIndexOf('/');
		if (n > 0) MW->setDir(Dir = Name.left(n+1));
		FitInFile2E->setText(Name);
	}
}

void MCFSettingsDialog::sFitInFileNameD()
{
	QString CurText = FitInFileNameE->text();
	QString Name = QFileDialog::getOpenFileName(this, "Please select file with name of fit input file", 
												(!CurText.isEmpty() ? CurText : (!Dir.isEmpty() ? Dir : MW->getDir())));
	if (!Name.isEmpty())
	{
		int n = Name.lastIndexOf('/');
		if (n > 0) MW->setDir(Dir = Name.left(n+1));
		FitInFileNameE->setText(Name.right(Name.length() - n - 1));
	}
}

void MCFSettingsDialog::sFitOutFileNameD()
{
	QString CurText = FitOutFileNameE->text();
	QString Name = QFileDialog::getOpenFileName(this, "Please select file with name of fit output file", 
												(!CurText.isEmpty() ? CurText : (!Dir.isEmpty() ? Dir : MW->getDir())));
	if (!Name.isEmpty())
	{
		int n = Name.lastIndexOf('/');
		if (n > 0) MW->setDir(Dir = Name.left(n+1));
		FitOutFileNameE->setText(Name.right(Name.length() - n - 1));
	}
}

void MCFSettingsDialog::sFitOutInFileD()
{
	QString CurText = FitOutInFileE->text();
	QString Name = QFileDialog::getOpenFileName(this, "Please select fit input file for getting data output", 
												(!CurText.isEmpty() ? CurText : (!Dir.isEmpty() ? Dir : MW->getDir())));
	if (!Name.isEmpty())
	{
		int n = Name.lastIndexOf('/');
		if (n > 0) MW->setDir(Dir = Name.left(n+1));
		FitOutInFileE->setText(Name);
	}
}

void MCFSettingsDialog::sFitProgD()
{
	QString CurText = FitProgE->text();
	QString Name = QFileDialog::getOpenFileName(this, "Please select fit program", 
												(!CurText.isEmpty() ? CurText : (!Dir.isEmpty() ? Dir : MW->getDir())), 
												"Programs (*.exe)");
	if (!Name.isEmpty())
	{
		int n = Name.lastIndexOf('/');
		if (n > 0) MW->setDir(Dir = Name.left(n+1));
		FitProgE->setText(Name);
	}
}

void MCFSettingsDialog::sInitPotDir()
{
	QString CurText = InitPotDirE->text();
	QString Name = QFileDialog::getExistingDirectory(this, 
									"Please select directory with initial potentials",
									(CurText.isEmpty() ? (Dir.isEmpty() ? MW->getDir() : Dir) : CurText));
	if (!Name.isEmpty())
	{
		int n = Name.lastIndexOf('/');
		if (n>0) MW->setDir(Dir = Name.left(n+1));
		InitPotDirE->setText(Name);
		potentialChanged();
	}
}

void MCFSettingsDialog::sPotentialName()
{
	QString CurText = PotentialNameE->text();
	QString Name = QFileDialog::getOpenFileName(this, 
				"Please select file with name of potential given in the input files.", 
				(!CurText.isEmpty() ? CurText : (!Dir.isEmpty() ? Dir : MW->getDir())));
	if (!Name.isEmpty())
	{
		int n = Name.lastIndexOf('/');
		if (n > 0) MW->setDir(Dir = Name.left(n+1));
		PotentialNameE->setText(Name.right(Name.length() - n - 1));
	}
}

void MCFSettingsDialog::sPotTempFileNameD()
{
	QString CurText = PotTempFileNameE->text();
	QString Name = QFileDialog::getOpenFileName(this, "Please select file with name of temporary potentials", 
												(!CurText.isEmpty() ? CurText : (!Dir.isEmpty() ? Dir : MW->getDir())));
	if (!Name.isEmpty())
	{
		int n = Name.lastIndexOf('/');
		if (n > 0) MW->setDir(Dir = Name.left(n+1));
		PotTempFileNameE->setText(Name.right(Name.length() - n - 1));
	}
}

void MCFSettingsDialog::sResultLogFile()
{
	QString CurText = ResultLogFileE->text();
	QString Name = QFileDialog::getSaveFileName(this, "Please select file name for logfile", 
												(!CurText.isEmpty() ? CurText : (!Dir.isEmpty() ? Dir : MW->getDir())));
	if (!Name.isEmpty())
	{
		int n = Name.lastIndexOf('/');
		if (n > 0) MW->setDir(Dir = Name.left(n+1));
		ResultLogFileE->setText(Name);
	}
}

void MCFSettingsDialog::statusChanged(bool newStatus)
{
	if (newStatus)
	{
		ServerIP->setEnabled(true);
		WorkDirE2->setEnabled(true);
		WorkDirE3->setEnabled(true);
		sWorkDirB2->setEnabled(true);
		sWorkDirB3->setEnabled(true);
		search2DPlane->setEnabled(false);
		searchMin->setEnabled(false);
		estimateUncert->setEnabled(false);
		ParameterB->setEnabled(false);
		ParameterB2->setEnabled(false);
		PrecE->setEnabled(false);
		RStartE->setEnabled(false);
		RStartE2->setEnabled(false);
		RStopE->setEnabled(false);
		RStopE2->setEnabled(false);
		StepE->setEnabled(false);
		StepE2->setEnabled(false);
	}
	else
	{
		bool s2DChecked = search2DPlane->isChecked();
		search2DPlane->setEnabled(true);
		searchMin->setEnabled(true);
		estimateUncert->setEnabled(true);
		WorkDirE2->setEnabled(s2DChecked);
		WorkDirE3->setEnabled(s2DChecked);
		sWorkDirB2->setEnabled(s2DChecked);
		sWorkDirB3->setEnabled(s2DChecked);
		PrecE->setEnabled(!s2DChecked);
		RStartE->setEnabled(true);
		RStartE2->setEnabled(s2DChecked);
		RStopE->setEnabled(true);
		RStopE2->setEnabled(s2DChecked);
		StepE->setEnabled(true);
		StepE2->setEnabled(s2DChecked);
		ParameterB->setEnabled(true);
		ParameterB2->setEnabled(s2DChecked);
	}
}

void MCFSettingsDialog::sWorkDir()
{
	QString CurText = WorkDirE->text();
	QString Name = QFileDialog::getExistingDirectory(this, "Please select work directory", 
												(!CurText.isEmpty() ? CurText : (!Dir.isEmpty() ? Dir : MW->getDir())));
	if (!Name.isEmpty())
	{
		int n = Name.lastIndexOf('/');
		if (n > 0) MW->setDir(Dir = Name.left(n+1));
		WorkDirE->setText(Name);
	}
}

void MCFSettingsDialog::sWorkDir2()
{
	QString CurText = WorkDirE2->text();
	QString Name = QFileDialog::getExistingDirectory(this, "Please select work directory 2", 
												(!CurText.isEmpty() ? CurText : (!Dir.isEmpty() ? Dir : MW->getDir())));
	if (!Name.isEmpty())
	{
		int n = Name.lastIndexOf('/');
		if (n > 0) MW->setDir(Dir = Name.left(n+1));
		WorkDirE2->setText(Name);
	}
}

void MCFSettingsDialog::sWorkDir3()
{
	QString CurText = WorkDirE3->text();
	QString Name = QFileDialog::getExistingDirectory(this, "Please select work directory 3", 
												(!CurText.isEmpty() ? CurText : (!Dir.isEmpty() ? Dir : MW->getDir())));
	if (!Name.isEmpty())
	{
		int n = Name.lastIndexOf('/');
		if (n > 0) MW->setDir(Dir = Name.left(n+1));
		WorkDirE3->setText(Name);
	}
}
