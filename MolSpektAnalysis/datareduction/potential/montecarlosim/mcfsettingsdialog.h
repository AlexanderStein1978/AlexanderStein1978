//
// C++ Interface: MCFSettingsDialog
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2011 - 2019
//
// Copyright: See README file that comes with this source code
//
//


#ifndef MCFSETTINGSDIALOG_H
#define MCFSETTINGSDIALOG_H


#include <QDialog>

class MainWindow;
class Potential;


class MCFSettingsDialog : public QDialog
{
	Q_OBJECT
	
public:
	MCFSettingsDialog(MainWindow *MW);
	void getResults(QString &WorkDir, QString &WorkDir2, QString &WorkDir3, QString &ResultLogFile, QString &FitProg, QString &FitData, 
					QString &FitInFile1, QString &FitInFile2, QString &FitOutInFile, QString &FitInFileName, QString &FitOutFileName, 
				    QString &BadListName, QString &InitPotDir, QString &PotentialName, QString &PotTempFileName, int &Method, 
				    double &Prec, double &RStart, double &RStart2, double &RStop, double &RStop2, double &Step, double &Step2, 
				    int &Parameter, int &Parameter2, QString &serverIP, bool &internalFit, bool &Improvement, int &numParIt, 
				    Potential *&Pot);
	
private slots:
	void statusChanged(bool newStatus);
	void sWorkDir();
	void sWorkDir2();
	void sWorkDir3();
	void sResultLogFile();
	void sInitPotDir();
	void sFitProgD();
	void sFitDataD();
	void sFitInFile1D();
	void sFitInFile2D();
	void sFitOutInFileD();
	void sFitInFileNameD();
	void sFitOutFileNameD();
	void sBadListNameD();
	void sPotentialName();
	void sPotTempFileNameD();
	void accept();
	void potentialChanged();
	void MethodChanged(bool s2DChecked);
	void intExtFitProgChanged(bool internal);
	void numParItChanged();
	
private:
	QLineEdit *FitProgE, *FitDataE, *FitInFile1E, *FitInFile2E, *FitOutInFileE, *WorkDirE2, *WorkDirE3;
	QLineEdit *FitInFileNameE, *FitOutFileNameE, *BadListNameE, *PotTempFileNameE;
	QLineEdit *PrecE, *RStartE, *RStopE, *StepE, *ResultLogFileE, *WorkDirE, *InitPotDirE;
	QLineEdit *RStartE2, *RStopE2, *StepE2, *PotentialNameE, *ServerIP, *numParIt;
	QRadioButton *searchMin, *estimateUncert, *search2DPlane, *internalFitRoutine, *externalProgram;
	QPushButton *sFitProgB, *sFitDataB, *sFitInFile1B, *sFitInFile2B, *sFitOutInFileB;
	QPushButton *sFitInFileNameB, *sFitOutFileNameB, *sBadListNameB, *sPotTempFileNameB;
	QPushButton *OK, *Cancel, *sResultLogFileB, *sWorkDirB, *sWorkDirB2, *sWorkDirB3;
	QPushButton *sInitPotDirB, *sPotentialNameB;
	QComboBox *ParameterB, *ParameterB2, *PotBox;
	QLabel *WorkDirL, *FitProgLabel, *DataFileLabel, *inputFile1Label, *inputFile2Label, *inffoutLabel, *inFileLabel, *FitOutFileLabel;
	QLabel *BadListNameLabel, *PotTempFileNameLabel, *PotentialLabel;
	QCheckBox *isClient, *withImprovement;
	MainWindow *MW;
	QString Dir;
	QIntValidator *NumParItValid;
	Potential *StartPot;
};

#endif
