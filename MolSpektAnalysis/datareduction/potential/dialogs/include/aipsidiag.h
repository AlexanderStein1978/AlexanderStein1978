//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//

#ifndef AIPSIDIAG_H
#define AIPSIDIAG_H

/**
	@author Alexander Stein <AlexanderStein@t-online.de>
*/

#include <QDialog>

class QCheckBox;
class QRadioButton;
class QLineEdit;
class QListWidget;
class QComboBox;
class QStringList;
class QString;

class MainWindow;

class AIPSIDialog : public QWidget
{
	Q_OBJECT
			
public:
    AIPSIDialog(MainWindow *MW);
    ~AIPSIDialog();
		
private slots:
	void showFileDialog();
	void PotFileBoxChanged(int i);
	void AtMolChecked(int S);
	void ImportData();
	
private:
	void openFile(QString FileName, bool &GF, int &anRow, int &anCol, double *&aR, double **&aE);
	
	MainWindow *MW;
	int *NCol, *NRow, NFiles;
	double **R, ***E;
	QComboBox *LUBox, *EUBox, *PotNumBox, *PotFileBox, *PotentialBox, *MoleculeBox;
	QListWidget *FileList;
	QLineEdit *Source, *Path;
	QRadioButton *PotOffsetMinSero, *PotOffsetASero, *PotOffsetMatchA, *PotOffsetConst;
	QCheckBox *AtMol;
	QStringList Files;
};

#endif
