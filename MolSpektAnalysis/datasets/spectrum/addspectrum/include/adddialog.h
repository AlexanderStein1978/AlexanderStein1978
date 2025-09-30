//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//

#ifndef ADDDIALOG_H
#define ADDDIALOG_H

/**
	@author Alexander Stein <AlexanderStein@t-online.de>
*/

#include "DiagWindow.h"
#include "Spektrum.h"

#include <QMenu>
#include <QToolButton>

class QComboBox;
class QLineEdit;
class QLabel;
class QPushButton;
class QRadioButton;
class QCheckBox;
class QImage;
class QTextEdit;

class MainWindow;
class DiagWindow;
class FitData;
class TermTable;
class AddDialog;
class AddSpectrum;
class ResidualFit;
struct TermEnergy;

enum autoUpdateResult {autoUpdateSuccess, autoUpdatePartialFailure, 
					   autoUpdateCompleteFailure, autoUpdateInsufficientData};
enum MarkedPeakLevelState {hiddenSaved = -1, hiddenNew, savedAssignment, newAssignment,
						   markedSaved, savedToDelete, newToDelete};


class AddDialog : public QWidget
{
	Q_OBJECT
	
	public:
    	AddDialog(MainWindow *MW, Spektrum *Spekt);
    	~AddDialog();
		void enableOptimize(bool enable);
		
	private slots:
		void Apply();
		void MolChanged(QString Name);
		void MethodChanged();
		void UpdateSpektList();

	private:
		QComboBox *MolBox, *USBox, *LSBox, *IsoBox, *RIsoBox, *CompBox, *SpektBox, *RSpektBox;
		QComboBox *RCompBox, *TUnit;
		QLineEdit *L1E, *L2E, *TempE, *RiE, *RaE, *ResE, *RJsE;
		QLabel *L1, *L2, *L14;
		QPushButton *ApplyB, *CloseB;
		QRadioButton *AddByBand, *AddByProg, *AddByProgFCF;
		QCheckBox *AutoAnalyse, *DrawImage, *UpdateAllVs;
		MainWindow *MW;
		AddSpectrum *Image;
		Spektrum *Spekt;
		TermTable *TT, *LT;
		double ****TData, **IData;
		int NJ, NC, JStart, JStep, NumIso;
		bool Marked;
		QTextEdit *LogWindow;
};

#endif
