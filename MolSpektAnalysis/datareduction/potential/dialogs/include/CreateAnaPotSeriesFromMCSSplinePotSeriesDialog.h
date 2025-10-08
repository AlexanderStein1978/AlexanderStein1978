//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//


#ifndef CREATEANAPOTSERIESFROMMCSSPLINEPOTSERIESDIALOG_H
#define CREATEANAPOTSERIESFROMMCSSPLINEPOTSERIESDIALOG_H


#include <QDialog>

class MainWindow;
class ElState;

class QCheckBox;
class QLineEdit;
class QComboBox;


class CreateAnaPotSeriesFromMCSSplinePotSeriesDialog : public QDialog
{
	Q_OBJECT
public:
    CreateAnaPotSeriesFromMCSSplinePotSeriesDialog(MainWindow* parent);
	void getData(QString &MolFN, QString &StateN, QString &FDDir, QString &SPDir, QString &APDir, 
				 bool &improveAnaPots, bool &UseSvd, bool &UseLeveMarq, int &MaxIt, double &Prec);

private slots:
	void showFDDirDialog();
	void showSplinePotDirDialog();
	void showAnaPotDirDialog();
	void molChanged(int i);
	void improveChanged(bool ckecked);
	void MaxItChanged(const QString &text);
	
private:
	MainWindow *MW;
	ElState *State;
	QLineEdit *FDDirE, *SplinePotDirE, *AnaPotDirE, *PrecE;
	QComboBox *MoleculeBox, *StateBox, *MaxIt;
	QCheckBox *improve, *useSVD, *useLM;
	QPushButton *SPDB;
};

#endif
