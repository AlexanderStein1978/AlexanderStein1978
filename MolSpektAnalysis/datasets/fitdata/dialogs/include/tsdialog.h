//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//


#ifndef TSDIALOG_H
#define TSDIALOG_H


#include <QDialog>
#include <QCheckBox>

class LineTable;
class Molecule;
class ElState;

class QListWidget;
class QComboBox;
class QPushButton;


class TSDialog : public QDialog
{
	Q_OBJECT
	
	public:
		TSDialog(QWidget *parent, Molecule *M, ElState *S);
		~TSDialog();
		
		int getNumSelected();
		LineTable *getSelected(int n);
		bool getKeepOther();
		
		inline bool getUpdateSA()
		{
			return UpdateSA->isChecked();
		}
		
	private slots:
		void addTable();
		void removeTable();
		void OK();
        void transBoxChanged(int StateIndex);
		
		inline void KeepOtherChanged(bool Checked)
		{
			UpdateSA->setEnabled(Checked);
		}
		
	private:
        LineTable **Tables;
        int NTables;
        Molecule *Mol;
		QListWidget *List;
        QComboBox *Box, *TransBox;
		QCheckBox *KeepOther, *UpdateSA;
		QPushButton *Add, *Remove, *OKB, *Cancel;
};

#endif

