//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//


#include "tablewidgetwindow.h"
#include "isotab.h"
#include "molecule.h"
#include "Spektrum.h"
#include "utils.h"
#include "heapsort.h"

#include <QFile>
#include <QTextStream>
#include <QMessageBox>
#include <QProgressBar>


TableWidgetWindow::TableWidgetWindow(Type typ, MainWindow *MW, Molecule *M) : TableWindow(typ, MW, M)
{
	int w = width(), h = height();
	int w3 = w / 3, w2 = w / 2, w4 = w / 4;
	QGridLayout *L, *L2;
	QBoxLayout *Layout;
	switch (Typ)
	{
		case TextTable1:
		case TextTable2:
			setFilter("Text tables (*.dat)");
			setFileExt(".dat");
			Tab = new QTableWidget(this);
			Layout = new QBoxLayout(QBoxLayout::LeftToRight, this);
			Layout->addWidget(Tab);
			break;
		case DunhamTable:
			Tab = new QTableWidget(MaxDunCoefficients, 5, this);
			Tab->setGeometry(0, 70, w, h - 70);

			vMLabel = new QLabel("Max v:", this);
			vMLabel->setGeometry(10, 40, 40, 20);
			vMax = new QLineEdit("", this);
			vMax->setGeometry(50, 40, w3 - 60, 20);
			connect(vMax, SIGNAL(textEdited(QString)), this, SLOT(vMaxChanged()));

			JMLabel = new QLabel("Max J:", this);
			JMLabel->setGeometry(w3 + 10, 40, 40, 20);
			JMax = new QLineEdit("", this);
			JMax->setGeometry(w3 + 50, 40, w3 - 60, 20);
			connect(JMax, SIGNAL(textEdited(QString)), this, SLOT(JMaxChanged()));

			errLabel = new QLabel("Error:", this);
			errLabel->setGeometry(2 * w3 + 10, 40, 40, 20);
			error = new QLineEdit("", this);
			error->setGeometry(2 * w3 + 50, 40, w3 - 60, 20);
			connect(error, SIGNAL(textEdited(QString)), this, SLOT(errorChanged()));

			setFilter("Dunhame coefficents (*.dun)");
			setFileExt(".dun");
			break;
		case TermEnergyView:
			Tab = new QTableWidget(this);
			Tab->setGeometry(0, 70, w, h - 70);

			ViewLabel = new QLabel("View:", this);
			ViewLabel->setGeometry(10, 40, 40, 20);
			View = new QComboBox(this);
			View->setGeometry(50, 40, w3 - 60, 20);
			View->setEditable(false);

			IsoLabel = new QLabel("Isotopologue:", this);
			IsoLabel->setGeometry(w3 + 10, 40, 80, 20);
			Iso = new QComboBox(this);
			Iso->setGeometry(w3 + 90, 40, w3 - 100, 20);
			Iso->setEditable(false);

			CompLabel = new QLabel("Component:", this);
			CompLabel->setGeometry(2 * w3 + 10, 70, 40, 20);
			Comp = new QComboBox(this);
			Comp->setGeometry(2 * w3 + 80, 40, w3 - 90, 20);
			Comp->setEditable(false);

			Source->setGeometry(220, 10, w - 380, 20);

			TLabel = new QLabel("T:", this);
			TLabel->setGeometry(w - 140, 10, 20, 20);
			Temp = new QLineEdit("20", this);
			Temp->setGeometry(w - 120, 10, 60, 20);
			Temp->setEnabled(false);
			TUnit = new QComboBox(this);
			TUnit->setGeometry(w - 60, 10, 50, 20);
			TUnit->setEditable(false);
			TUnit->setEnabled(false);

			Name->setReadOnly(true);
			Source->setReadOnly(true);
			break;
		case FranckCondonView:
			Tab = new QTableWidget(this);
			Tab->setGeometry(0, 70, w, h - 70);

			MolLabel = new QLabel("Molecule:", this);
			MolLabel->setGeometry(10, 10, 70, 20);
			MolBox = new QComboBox(this);
			MolBox->setGeometry(80, 10, w4 - 100, 20);
			MolBox->setEditable(false);

			ViewLabel = new QLabel("View:", this);
			ViewLabel->setGeometry(10, 40, 40, 20);
			View = new QComboBox(this);
			View->setGeometry(50, 40, w4 - 60, 20);
			View->setEditable(false);

			IsoLabel = new QLabel("Isotopologue:", this);
			IsoLabel->setGeometry(w4 + 10, 10, 70, 20);
			Iso = new QComboBox(this);
			Iso->setGeometry(w4 + 80, 10, w4 - 100, 20);
			Iso->setEditable(false);

			TLabel = new QLabel("T:", this);
			TLabel->setGeometry(w4 + 10, 40, 20, 20);
			Temp = new QLineEdit("20", this);
			Temp->setGeometry(w4 + 30, 40, 60, 20);
			Temp->setEnabled(false);
			TUnit = new QComboBox(this);
			TUnit->setGeometry(w4 + 90, 40, 50, 20);
			TUnit->setEditable(false);
			TUnit->setEnabled(false);

			USLabel = new QLabel("Upper state:", this);
			USLabel->setGeometry(w2 + 10, 10, 70, 20);
			USBox = new QComboBox(this);
			USBox->setGeometry(w2 + 80, 10, w4 - 100, 20);
			USBox->setEditable(false);

			LSLabel = new QLabel("Lower state:", this);
			LSLabel->setGeometry(3 * w4 + 10, 10, 70, 20);
			LSBox = new QComboBox(this);
			LSBox->setGeometry(3 * w4 + 80, 10, w4 - 100, 20);
			LSBox->setEditable(false);

			JsLabel = new QLabel("J':", this);
			JsLabel->setGeometry(w2 + 10, 40, 30, 20);
			Js = new QLineEdit("0", this);
			Js->setGeometry(w2 + 40, 40, 50, 20);

			JssLabel = new QLabel("J'':", this);
			JssLabel->setGeometry(w2 + 110, 40, 30, 20);
			Jss = new QLineEdit("0", this);
			Jss->setGeometry(w2 + 140, 40, 50, 20);

			Calc = new QPushButton("Calculate", this);
			Calc->setGeometry(3.25 * w4, 40, 0.75 * w4 - 10, 20);
			break;
		case PotData:
			Tab = new QTableWidget(this);
			Tab->setGeometry(0, 40, w, h - 40);
			setFilter("Measured term energies (*.mterm)");
			setFileExt(".mterm");
			break;
		case FranckCondonTable:
			setFilter("FCF tables (*.fcf)");
			setFileExt(".fcf");
			L = new QGridLayout(this);
			L->addWidget(new QLabel("Name:", this), 0, 0);
			L->addWidget(Name = new QLineEdit("", this), 0, 1);
			connect(Name, SIGNAL(editingFinished()), this, SLOT(setName()));
			L->setColumnStretch(1, 1);
			L->addWidget(new QLabel("Upper potential:", this), 0, 2);
			L->addWidget(Pot1 = new QLineEdit("", this), 0, 3);
			Pot1->setReadOnly(true);
			L->setColumnStretch(3, 1);
			L->addWidget(new QLabel("Lower potential:", this), 0, 4);
			L->addWidget(Pot2 = new QLineEdit("", this), 0, 5);
			Pot2->setReadOnly(true);
			L->setColumnStretch(5, 1);
			L->addWidget(new QLabel("Date calculated:", this), 1, 0);
			L->addWidget(Date = new QLineEdit("", this), 1, 1);
			L->addWidget(new QLabel("View isotopologue:", this), 1, 2);
			L->addWidget(Iso = new QComboBox(this), 1, 3);
			Iso->setEditable(false);
			L2 = new QGridLayout;
			L2->addWidget(new QLabel("v':", this), 0, 0);
			L2->addWidget(JsB = new QComboBox(this), 0, 1);
			JsB->setEditable(false);
			L2->addWidget(new QLabel("v'':", this), 0, 2);
			L2->addWidget(JssB = new QComboBox(this), 0, 3);
			JssB->setEditable(false);
			L->addLayout(L2, 1, 4, 1, 2);
			L->addWidget(Tab = new QTableWidget(this), 2, 0, 1, 6);
			break;
		case FitSeriesResultTable:
			L = new QGridLayout(this);
			L->addWidget(new QLabel("Number of parallel fits:", this), 0, 0);
			L->addWidget(NumParFits = new QLineEdit(this), 0, 1);
			L->addWidget(Progress = new QProgressBar(this), 0, 2);
			Progress->setFormat("%v iterations of %m finished");
			Progress->setMinimum(0);
			L->addWidget(Tab = new QTableWidget(this), 1, 0, 1, 3);
			break;
        case External:
            break;
		default:
			printf("TableWindow::TableWindow: Error: The type %d is not a valid type for a tablewindow!", Typ);
			break;
	}
	if (Typ != TermEnergyTable && Typ != FitDataSet)
	{
		connect(Tab, SIGNAL(itemChanged(QTableWidgetItem*)), this, SLOT(Changed()));
		connect(Tab, SIGNAL(itemChanged(QTableWidgetItem*)), this, SLOT(tabItemChanged(QTableWidgetItem*)));
	}
}

TableWidgetWindow::~TableWidgetWindow()
{
}

void TableWidgetWindow::resizeHelper(QRect& G)
{
	int w = width();
	int w3 = w / 3, w2 = w / 2, w4 = w / 4;
	switch(Typ)
	{
		case DunhamTable:
			vMax->setGeometry(50, 40, w3 - 60, 20);
			JMLabel->setGeometry(w3 + 10, 40, 40, 20);
			JMax->setGeometry(w3 + 50, 40, w3 - 60, 20);
			errLabel->setGeometry(2 * w3 + 10, 40, 40, 20);
			error->setGeometry(2 * w3 + 50, 40, w3 - 60, 20);
			G = Tab->geometry();
			G.setWidth(width());
			G.setHeight(height() - 70);
			Tab->setGeometry(G);
			break;
		case TermEnergyView:
			View->setGeometry(50, 40, w3 - 60, 20);
			IsoLabel->setGeometry(w3 + 10, 40, 80, 20);
			Iso->setGeometry(w3 + 90, 40, w3 - 100, 20);
			CompLabel->setGeometry(2 * w3 + 10, 40, 70, 20);
			Comp->setGeometry(2 * w3 + 80, 40, w3 - 90, 20);
			G = Tab->geometry();
			G.setWidth(width());
			G.setHeight(height() - 70);
			Tab->setGeometry(G);
			Source->setGeometry(220, 10, w - 380, 20);
			TLabel->setGeometry(w - 140, 10, 20, 20);
			Temp->setGeometry(w - 120, 10, 60, 20);
			TUnit->setGeometry(w - 60, 10, 50, 20);
			break;
		case FranckCondonView:
			Tab->setGeometry(0, 70, w, height() - 70);
			MolLabel->setGeometry(10, 10, 70, 20);
			MolBox->setGeometry(80, 10, w4 - 100, 20);
			ViewLabel->setGeometry(10, 40, 40, 20);
			View->setGeometry(50, 40, w4 - 60, 20);
			IsoLabel->setGeometry(w4 + 10, 10, 70, 20);
			Iso->setGeometry(w4 + 80, 10, w4 - 100, 20);
			TLabel->setGeometry(w4 + 10, 40, 20, 20);
			Temp->setGeometry(w4 + 30, 40, 60, 20);
			TUnit->setGeometry(w4 + 90, 40, 50, 20);
			USLabel->setGeometry(w2 + 10, 10, 70, 20);
			USBox->setGeometry(w2 + 80, 10, w4 - 100, 20);
			LSLabel->setGeometry(3 * w4 + 10, 10, 70, 20);
			LSBox->setGeometry(3 * w4 + 80, 10, w4 - 100, 20);
			JsLabel->setGeometry(w2 + 10, 40, 30, 20);
			Js->setGeometry(w2 + 40, 40, 50, 20);
			JssLabel->setGeometry(w2 + 110, 40, 30, 20);
			Jss->setGeometry(w2 + 140, 40, 50, 20);
			Calc->setGeometry(3.25 * w4, 40, 0.75 * w4 - 10, 20);
			break;
		case PotData:
			G = Tab->geometry();
			G.setWidth(width());
			G.setHeight(height() - 40);
			Tab->setGeometry(G);
			break;
		default:
			break;
	}
}

bool TableWidgetWindow::readData(QTextStream& S)
{
	int r, cc, lc, n;
	QString Buffer;
	QStringList L;
	bool Success;
	QRegExp specialSectionStart = GetStartSpecialPartRegExp();
	Tab->blockSignals(true);
	for (r=0, cc = Tab->columnCount(); !S.atEnd(); ++r)
	{
		if (Tab->rowCount() == r) Tab->setRowCount(r + 100);
		if ((Buffer = S.readLine()).left(15) == "Column titles: ")
		{
			if (Typ == -1)
			{
				L = Buffer.right(Buffer.length() - 15).split(Spacer);
				Tab->setColumnCount(cc = L.count());
				Tab->setHorizontalHeaderLabels(L);
			}
			Buffer = S.readLine();
		}
		else if (!specialSectionStart.isEmpty() && Buffer.indexOf(specialSectionStart) >= 0 && (!(Success = ReadSpecialPart(S, Buffer)) || S.atEnd())) break;
		L = Buffer.split(Spacer);
		if ((lc = L.count()) > cc) Tab->setColumnCount(cc = lc);
		for (n=0; n < lc; ++n) Tab->setItem(r, n, new QTableWidgetItem(L[n]));
		while (n < cc) Tab->setItem(r, ++n, new QTableWidgetItem(""));
	}
	if (Typ != 4 && lc < cc) --r;
	Tab->setRowCount(r);
	Tab->blockSignals(false);
	return true;
}

void TableWidgetWindow::writeData(QTextStream& S)
{
	int c, r, R = Tab->rowCount(), C = Tab->columnCount();
	QTableWidgetItem *I;
	QStringList L;
	for (c=0; c<C; ++c) L << ((I = Tab->horizontalHeaderItem(c)) != 0 ? I->text() : "?");
	S << "Column titles: " << L.join(Spacer) << "\n";
	for (r=0; r < R; ++r)
	{
		L.clear();
		for (c=0; c<C; ++c) L << ((I = Tab->item(r, c)) != 0 ? I->text() : "");
		S << L.join(Spacer).remove('\n') << "\n";
	}
}

void TableWidgetWindow::setTabDimensions(int NRows, int NCols)
{
	Tab->setRowCount(NRows);
	Tab->setColumnCount(NCols);
}

void TableWidgetWindow::getTabDimensions(int& NRows, int& NCols)
{
	NRows = Tab->rowCount();
	NCols = Tab->columnCount();
}

void TableWidgetWindow::setRowData(int Row, QString* Data)
{
	int n, N = Tab->columnCount();
	for (n=0; n<N; n++) Tab->setItem(Row, n, new QTableWidgetItem(Data[n]));
}

void TableWidgetWindow::AddRow()
{
	if (Tab == 0) return;
	int r, c, nr = Tab->rowCount(), nc = Tab->columnCount(), cr = Tab->currentRow();
    if (cr == -1) cr = 0;
	Tab->blockSignals(true);
	Tab->setRowCount(nr + 1);
	for (r = nr; r > cr; r--) for (c=0; c < nc; c++) Tab->setItem(r, c, Tab->takeItem(r-1, c));
    for (c=0; c < nc; ++c) Tab->setItem(cr, c, new QTableWidgetItem(""));
	Tab->blockSignals(false);
	Changed();
}

void TableWidgetWindow::DeleteRows()
{
    if (Tab == 0) return;
	QList<QTableWidgetSelectionRange> SR = Tab->selectedRanges();
	Tab->blockSignals(true);
	int j, k, n = Tab->rowCount(), c, nc = Tab->columnCount();
	for (j=0; j < SR.size(); j++) for (k=SR[j].topRow(); k <= SR[j].bottomRow(); k++)
			for (c=0; c < nc; c++) Tab->setItem(k, c, 0);
	for (j=0; (j < n ? Tab->item(j, 0) != 0 : false); j++) ;
    for (k=j; k < n; k++) if (Tab->item(k, 0) != 0)
	{
		for (c=0; c < nc; c++) Tab->setItem(j, c, Tab->takeItem(k, c));
		j++;
	}
    Tab->setRowCount(j);
	Tab->blockSignals(false);
	Changed();
}

bool TableWidgetWindow::checkAllConnections(int FileColumn)
{
    bool ret = true, changed = false;
    Spektrum *Spectrum = new Spektrum(MW);
    QString SpektPath, NewPath;
    QStringList CheckedPaths, Replacements;
    int n, i;
    Tab->blockSignals(true);
    for (n=0; n < Tab->rowCount(); ++n) if (Tab->item(n, FileColumn) != 0 && (SpektPath = Tab->item(n, FileColumn)->text()).indexOf('(') == -1)
    {
        if ((i = CheckedPaths.indexOf(SpektPath)) == -1)
        {
            if (Spectrum->readData(SpektPath)) NewPath = Spectrum->getFileName();
            else
            {
                NewPath = SpektPath;
                ret = false;
            }
            CheckedPaths << SpektPath;
            Replacements << NewPath;
        }
        else NewPath = Replacements[i];
        if (NewPath != SpektPath)
        {
            Tab->item(n, FileColumn)->setText(NewPath);
            changed = true;
        }
    }
    Tab->blockSignals(false);
    if (changed) Changed();
    delete Spectrum;
    return ret;
}

void TableWidgetWindow::copyRows(int& nR, int& nC, QString**& Data)
{
	if (Tab == 0) return;
	int i, j, r, c;
	QList<QTableWidgetSelectionRange> SR = Tab->selectedRanges();
	if (Data != 0) Destroy(Data, nR);
	nC = Tab->columnCount();
	for (nR = i = 0; i < SR.size(); i++) nR += SR[i].rowCount();
	Data = CreateQString(nR, nC);
	for (i=j=0; i < SR.size(); i++) for (r = SR[i].topRow(); r <= SR[i].bottomRow(); r++, j++)
		for (c=0; c < nC; c++) Data[j][c] = (Tab->item(r, c) != 0 ? Tab->item(r, c)->text() : "");
}

void TableWidgetWindow::cutRows(int &nR, int &nC, QString **&Data)
{
	if (Tab == 0) return;
	int i, j, r, c, NR = Tab->rowCount();
	QList<QTableWidgetSelectionRange> SR = Tab->selectedRanges();
	if (Data != 0) Destroy(Data, nR);
	nC = Tab->columnCount();
	for (nR = i = 0; i < SR.size(); i++) nR += SR[i].rowCount();
	Data = CreateQString(nR, nC);
	Tab->blockSignals(true);
	for (i=j=0; i < SR.size(); i++) for (r = SR[i].topRow(); r <= SR[i].bottomRow(); r++, j++)
	{
		for (c=0; c < nC; c++)
			Data[j][c] = (Tab->item(r, c) != 0 ? Tab->item(r, c)->text() : "");
		Tab->setItem(r, 0, 0);
	}
	for (r=0; (r < NR ? Tab->item(r, 0) != 0 : false); r++) ;
	for (j=r; j < NR; j++) if (Tab->item(j, 0) != 0)
	{
		for (c=0; c < nC; c++) Tab->setItem(r, c, Tab->takeItem(j, c));
		r++;
	}
	Tab->setRowCount(r);
	Tab->blockSignals(false);
	Changed();
}

void TableWidgetWindow::insertRows(int nR, int nC, QString **Data)
{
	int r = Tab->rowCount(), c, n, NC = Tab->columnCount();
	Tab->blockSignals(true);
	Tab->setRowCount(r + nR);
	for (n=0; n < nR; n++, r++) for (c=0; c < NC; c++)
		Tab->setItem(r, c, new QTableWidgetItem(c < nC ? Data[n][c] : ""));
	Tab->blockSignals(false);
	Changed();
}

void TableWidgetWindow::MarkLines(int* rN, int N)
{
	if (Tab == 0) return;
	int n=0, r1, MC = Tab->columnCount() - 1;
	QList<QTableWidgetSelectionRange> SL = Tab->selectedRanges();
	for (r1 = 0; r1 < SL.count(); r1++) Tab->setRangeSelected(SL[r1], false);
	Tab->scrollToItem(Tab->item(rN[0], 0), QAbstractItemView::PositionAtTop);
	Tab->setCurrentItem(Tab->item(rN[0], 0));
	for (r1 = 1; r1 != 0; ) for (r1 = 0, n=1; n<N; n++) if (rN[n] < rN[n-1])
	{
		r1 = rN[n];
		rN[n] = rN[n-1];
		rN[n-1] = r1;
	}
	for (n=0; n<N; )
	{
		for (r1 = rN[n++]; n<N  && rN[n] == rN[n-1] + 1; n++) ;
		Tab->setRangeSelected(QTableWidgetSelectionRange(r1, 0, rN[n-1], MC), true);
	}
	if (!isVisible()) show();
	activateWindow();
	Tab->setFocus();
}

void TableWidgetWindow::setIsoIcon(int Col, int c)
{
	if (Tab == nullptr || molecule == nullptr) return;
	IsoTab* Iso = molecule->getIso();
	int r, n, N = Tab->rowCount();
	for (r=0; r<N; r++)
	{
		n = (c==0 ? Tab->item(r, Col)->text().toInt() : (Tab->item(r, Col)->text().toInt() - 1) / 10);
		Tab->item(r, Col)->setIcon(Iso->IsoImage[n]);
	}
}

void TableWidgetWindow::shiftCellValue(int v)
{
	if (Tab == 0)
	{
		printf(("TableWidgetWindow::shiftCellValue error: this function is not implemented for "
			   + getTypeString()).toLatin1().data());
		return;
	}
	int n, r, c;
	QList<QTableWidgetSelectionRange> L = Tab->selectedRanges();
	for (n=0; n < L.count(); n++)
		for (r = L[n].topRow(); r <= L[n].bottomRow(); r++)
			for (c = L[n].leftColumn(); c <= L[n].rightColumn(); c++)
				Tab->item(r, c)->setText(QString::number(Tab->item(r, c)->text().toInt() + v));
}

void TableWidgetWindow::tabItemChanged(QTableWidgetItem* Item)
{
	int n, N = Tab->columnCount(), r = Item->row();
	QString *IT = new QString[N];
	for (n=0; n<N; n++) IT[n] = (Tab->item(r, n) != 0 ? Tab->item(r, n)->text() : "");
	emit TabRowChanged(IT, N);
	delete[] IT;
}

void TableWidgetWindow::exportTableData(QString FileName, bool selectedCells, bool exchangeRowsColumns)
{
	if (Tab == 0) return;
	int n, r, c;
	bool VHI = (Tab->verticalHeaderItem(0) != 0 ?
					Tab->verticalHeaderItem(0)->text() != "1" : false);
	QFile F(FileName);
	F.open(QIODevice::WriteOnly);
	QTextStream S(&F);
	if (selectedCells)
	{
		QList<QTableWidgetSelectionRange> L = Tab->selectedRanges();
		if (exchangeRowsColumns) for (n=0; n < L.count(); n++)
		{
			if (VHI) for (r = L[n].topRow(); r <= L[n].bottomRow(); r++)
				S << '\t' << Tab->verticalHeaderItem(r)->text();
			S << '\n';
			for (c = L[n].leftColumn(); c <= L[n].rightColumn(); c++)
			{
				S << Tab->horizontalHeaderItem(c)->text() << '\t';
				for (r = L[n].topRow(); r <= L[n].bottomRow(); r++)
					S << (Tab->item(r, c) != 0 ? Tab->item(r, c)->text() : "")
					  << (r < L[n].bottomRow() ? '\t' : '\n');
			}
		}
		else for (n=0; n < L.count(); n++)
		{
			for (c = L[n].leftColumn(); c <= L[n].rightColumn(); c++)
				S << '\t' << Tab->horizontalHeaderItem(c)->text();
			S << '\n';
			for (r = L[n].topRow(); r <= L[n].bottomRow(); r++)
			{
				if (VHI) S << Tab->verticalHeaderItem(r)->text() << '\t';
				for (c = L[n].leftColumn(); c <= L[n].rightColumn(); c++)
					S << (Tab->item(r, c) != 0 ? Tab->item(r, c)->text() : "")
					  << (c < L[n].rightColumn() ? '\t' : '\n');
			}
		}
	}
	else
	{
		if (exchangeRowsColumns)
		{
			if (VHI) for (r=0; r < Tab->rowCount(); r++)
				S << '\t' << Tab->verticalHeaderItem(r)->text();
			S << '\n';
			for (c=0; c < Tab->columnCount(); c++)
			{
				S << Tab->horizontalHeaderItem(c)->text() << '\t';
				for (r=0; r < Tab->rowCount(); r++)
					S << (Tab->item(r, c) != 0 ? Tab->item(r, c)->text() : "")
					  << (r < Tab->rowCount() - 1 ? '\t' : '\n');
			}
		}
		else
		{
			for (c=0; c < Tab->columnCount(); c++)
				S << '\t' << Tab->horizontalHeaderItem(c)->text();
			S << '\n';
			for (r=0; r < Tab->rowCount(); r++)
			{
				if (VHI) S << Tab->verticalHeaderItem(r)->text() << '\t';
				for (c=0; c < Tab->columnCount(); c++)
					S << (Tab->item(r, c) != 0 ? Tab->item(r, c)->text() : "")
					  << (c < Tab->columnCount() - 1 ? '\t' : '\n');
			}
		}
	}
}

QString **TableWidgetWindow::getData(int &NR, int &NC)
{
	int r, c;
	if (NR <= 0 || NR > Tab->rowCount()) NR = Tab->rowCount();
	if (NC <= 0 || NC > Tab->columnCount()) NC = Tab->columnCount();
	QString **Data = new QString*[NR];
	for (r=0; r < NR; r++)
	{
		Data[r] = new QString[NC];
		for (c=0; c < NC; c++) if (Tab->item(r, c) != 0) Data[r][c] = Tab->item(r, c)->text();
	}
	return Data;
}

void TableWidgetWindow::setCellText(QString Text)
{
    QList<QTableWidgetSelectionRange> selectedRanges = Tab->selectedRanges();
    int n, r, c;
	Tab->blockSignals(true);
    for (n=0; n < selectedRanges.size(); ++n) for (r = selectedRanges[n].topRow(); r <= selectedRanges[n].bottomRow(); ++r)
        for (c = selectedRanges[n].leftColumn(); c <= selectedRanges[n].rightColumn(); ++c) Tab->item(r, c)->setText(Text);
	Tab->blockSignals(false);
}

void TableWidgetWindow::setData(QString **Data, int NR, int NC)
{
	//printf("TableWidgetWindow::setData\n");
	if (Typ == 0)
	{
		printf("Error: For TableWindows of typ 0 the function 'setData' is not implemented.");
		return;
	}
	int r, c;
	QTableWidgetItem *I;
	Tab->blockSignals(true);
	Tab->setRowCount(NR);
	Tab->setColumnCount(NC);
	for (r=0; r < NR; r++) for (c=0; c < NC; c++)
	{
		if ((I = Tab->item(r, c)) == 0) Tab->setItem(r, c, new QTableWidgetItem(Data[r][c]));
		else I->setText(Data[r][c]);
	}
	Tab->blockSignals(false);
	Changed();
}

void TableWidgetWindow::setEditable(bool Editable)
{
	int r, c, NR = Tab->rowCount(), NC = Tab->columnCount();
	if (Editable) for (r=0; r < NR; r++) for (c=0; c < NC; ++c)
		Tab->item(r, c)->setFlags(Tab->item(r, c)->flags() | Qt::ItemIsEditable);
	else for (r=0; r < NR; r++) for (c=0; c < NC; ++c)
		Tab->item(r, c)->setFlags(Tab->item(r, c)->flags() & (~Qt::ItemIsEditable));
}

void TableWidgetWindow::setHorizontalHeader(QStringList &Labels)
{
	if (Tab != 0)
	{
		if (Tab->columnCount() < Labels.count()) Tab->setColumnCount(Labels.count());
		Tab->setHorizontalHeaderLabels(Labels);
	}
	else printf("Error: For TableWindows of typ 0 the function 'setHorizontalHeader' is not implemented.");
}

void TableWidgetWindow::setVerticalHeader(QStringList &Labels)
{
	if (Tab != 0) Tab->setVerticalHeaderLabels(Labels);
	else printf("Error: For TableWindows of typ 0 the funtion 'setVerticalHeader' is not implemented.");
}

int *TableWidgetWindow::heapSort(bool sortFuncs(const QTableWidget *const, const int, const int)) const
{
    return heapsort::heapSort(TabSortFunctor(Tab, sortFuncs), getNumLines());
}

void TableWidgetWindow::shrinkAllSpectRefs(int FileColumn)
{
    int n, m;
    QString FileName;
    for (n=0; n < Tab->rowCount(); ++n) if (Tab->item(n, FileColumn) != 0)
    {
        FileName = Tab->item(n, FileColumn)->text();
        if ((m = FileName.lastIndexOf(QRegExp("[\\/]"))) >= 0) Tab->item(n, FileColumn)->setText(FileName.right(FileName.length() - m - 1));
    }
}

void TableWidgetWindow::sortTab(int* S2)
{
	int i, P1=0, n, N = Tab->rowCount(), c, C = Tab->columnCount();
	QString AIt[C][2];
	Tab->blockSignals(true);
	for (i=0; i<N; ++i) if (S2[i] != i)
	{
		for (c=0; c<C; ++c) AIt[c][P1] = Tab->item(i, c)->text();
		while (S2[i] != i)
		{
			for (c=0; c<C; ++c)
			{
				AIt[c][1-P1] = Tab->item(S2[i], c)->text();
				Tab->item(S2[i], c)->setText(AIt[c][P1]);
			}
			P1 = 1 - P1;
			n = S2[S2[i]];
			S2[S2[i]] = S2[i];
			S2[i] = n;
		}
		for (c=0; c<C; ++c) Tab->item(i, c)->setText(AIt[c][P1]);
	}
	Tab->blockSignals(false);
	delete[] S2;
	Changed();
}

QStringList TableWidgetWindow::getHorizontalHeaderLabels()
{
	int c, C = Tab->columnCount();
	QStringList R;
	for (c=0; c<C; c++) R << Tab->horizontalHeaderItem(c)->text();
	return R;
}

void TableWidgetWindow::search(int column, int value, int smeqla)
{
    int n, rStart = 0, N = Tab->rowCount();
    QList<QTableWidgetSelectionRange> Selected = Tab->selectedRanges();
    for (n=0; n < Selected.count(); ++n) if (Selected[n].bottomRow() > rStart) rStart = Selected[n].bottomRow();
	switch (smeqla)
	{
		case 0:
            for (n = rStart; (n<N && Tab->item(n, column)->text().toInt() != value) || (n==N && rStart > 0); n++) if (n==N && rStart > 0) n = rStart = -1;
			break;
		case 1:
            for (n = rStart; (n<N && Tab->item(n, column)->text().toInt() >= value) || (n==N && rStart > 0); n++) if (n==N && rStart > 0) n = rStart = -1;
			break;
		case 2:
            for (n = rStart; (n<N && Tab->item(n, column)->text().toInt() <= value) || (n==N && rStart > 0); n++) if (n==N && rStart > 0) n = rStart = -1;
			break;
		case 3:
            for (n = rStart; (n<N && abs(Tab->item(n, column)->text().toInt()) >= value) || (n==N && rStart > 0); n++) if (n==N && rStart > 0) n = rStart = -1;
			break;
		case 4:
            for (n = rStart; (n<N && abs(Tab->item(n, column)->text().toInt()) <= value) || (n==N && rStart > 0); n++) if (n==N && rStart > 0) n = rStart = -1;
			break;
	}
	if (n<N)
	{
		Tab->item(n, column)->setSelected(true);
		Tab->scrollToItem(Tab->item(n, column));
	}
	else QMessageBox::information(this, "MolSpektAnalysis", "A number fullfilling the condition could not be found.");
}

void TableWidgetWindow::search(int column, double value, int smeqla)
{
    int n, rStart = 0, N = Tab->rowCount();
    QList<QTableWidgetSelectionRange> Selected = Tab->selectedRanges();
    for (n=0; n < Selected.count(); ++n) if (Selected[n].bottomRow() > rStart) rStart = Selected[n].bottomRow();
	switch (smeqla)
	{
		case 0:
            for (n = rStart; (n<N && Tab->item(n, column)->text().toDouble() != value) || (n==N && rStart > 0); n++) if (n==N && rStart > 0) n = rStart = -1;
			break;
		case 1:
            for (n = rStart; (n<N && Tab->item(n, column)->text().toDouble() >= value) || (n==N && rStart > 0); n++) if (n==N && rStart > 0) n = rStart = -1;
			break;
		case 2:
            for (n = rStart; (n<N && Tab->item(n, column)->text().toDouble() <= value) || (n==N && rStart > 0); n++) if (n==N && rStart > 0) n = rStart = -1;
			break;
		case 3:
            for (n = rStart; (n<N && abs(Tab->item(n, column)->text().toDouble()) >= value) ||(n==N && rStart > 0); n++) if (n==N && rStart > 0) n = rStart = -1;
			break;
		case 4:
            for (n = rStart; (n<N && abs(Tab->item(n, column)->text().toDouble()) <= value) ||(n==N && rStart > 0); n++) if (n==N && rStart > 0) n = rStart = -1;
			break;
	}
	if (n<N)
	{
		Tab->item(n, column)->setSelected(true);
		Tab->scrollToItem(Tab->item(n, column));
	}
	else QMessageBox::information(this, "MolSpektAnalysis", "A number fullfilling the condition could not be found.");
}

void TableWidgetWindow::search(QString Text, int column, bool completeCell)
{
    int n, N = Tab->rowCount(), c, C = Tab->columnCount(), rStart = 0;
    QList<QTableWidgetSelectionRange> Selected = Tab->selectedRanges();
    for (n=0; n < Selected.count(); ++n) if (Selected[n].bottomRow() > rStart) rStart = Selected[n].bottomRow();
	if (column >= 0 && column < C)
	{
        if (completeCell)
        {
            for (n = rStart, c = column; (n<N && Tab->item(n, column)->text() != Text) || (n==N && rStart > 0); n++) if (n==N && rStart > 0) n = rStart = -1;
        }
        else for (n = rStart, c = column; (n<N && Tab->item(n, column)->text().indexOf(Text) == -1) || (n==N && rStart > 0); n++) if (n==N && rStart > 0) n = rStart = -1;
	}
	else
	{
        if (completeCell) for (n = rStart, c=C; (n<N && c==C) || (n==N && rStart > 0); n++)
        {
            if (n==N && rStart > 0) n = rStart = -1;
            else for (c=0; c<C && Tab->item(n, c)->text() != Text; n++) ;
        }
        else for (n = rStart, c=C; (n<N && c==C) || (n==N && rStart > 0); n++)
        {
            if (n==N && rStart > 0) n = rStart = -1;
            else for (c=0; c<C && Tab->item(n, c)->text().indexOf(Text) == -1; n++) ;
        }
	}
	if (n<N)
	{
		Tab->item(n, c)->setSelected(true);
		Tab->scrollToItem(Tab->item(n, c));
	}
	else QMessageBox::information(this, "MolSpektAnalysis", "The text \'" + Text + "\' could not be found in the table.");
}
