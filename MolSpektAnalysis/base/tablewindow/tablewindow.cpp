//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//


#include <QResizeEvent>
#include <QFile>
#include <QTextStream>
#include <QTableView>
#include <QBoxLayout>
#include <QGridLayout>
#include <QPainter>
#include <QMessageBox>
#include <QProgressBar>

#include <math.h>

#include "tablewindow.h"
#include "molecule.h"
#include "utils.h"
#include "Spektrum.h"
#include "heapsort.h"


TableWindow::TableWindow(Type typ, MainWindow *mw, Molecule *M) : MDIChild(typ, mw)
{
	//printf("TableWindow::TableWindow: typ=%d\n", typ);
	QLabel *NL = 0, *SL = 0;
	molecule = 0;
	setMolecule(M);
	Typ = typ;
	MW = mw;
	vMLabel = 0;
	vMax = 0;
	JMLabel = 0;
	JMax = 0;
	errLabel = 0;
	error = 0;
	Tab = 0;
	table = 0;
	Iso = 0;
	IsoLabel = 0;
	Comp = 0;
	CompLabel = 0;
	View = 0;
	ViewLabel = 0;
	Js = 0;
	JsLabel = 0;
	Jss = 0;
	JssLabel = 0;
	LSBox = 0;
	LSLabel = 0;
	MolBox = 0;
	MolLabel = 0;
	USBox = 0;
	USLabel = 0;
	Source = 0;
	Name = 0;
	JssB = JsB = 0;
	Pot1 = Pot2 = Date = 0;
	switch (typ)
	{
		case TermEnergyView:
			setMinimumSize(500, 300);
			break;
		case FranckCondonView:
			setMinimumSize(800, 600);
			break;
		default:
			setMinimumSize(300, 300);
			break;
	}
    if (Typ != FranckCondonView && Typ != TextTable1 && Typ != TextTable2 && Typ != FitSeriesResultTable && Typ != FranckCondonTable
            && Typ != External)
	{
		NL = new QLabel("Name: ", this);
		NL->setGeometry(10, 10, 50, 20);
		Name = new QLineEdit("", this);
		Name->setGeometry(60, 10, 80, 20);
		connect(Name, SIGNAL(editingFinished()), this, SLOT(setName()));
		
		SL = new QLabel("Source: ", this);
		SL->setGeometry(160, 10, 60, 20);
		Source = new QLineEdit("", this);
		Source->setGeometry(230, 10, width() - 240, 20);
		connect(Source, SIGNAL(editingFinished()), this, SLOT(sourceChanged()));
	}
	int w = width(), h = height();
	int w3 = w / 3, w2 = w / 2, w4 = w / 4;
	QGridLayout *L, *L1, *L2;
	QBoxLayout *Layout;
	switch (typ)
	{
		case FitDataSet:
			L = new QGridLayout(this);
			L1 = new QGridLayout;
			L1->addWidget(NL, 0, 0);
			L1->addWidget(Name, 0, 1);
			L1->addWidget(SL, 0, 2);
			L1->addWidget(Source, 0, 3);
			L1->setColumnStretch(3, 3);
			L1->setColumnStretch(1, 1);
			L->addLayout(L1, 0, 0, 1, 4);
			vMLabel = new QLabel("Max v:", this);
			L->addWidget(vMLabel, 1, 0);
			vMax = new QLineEdit(this);
			L->addWidget(vMax, 1, 1);
			JMLabel = new QLabel("Max J:", this);
			L->addWidget(JMLabel, 1, 2);
			JMax = new QLineEdit(this);
			L->addWidget(JMax, 1, 3);
			table = new MTable(this);
			L->addWidget(table, 2, 0, 1, 4);
			setFilter("Fit datasets (*.fdat)");
			setFileExt(".fdat");
			connect(vMax, SIGNAL(editingFinished()), this, SLOT(Changed()));
			connect(JMax, SIGNAL(editingFinished()), this, SLOT(Changed()));
			break;
		case TextTable1:
		case TextTable2: 
			setFilter("Text tables (*.dat)");
			setFileExt(".dat");
			Tab = new QTableWidget(this);
			Layout = new QBoxLayout(QBoxLayout::LeftToRight, this);
			Layout->addWidget(Tab); 
			break;
		case TermEnergyTable:
			table = new MTable(this);
			table->setGeometry(0, 40, w, h - 40);
			setFilter("Term energy tables (*.term)");
			setFileExt(".term");
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
		case LineTab:
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
			printf("TableWindow::TableWindow: Error: The type %d is not a valid type for a tablewindow!", typ);
			break;
	}
	if (Typ != TermEnergyTable && Typ != FitDataSet) 
	{
		connect(Tab, SIGNAL(itemChanged(QTableWidgetItem*)), this, SLOT(Changed()));
		connect(Tab, SIGNAL(itemChanged(QTableWidgetItem*)), this, SLOT(tabItemChanged(QTableWidgetItem*)));
	}
    if (Typ != -1 && Typ != -2 && Typ != 3 && Typ != 5 && Typ != External)
        connect(Source, SIGNAL(editingFinished()), this, SIGNAL(SourceChanged()));
	Saved();
}

TableWindow::~TableWindow()
{
	int n;
	for (n=0; n < ViewLists.count(); n++) delete[] ViewLists[n].ViewnRows;
}

void TableWindow::resizeEvent(QResizeEvent *E)
{
	//printf("Beginn TableWindow::resizeEvent\n");
	if (Typ <= -1 || Typ == 5) 
	{
		E->ignore();
		return;
	}
	QRect G;
	if (Typ <= 1 || Typ == 4)
	{
		G = Source->geometry();
		G.setWidth(width() - 240);
		Source->setGeometry(G);
	}
	//printf("ResizeEvent: Typ=%d\n", Typ);
	int w = width();
	int w3 = w / 3, w2 = w / 2, w4 = w / 4;
	switch(Typ)
	{
		case TermEnergyTable:
			G = table->geometry();
			G.setWidth(width());
			G.setHeight(height() - 40);
			table->setGeometry(G);
			break;
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
		case LineTab:
		case PotData:
			G = Tab->geometry();
			G.setWidth(width());
			G.setHeight(height() - 40);
			Tab->setGeometry(G);
			break;
		default:
			break;
	}
	E->ignore();
	//printf("Ende TableWindow::resizeEvent\n");
}

bool TableWindow::readData(QString Filename)
{
	//printf("TableWindow::readData\n");
	if (Typ == 0) 
	{
		printf("TableWindow::readData error: function not suited for Typ==0!\n");
		return false;
	}
	QFile Datei(Filename);
	if (!read(&Datei)) return false;
	int n, r, cc, lc = 0;
	QString Buffer;
	QTextStream S(&Datei);
	QString Spacer = (Typ > 0 ? " | " : "\t");
	QStringList L;
	if (Typ > 0 || Typ == -3)
	{
		if ((Buffer = S.readLine()).left(8) != "Source: ") 
		{
			if (Typ != -3 || Buffer.left(5).toInt() > 0) return false;
            else
            {
                Buffer = S.readLine();
                if (Buffer.left(5).toInt() > 0 && Buffer.mid(15, 5).toInt() > 0 && Buffer.mid(20, 5).toInt() > 0) return false;
                else setWindowTitle(Buffer);
            }
		}
		else
		{
			setSource(Buffer.right(Buffer.length() - 8));
			if ((Buffer = S.readLine()).left(6) != "Name: ") return false;
			setName(Buffer.right(Buffer.length() - 6));
			//printf("Nach setName\n");
		}
	}
	else setWindowTitle(S.readLine());
	if (Typ == -3) for (n=0; n<2; n++)
	{
		Buffer = S.readLine();
		if (Buffer.left(7) == "Max v: ") vMax->setText(Buffer.right(Buffer.length() - 7));
		else if (Buffer.left(7) == "Max J: ") JMax->setText(Buffer.right(Buffer.length() - 7));
	}
    QRegExp specialSectionStart = GetStartSpecialPartRegExp();
	bool Success = true;
	if (Typ == FitDataSet)
	{
		Buffer = readData(S);
		if (!Buffer.isEmpty()) Success = ReadSpecialPart(S, Buffer);
	}
	else
	{
		Tab->blockSignals(true);
		for (r=0, cc = Tab->columnCount(); !S.atEnd(); r++)
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
			for (n=0; n < lc; n++) Tab->setItem(r, n, new QTableWidgetItem(L[n]));
			while (n < cc) Tab->setItem(r, n++, new QTableWidgetItem(""));
		}
		if (Typ != 4 && lc < cc) r--;
		Tab->setRowCount(r);
		Tab->blockSignals(false);
	}
	Saved();
	//printf("Ende readData\n");
    return Success;
}

bool TableWindow::readData(QTextStream&)
{
	return false;
}

bool TableWindow::writeData(QString Filename)
{
	QFile Datei(Filename);
	if (!write(&Datei)) return false;
	int r, c, R = Tab->rowCount(), C = Tab->columnCount();
	QTextStream S(&Datei);
	QStringList L;
	QTableWidgetItem *I;
	QString Spacer = (Typ >= 0 && Typ != FitSeriesResultTable ? " | " : "\t");
	if ((Typ >= 0 && Typ != FitSeriesResultTable) || Typ == -3)
	{
		S << "Source: " << Source->text() << "\n";
		S << "Name: " << Name->text() << "\n";
	}
	else S << windowTitle() << "\n";
	if (Typ == -3)
	{
		S << "Max v: " << vMax->text() << "\n";
		S << "Max J: " << JMax->text() << "\n";
	}
	for (c=0; c<C; c++) L << ((I = Tab->horizontalHeaderItem(c)) != 0 ? I->text() : "?");
	S << "Column titles: " << L.join(Spacer) << "\n";
	if (Typ == FitDataSet) writeData(S);
	else for (r=0; r < R; r++)
	{
		L.clear();
		for (c=0; c < C; c++) L << ((I = Tab->item(r, c)) != 0 ? I->text() : "");
		S << L.join(Spacer).remove('\n') << "\n";
	}
	Saved();
	return true;
}

bool TableWindow::isAssigned()
{
	if (molecule == 0) return false;
	return true;
}

Molecule *TableWindow::getMolecule()
{
	return molecule;
}

QString TableWindow::getSource() const
{
	return Source->text();
}

int TableWindow::getvMax()
{
	QString T;
	if ((Typ == 1 || Typ == -3) && !(T = vMax->text()).isEmpty()) return T.toInt();
	return -1;
}

int TableWindow::getJMax()
{
	QString T;
	if ((Typ == 1 || Typ == -3) && !(T = JMax->text()).isEmpty()) return T.toInt();
	return -1;
}

void TableWindow::setJMax(int JM)
{
	if (JMax != 0) JMax->setText(QString::number(JM));
}

void TableWindow::setvMax(int vM)
{
	if (vMax != 0) vMax->setText(QString::number(vM));
}

void TableWindow::setViewnRows(MDIChild* Viewer, int NRows, int* Rows)
{
	int n;
	for (n=0; (n < ViewLists.count() ? Viewer != ViewLists[n].Viewer : false); n++) ;
	if (n == ViewLists.count())
	{
		ViewList B;
		B.Viewer = Viewer;
		ViewLists.append(B);
	}
	else delete[] ViewLists[n].ViewnRows;
	ViewLists[n].NumViewnRows = NRows;
	ViewLists[n].ViewnRows = Rows;
	emit SelChanged();
}

void TableWindow::getViewnE(int*& /*Js*/, double*& /*E*/, int& /*N*/)
{
	printf(("TableWindow::getViewnE: Error, this function is not implemented for " + getTypeString() + '!').toLatin1().data());
}

void TableWindow::getViewnRows(bool *RB)
{
	int N = getNumLines(), n, m;
	for (n=0; n<N; n++) RB[n] = false;
	for (n=0; n < ViewLists.count(); n++) for (m=0; m < ViewLists[n].NumViewnRows; m++)
		if (ViewLists[n].ViewnRows[m] >= 0 && ViewLists[n].ViewnRows[m] < N) RB[ViewLists[n].ViewnRows[m]] = true;
}

double TableWindow::getError()
{
	if (Typ == 1) return error->text().toDouble();
	return -1.0;
}

void TableWindow::setMolecule(Molecule *mol)
{
	if (mol != molecule)
	{
		molecule = mol;
		Changed();
	}
}

void TableWindow::setName()
{
	setName(Name->text());
}

void TableWindow::setName(QString name)
{
	//printf("TableWindow::setName\n");
	MDIChild::setName(name);
	name = getName();
	if (Name != 0) Name->setText(name);
	if (Filename.isEmpty()) setWindowTitle(name);
	if (molecule != 0) molecule->ChangeObjName(getType());
}

void TableWindow::setSource(QString source)
{
	if (source == Source->text()) return;
	Source->setText(source);
	Changed();
}

void TableWindow::sourceChanged()
{
	if (Source->text() != tSource)
	{
		tSource = Source->text();
		Changed();
	}
}

void TableWindow::vMaxChanged()
{
	if (vMax->text() != tvMax)
	{
		tvMax = vMax->text();
		Changed();
	}
}

void TableWindow::JMaxChanged()
{
	if (JMax->text() != tJMax)
	{
		tJMax = JMax->text();
		Changed();
	}
}

void TableWindow::errorChanged()
{
	if (error->text() != terror)
	{
		terror = error->text();
		Changed();
	}
}

void TableWindow::RemoveDoubled()
{
	printf(("The function RemoveDoubled is not implemented in " + getTypeString()
				+ "s!").toLatin1().data());
	return;
}

int TableWindow::getNumParIt()
{
	if (Typ != FitSeriesResultTable) return -1;
	return NumParFits->text().toInt();
}

void TableWindow::setMaxParFits(int Max)
{
	if (Typ == FitSeriesResultTable) NumParFits->setValidator(new QIntValidator(0, Max, NumParFits));
}

void TableWindow::setNumIterations(int Finished, int Max)
{
	if (Typ != FitSeriesResultTable) return;
	if (Max >= 0) Progress->setMaximum(Max);
	Progress->setValue(Finished);
}

void TableWindow::setNumParIt(int N)
{
	if (Typ == FitSeriesResultTable) NumParFits->setText(QString::number(N));
}

