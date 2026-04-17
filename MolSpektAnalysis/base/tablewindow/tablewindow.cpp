//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2025
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


TableWindow::TableWindow(Type typ, MainWindow *mw, Molecule *M) : MDIChild(typ, mw), Spacer(typ >= 0 && typ != FitSeriesResultTable ? " | " : "\t")
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
	resizeHelper(G);
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
	int n;
	QString Buffer;
	QTextStream S(&Datei);
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
	if (Typ == -3) for (n=0; n<2; ++n)
	{
		Buffer = S.readLine();
		if (Buffer.left(7) == "Max v: ") vMax->setText(Buffer.right(Buffer.length() - 7));
		else if (Buffer.left(7) == "Max J: ") JMax->setText(Buffer.right(Buffer.length() - 7));
	}
<<<<<<< HEAD
	bool Success = readData(S);
=======
	bool Success = true;
	if (Typ == FitDataSet)
	{
		Buffer = readData(S);
		if (!Buffer.isEmpty()) Success = ReadSpecialPart(S, Buffer);
	}
	else Success = readData(S);
>>>>>>> f91263c093dfe7215b3249af2e1113f12a7a6877
	Saved();
	//printf("Ende readData\n");
    return Success;
}

bool TableWindow::writeData(QString Filename)
{
	QFile Datei(Filename);
	if (!write(&Datei)) return false;
	QTextStream S(&Datei);
	QStringList L;
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
	if (Typ == FitDataSet) writeData(S);
	else writeData(S);
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

