//
// C++ Implementation: linetable
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2008 - 2019
//
// Copyright: See README file that comes with this source code
//
//


#include <QStringList>
#include <QString>
#include <QFile>
#include <QTextStream>
#include <QTableWidget>
#include <QTableWidgetItem>
#include <QTableWidgetSelectionRange>
#include <QList>
#include <QMessageBox>
#include <QPainter>
#include <QFileDialog>
#include <QPushButton>
#include <QComboBox>
#include <QCheckBox>

#include "linetable.h"
#include "elstate.h"
#include "molecule.h"
#include "duntable.h"
#include "MainWindow.h"
#include "Spektrum.h"
#include "tools.h"
#include "termtable.h"
#include "utils.h"
#include "potentialplot.h"
#include "potential.h"
#include "atom.h"
#include "tableline.h"
#include "isotab.h"
#include "termenergy.h"
#include "progression.h"
#include "line.h"
#include "questionbox.h"
#include "vsolistelement.h"
#include "matrix.h"

#include <math.h>


LineTable::LineTable(MainWindow *MW, Molecule *M, Transition *T) : TableWindow(LineTab, MW, M)
{	
	setFilter("Line table (*.lines)");
	setFileExt(".lines");
	transition = T;
	NR = 0;
	lRow = -1;
	Error = 0.005;
	OvError = 0.02;
	termTable = 0;
	mvs = mJs = 0;
	SelE = 0;
	SelJs = 0;
	NSel = 0;
	SO = 0;
	NSO = 0;
	MaxPN = NpProg = NpL = 0;
	if (M!=0) Iso = M->getIso();
	else Iso = 0;
	setName("newLineTable");
	setSource("own work");
	setWindowTitle("New line table");
	Tab->setColumnCount(TableNormCols);
	Tab->setRowCount(0);
	Tab->setHorizontalHeaderLabels(HeaderLabels << "PN" << "v'" << "J'" << "v''" << "J''" 
			<< "FC" << "wave number" << "error" << "isotope" << "file name" << "SNR" 
			<< "Deviation" << "Comment");
	Tab->setColumnWidth(CPN, 50);
	Tab->setColumnWidth(Cvs, 50);
	Tab->setColumnWidth(CJs, 50);
	Tab->setColumnWidth(Cvss, 50);
	Tab->setColumnWidth(CJss, 50);
	Tab->setColumnWidth(CF, 30);
	Tab->setColumnWidth(CWN, 150);
    Tab->setColumnWidth(CFile, 250);
	resize(600, 600);
	//connect(Tab, SIGNAL(itemChanged(QTableWidgetItem*)), this, SLOT(Changed()));
	connect(Tab, SIGNAL(itemSelectionChanged()), this, SLOT(TabSelChanged()));
	Saved();	
}

LineTable::~LineTable()
{
	//printf("LineTable::~LineTable\n");
	if (Iso != 0) delete Iso;
	if (SelE != 0)
	{
		delete[] SelE;
		delete[] SelJs;
	}
	if (SO != 0) delete[] SO;
	//printf("Ende ~LineTable\n");
}

void LineTable::setMolecule(Molecule *Mol)
{
	TableWindow::setMolecule(Mol);
	if (Iso != 0) delete Iso;
	if (Mol != 0) Iso = Mol->getIso();
	else Iso = 0;
	//printf("Ende von LineTable::setMolecule\n");
}

void LineTable::updateMarker(Spektrum *Spectrum)
{
	UpdateMarker(Spectrum, 0, 0, false);
}

void LineTable::UpdateMarker(Spektrum *Spectrum, int nLines, int *Lines, bool remove)
{
	//printf("LineTable::UpdateMarker()\n");
	Marker *marker; 
	int i, j, nr = Tab->rowCount(), s, *n, N, iB = 0, **p, *lines, I, AM;
	if (nr < 2) return;
	bool L = false;
	double dB, **F;
	QString T, lSN, uSN;
	Spektrum **Spectra;
	ElState *lState, *uState;
	if (transition != 0)
	{
		if ((lState = transition->getLowerState()) != 0) lSN = lState->getName();
		if ((uState = transition->getUpperState()) != 0) uSN = uState->getName();
	}
	if (Spectrum != 0) 
	{
		Spectra = new Spektrum*[N=1];
		Spectra[0] = Spectrum;
	}
	else if (MW != 0)
	{
		N = MW->getNumSpectra();
		Spectra = new Spektrum*[N];
		for (s=0; s<N; s++) Spectra[s] = MW->getSpectrum(s);
	}
	else return;
	QString SpektFile[N];
	p = new int*[N];
	F = new double*[N];
	n = new int[N];
	for (s=0; s<N; s++) 
	{
		SpektFile[s] = Spectra[s]->getFName();
		if ((i = SpektFile[s].lastIndexOf('.')) > -1)
			SpektFile[s] = SpektFile[s].left(SpektFile[s][i-1] == 'm' ? i-1 : i);
		n[s] = 0;
	}
	if (Lines == 0) 
	{
		lines = new int[nLines = nr];
		for (i=0; i<nr; i++) lines[i] = i;
	}
	else lines = Lines;
	for (i=0; i<nLines; i++) 
	{
		T = Tab->item(lines[i], CFile)->text();
		//if (T == "" && L) n[s]++;
		for (s=0; s<N; s++) if (T.indexOf(SpektFile[s]) > -1) n[s]++;
	}
	for (s=0; s<N; s++) 
	{
		p[s] = new int[n[s]];
		F[s] = new double[n[s]];
	}
	//printf("N=%d, n[0]=%d, l[0]=%d, SpektFile[0]=%s\n", N, n[0], l[0], SpektFile[0].ascii());
	for (s=0; s<N; s++) n[s]=0;
	for (i=0, L = false; i < nLines; i++)
	{
		T = Tab->item(lines[i], CFile)->text();
		if (T != "" || !L) 
		{
			for (s=0, L = false; (s<N ? T.indexOf(SpektFile[s]) == -1 : false); s++) ;
			if (s<N) L = true;
		}
		if (L)
		{
			//printf("i=%d, n[0]=%d\n", lines[i], n[0]);
			F[s][n[s]] = Tab->item(lines[i], CWN)->text().toDouble();
			p[s][n[s]++] = lines[i];
		}
	}
	//printf("Vor letzter Schleife\n");
	for (s=0, iB=1; s<N; s++) 
	{
		while (iB != -1) for (i=1, iB = -1; i < n[s]; i++) 
				if (F[s][i-1] > F[s][i])
		{
			dB = F[s][i-1];
			F[s][i-1] = F[s][i];
			F[s][i] = dB;
			iB = p[s][i];
			p[s][i] = p[s][i-1];
			p[s][i-1] = iB;
		}
		Spectra[s]->GetMarker(AM, marker);
		if (AM == 0) continue;
		for (i=0, j=0; i<n[s]; i++) 
		{
			while(j < AM ? marker[j].Line[0] < F[s][i] : false) j++;
			if (j == AM) j--;
			if (j > 0 ? marker[j].Line[0] - F[s][i] > F[s][i] - marker[j-1].Line[0] : false) j--;
			if (!marker[j].Marked)
			{
				if (remove)
				{
					I = (Tab->item(p[s][i], CIso)->text().toInt() - 1)/10;
					if (I < 0 || I >= Iso->numIso) I = 0;
					if (marker[j].IsoName == Iso->texName[I]) marker[j].DisplayData = false;
				}
				else
				{
					//printf("i=%d, j=%d\n", i, j);
					//printf("p[%d][%d]=%d, Iso->numIso=%d\n", s, i, p[s][i], Iso->numIso);
					marker[j].vs = Tab->item(p[s][i], Cvs)->text().toInt();
					marker[j].Js = Tab->item(p[s][i], CJs)->text().toInt();
					marker[j].vss = Tab->item(p[s][i], Cvss)->text().toInt();
					marker[j].Jss = Tab->item(p[s][i], CJss)->text().toInt();
					marker[j].Iso = (Tab->item(p[s][i], CIso)->text().toInt() - 1)/10;
					marker[j].DD = Tab->item(p[s][i], CDev)->text().toDouble();
					marker[j].FC = Tab->item(p[s][i], CF)->text().toInt();
					marker[j].Mol = molecule;
					if (Iso != 0) 
					{
						if (marker[j].Iso < 0 || marker[j].Iso >= Iso->numIso ) marker[j].Iso = 0;
						marker[j].IsoName = Iso->texName[marker[j].Iso];
					}
					else marker[j].IsoName = "";
					//printf("Vor states\n");
					if (transition != 0) 
					{
						if ((marker[j].LState = transition->getLowerState()) != 0) 
							marker[j].lState = marker[j].LState->getName();
						else marker[j].lState = "";
						if ((marker[j].UState = transition->getUpperState()) != 0) 
							marker[j].uState = marker[j].UState->getName();
						else marker[j].uState = "";
					}
					else marker[j].lState = marker[j].uState = "";
					marker[j].DisplayData = true;
					//printf("Nach markieren\n");
				}
			}
		}
	}
	//printf("vor delete\n");
	delete[] Spectra;
	Destroy(p, N);
	Destroy(F, N);
	delete[] n;
	if (Lines == 0) delete[] lines;
	//printf("Ende UpdateMarker\n");
}

void LineTable::setTransition(Transition *T)
{
	transition = T;
}

Transition *LineTable::getTransition()
{
	return transition;
}

int LineTable::getAnzahlLinien()
{
    int N;
	//printf("LineTable::getAnzahlLinien: name = %s\n", getName().toAscii().data());
    for (N=0; (N < Tab->rowCount() ? Tab->item(N, 0)->icon().isNull() || 
		 Tab->item(N, Cvs)->text().isEmpty() : false); N++) ;// printf("N=%d\n", N);
    return N;
}

void LineTable::getLines(int **Zuordnung, double *Energien, double *Unc)
{
	//printf("LineTable::getLines()\n");
	int N;
    for (N=0; (N < Tab->rowCount() ? Tab->item(N, 0)->icon().isNull() : false); N++)
    {
		Zuordnung[N][0] = (Tab->item(N, CIso)->text().toInt() - 1) / 10;
		Zuordnung[N][1] = Tab->item(N, Cvs)->text().toInt();
		Zuordnung[N][2] = Tab->item(N, CJs)->text().toInt();
		Zuordnung[N][3] = Tab->item(N, Cvss)->text().toInt();
		Zuordnung[N][4] = Tab->item(N, CJss)->text().toInt();
		Zuordnung[N][5] = Tab->item(N, CF)->text().toInt();
		Energien[N] = Tab->item(N, CWN)->text().toDouble();
		Unc[N] = Tab->item(N, Cerr)->text().toDouble();
    }
	//printf("Ende von getLines\n");
}

void LineTable::getLines(const QString &Filename, double **Lines, int *numLines)
{
	int i, nr = Tab->rowCount(), n=0, l=Filename.length();
	bool L = false;
	QString T;
	for (i=0; i<nr; i++) 
	{
		T = Tab->item(i, 7)->text();
		if ((T.isEmpty() && L) || T.right(l) == Filename)
		{
			n++;
			L = true;
		}
		else L = false;
	}
	double *R = new double[n];
	//printf("getLines:n=%d\n", n);
	L = false;
	for (i=0, n=0; i<nr; i++)
	{
		T = Tab->item(i, 7)->text();
		if ((T.isEmpty() && L) || T.right(l) == Filename)
		{
			R[n++] = Tab->item(i, 4)->text().toDouble();
			L = true;
		}
		else L = false;
	}
	//printf("n=%d\n", n);
	*Lines = R;
	*numLines = n;
}

void LineTable::getLines(TableLine*& L, int& N)
{
	int r;
	if ((N = Tab->rowCount()) > 0) L = new TableLine[N];
	for (r=0; r<N; r++)
	{
		L[r].dev = Tab->item(r, CDev)->text().toDouble();
		L[r].err = Tab->item(r, Cerr)->text().toDouble();
		L[r].DevR = L[r].dev / L[r].err;
		L[r].FC = Tab->item(r, CF)->text().toInt();
		L[r].File = Tab->item(r, CFile)->text();
		L[r].Iso = (Tab->item(r, CIso)->text().toInt() - 1) / 10;
		L[r].isTE = false;
		L[r].Js = Tab->item(r, CJs)->text().toInt();
		L[r].Jss = Tab->item(r, CJss)->text().toInt();
		L[r].LTab = this;
		L[r].PN = Tab->item(r, CPN)->text().toInt();
		L[r].Row = r;
		L[r].SourceN = 0;
		L[r].vs = Tab->item(r, Cvs)->text().toInt();
		L[r].vss = Tab->item(r, Cvss)->text().toInt();
		L[r].WN = Tab->item(r, CWN)->text().toDouble();
	}
}

void LineTable::getSortedLines(TableLine*& L, int& N, int SortOrder)
{
	if (SortOrder != 0 || (N = Tab->rowCount()) == 0)
	{
		N=0;
		L=0;
		return;
	}
	int r, *SO = heapSort(sortByFrequency);
	L = new TableLine[N];
	for (r=0; r<N; r++)
	{
		L[SO[r]].dev = Tab->item(r, CDev)->text().toDouble();
		L[SO[r]].err = Tab->item(r, Cerr)->text().toDouble();
		L[SO[r]].DevR = L[r].dev / L[r].err;
		L[SO[r]].FC = Tab->item(r, CF)->text().toInt();
		L[SO[r]].File = Tab->item(r, CFile)->text();
		L[SO[r]].Iso = (Tab->item(r, CIso)->text().toInt() - 1) / 10;
		L[SO[r]].isTE = false;
		L[SO[r]].Js = Tab->item(r, CJs)->text().toInt();
		L[SO[r]].Jss = Tab->item(r, CJss)->text().toInt();
		L[SO[r]].LTab = this;
		L[SO[r]].PN = Tab->item(r, CPN)->text().toInt();
		L[SO[r]].Row = r;
		L[SO[r]].SourceN = 0;
		L[SO[r]].vs = Tab->item(r, Cvs)->text().toInt();
		L[SO[r]].vss = Tab->item(r, Cvss)->text().toInt();
		L[SO[r]].WN = Tab->item(r, CWN)->text().toDouble();
		L[SO[r]].SNR = Tab->item(r, CSNR)->text().toDouble();
	}
	delete[] SO;
}

int LineTable::getNgTE(int *mv, int mJ)
{
	int n, I, N, P, p, NR = Tab->rowCount(), v, J;
	if (NSO != NR) ShowUpTerm();
	for (n=N=0, P=-1; n < NR; n++) if ((p = Tab->item(SO[n], CPN)->text().toInt()) != P)
	{
		if (Tab->item(SO[n], CEUp)->text().toDouble() == 0.0) continue;
		v = Tab->item(SO[n], Cvs)->text().toInt();
		J = Tab->item(SO[n], CJs)->text().toInt();
		if (v >= 0 && J >= 0 ? (mv != 0 ? (J <= mJ ? v > mv[J] : true) : false) : true) 
			continue;
		I = Tab->item(SO[n], CIso)->text().toInt();
		while (10 * (I / 10) != I && n < NR) I = Tab->item(SO[n++], CIso)->text().toInt();
		if (n == NR) break;
		P = p;
		N++;
	}
	return N;
}

void LineTable::getgoodTE(int &N, TermEnergy *&E, int *mv, int mJ)
{
	//printf("LineTable::getgoodTE\n");
	int n, I, P = -1, p, NR = Tab->rowCount(), J, v;
	if (NSO != NR) ShowUpTerm();
	int gli[NR][2];
	//printf("Vor Schleife1\n");
	if (NSO != NR)
	{
		N = mJ = 0;
		mv = 0;
		E = 0;
		return;
	}
	for (n=N=0; n < NR; n++) if ((p = Tab->item(SO[n], CPN)->text().toInt()) != P)
	{
		//printf("N=%d, n=%d\n", N, n);
		if (Tab->item(SO[n], CEUp)->text().toDouble() == 0.0) continue;
		v = Tab->item(SO[n], Cvs)->text().toInt();
		J = Tab->item(SO[n], CJs)->text().toInt();
		if (mv != 0 ? (v >= 0 && J >= 0 ? (J <= mJ ? v > mv[J] : true) : true) : J < 0) 
			continue;
		I = Tab->item(SO[n], CIso)->text().toInt();
		while (10 * (gli[N][0] = I / 10) != I && n < NR) 
			I = Tab->item(SO[n++], CIso)->text().toInt();
		if (n == NR) break;
		P = p;
		gli[N][1] = SO[n];
		gli[N++][0]--;
	}
	//printf("Vor Schleife2\n");
	E = new TermEnergy[N];
	for (n=0; n<N; n++)
	{
		E[n].v = Tab->item(gli[n][1], Cvs)->text().toInt();
		E[n].J = Tab->item(gli[n][1], CJs)->text().toInt();
		E[n].Iso = gli[n][0];
		E[n].err = 0.01;
		E[n].E = Tab->item(gli[n][1], CEav)->text().toDouble();
		E[n].PN = Tab->item(gli[n][1], CPN)->text().toInt();
		E[n].File = Tab->item(gli[n][1], CFile)->text();
		E[n].ef = Tab->item(gli[n][1], CJss)->text().toInt() != E[n].J;
		E[n].dev = E[n].DevR = 0.0;
		E[n].FC = Tab->item(gli[n][1], CF)->text().toInt();
	}
	//printf("Ende getgoodTE");
}

int LineTable::getNgL(int *mv, int mJ)
{
	int n, I, N, NR = Tab->rowCount(), J;
	for (n=N=0; n < NR; n++)
	{
		if (mv != 0) if ((J = Tab->item(n, CJss)->text().toInt()) <= mJ) 
				if (Tab->item(n, Cvss)->text().toInt() > mv[J]) continue;
		I = Tab->item(n, CIso)->text().toInt();
		if (10 * (I / 10) == I) N++;
	}
	return N;
}

void LineTable::getgoodLines(int &N, TableLine *&L, int *mv, int mJ, bool SortFunction(const QTableWidget *const, const int, const int))
{
	//printf("Beginn LineTable::getgoodLines\n");
	int n, i, NR = Tab->rowCount(), J;
	//printf("NR=%d\n", NR);
    int *SO = new int[NR], *S1 = heapSort(SortFunction);
	int gli[NR][2];
	for (n=0; n < NR; n++) SO[S1[n]] = n;
	for (n = N = 0; n < NR; n++)
	{
		if (mv != 0) if ((J = Tab->item(SO[n], CJss)->text().toInt()) <= mJ) 
				if (Tab->item(SO[n], Cvss)->text().toInt() > mv[J]) continue;
		//printf("n=%d, NR=%d\n", n, NR);
		i = Tab->item(SO[n], CIso)->text().toInt();
		if (10 * (gli[N][0] = i / 10) == i) 
		{
			gli[N][1] = SO[n];
			gli[N++][0]--;
		}
	}
	L = new TableLine[N];
	//printf("n=%d, N=%d\n", n, N);
	for (n=0; n<N; n++)
	{
		L[n].FC = Tab->item(gli[n][1], CF)->text().toInt();
		L[n].PN = Tab->item(gli[n][1], CPN)->text().toInt();
		L[n].vs = Tab->item(gli[n][1], Cvs)->text().toInt();
		L[n].Js = Tab->item(gli[n][1], CJs)->text().toInt();
		L[n].vss = Tab->item(gli[n][1], Cvss)->text().toInt();
		L[n].Jss = Tab->item(gli[n][1], CJss)->text().toInt();
		L[n].Iso = gli[n][0];
		L[n].WN = Tab->item(gli[n][1], CWN)->text().toDouble();
		L[n].err = Tab->item(gli[n][1], Cerr)->text().toDouble();
		L[n].File = Tab->item(gli[n][1], CFile)->text();
		L[n].Row = gli[n][1];
	}
	delete[] SO;
	delete[] S1;
}

int LineTable::getMaxvs()
{
	int i, n=0, nr = Tab->rowCount();
	for (i=0; i<nr; i++) if (Tab->item(i, Cvs)->text().toInt() > n) n = Tab->item(i, Cvs)->text().toInt();
	return n;
}

int LineTable::getMaxvss()
{
	int i, n=0, nr = Tab->rowCount();
	for (i=0; i < nr; i++) if (Tab->item(i, Cvss)->text().toInt() > n) n = Tab->item(i, Cvss)->text().toInt();
	return n;
}

int LineTable::getMaxJs()
{
	int i, n=0, nr = Tab->rowCount();
	for (i=0; i < nr; i++) if (Tab->item(i, CJs)->text().toInt() > n) n = Tab->item(i, CJs)->text().toInt();
	return n;
}

int LineTable::getMaxJss()
{
	int i, n=0, nr = Tab->rowCount();
	for (i=0; i < nr; i++) if (Tab->item(i, CJss)->text().toInt() > n) n = Tab->item(i, CJss)->text().toInt();
	return n;
}

void LineTable::getObsIso(bool *O, int N)
{
	int I, r, nr = Tab->rowCount();
	for (r=0; r<N; r++) O[r] = false;
	for (r=0; r < nr; r++) if ((I = Tab->item(r, CIso)->text().toInt() / 10 - 1) < N) if (I >= 0) 
				O[I] = true;
}

void LineTable::getProgressions(int &N, Progression *&P, bool SortFunction(const QTableWidget *const, const int, const int))
{
	//printf("LineTable::getProgressions\n");
	TableLine *L;
	QString S, aS;
	int vs = -1, Js = -1, Iso = -1, n, m, NP=0, NR, PN = -1;
    getgoodLines(NR, L, 0, 0, SortFunction);
	int lpb[NR];
	for (n=0; n < NR; n++) 
		if (L[n].vs != vs || L[n].Js != Js || L[n].Iso != Iso || L[n].File != S 
			|| L[n].PN != PN) 
	{
		//printf("PN=%d, N=%d\n", PN, NP);
		vs = L[n].vs;
		Js = L[n].Js;
		S = L[n].File;
		Iso = L[n].Iso;
		PN = L[n].PN;
		lpb[NP++] = n;
	}
	lpb[NP] = NR;
	P = new Progression[N = NP];
	for (n=0; n<N; n++)
	{
		P[n].N = lpb[n+1] - lpb[n];
		P[n].vs = L[lpb[n]].vs;
		P[n].Js = L[lpb[n]].Js;
		P[n].Iso = L[lpb[n]].Iso;
		P[n].L = new Line[P[n].N];
        P[n].PNum = L[lpb[n]].PN;
		//printf("N=%d, P[%d].N=%d\n", N, n, P[n].N);
		for (m=0; m < P[n].N; m++)
		{
			P[n].L[m].vss = L[lpb[n] + m].vss;
			P[n].L[m].Jss = L[lpb[n] + m].Jss;
			P[n].L[m].E = L[lpb[n] + m].WN;
			P[n].L[m].err = L[lpb[n] + m].err;
			P[n].L[m].row = L[lpb[n] + m].Row;
		}
	}
	//printf("Ende getProgressions\n");
}

void LineTable::setUncertainty(int* RowNumbers, double* NewUncertainty, int NLines)
{
	int n;
	Tab->blockSignals(true);
	for (n=0; n < NLines; n++) Tab->item(RowNumbers[n], Cerr)->setText(QString::number(NewUncertainty[n], 'f', 4));
	Tab->blockSignals(false);
	Changed();
}

bool LineTable::readData(QString Filename)
{
    //printf("LineTable::readData\n");
	int n, s, i, N = 0, d=0;
	bool FCA = false, Success = true, COK;
    QPixmap P;
    QFile Datei(Filename);
    if (!read(&Datei)) return false;
    QTextStream S1(&Datei);
    QString B, Comm, Buffer, B2 = S1.readLine();
	QStringList L;
	MaxPN = 0;
	if (B2.indexOf("Source", 0, Qt::CaseInsensitive) != -1)
	{
		for (n = B2.indexOf(":") + 1; B2[n].isSpace(); n++) ;
		setSource(B2.right(B2.length() - n));
		B2 = S1.readLine();
		if (B2.indexOf("Name", 0, Qt::CaseInsensitive) != -1)
		{
			for (n = B2.indexOf(":") + 1; B2[n].isSpace(); n++) ;
			setName(B2.right(B2.length() - n));
		}
		B2 = S1.readLine();
		if (B2.indexOf("FC") > 0) FCA = true;
	}
	else 
	{
		setName("importedLineTable");
		d=1;
		if (B2.indexOf("F1/F2") > 0) FCA = true;
	}
    Buffer = S1.readAll();
    QTextStream S(&Buffer, QIODevice::ReadOnly);
	while (!S.atEnd())
    {
		S.readLine();
		N++;
		//printf("N=%d\n", N);
    }
    Tab->blockSignals(true);
    Tab->setRowCount(N);
    QTextStream S2(&Buffer, QIODevice::ReadOnly);
	//printf("Ende: N=%d\n", N);
    for (n=0; n<N; n++)
    {
		for (i=0; i < TableNormCols; i++) if (Tab->item(n, i) == 0) 
				Tab->setItem(n, i, new QTableWidgetItem(""));
		B = S2.readLine();
		L = B.split(QRegExp("\\s+"), QString::SkipEmptyParts);
		if ((s = L.size()) < 8) 
		{
			for (i=0; i < TableNormCols; i++) Tab->item(n, i)->setText("");
			printf("Fehler beim Lesen von %s, Zeile %d ist zu kurz!\n", Filename.toLatin1().data(), n);
			Success = false;
			continue;
		}
		//printf("d=%d, n=%d, s=%d, NC=%d\n", d, n, s, TableNormCols);
		if (FCA) 
		{
			if (s < CC)
			{
				if (s==8)
				{
					Tab->item(n, Cvs)->setText(L[1]);
					Tab->item(n, CJs)->setText(L[2]);
					Tab->item(n, Cvss)->setText(L[3]);
					Tab->item(n, CJss)->setText(L[4]);
					Tab->item(n, CWN)->setText(L[5]);
					Tab->item(n, Cerr)->setText(L[6]);
					Tab->item(n, CF)->setText(QString::number(L[7].toInt() + 1));
					Tab->item(n, CIso)->setText("10");
					continue;
				}
				else if (s==11)
				{
					for (i=0; i<9; i++) Tab->item(n, i)->setText(L[i]);
					continue;
				}
				else
				{
					printf("Error, line %d is too short!", n);
					Success = false;
					continue;
				}
			}
			if ((i = L[0].toInt()) > MaxPN) MaxPN = i;
			for (i=0; i < CSNR; i++) Tab->item(n, i)->setText(L[i]);
			for (L[i].toDouble(&COK); !COK && !L[i].isEmpty(); L[++i].toDouble(&COK))
			{
				L[CFile] += ' ' + L[i];
				if (i==s-1) 
				{
					i++;
					break;
				}
			}
			if (i > CSNR) Tab->item(n, CFile)->setText(L[CFile]);
			if (COK) Tab->item(n, CSNR)->setText(L[i++]);
			if (i<s ? L[i].toDouble() != 0.0 : false) Tab->item(n, CDev)->setText(L[i++]);
			if (s >= TableNormCols && i<s) Tab->item(n, CC)->setText(L[i]);
		}
		else
		{
			if (d==0) if ((i = L[0].toInt()) > MaxPN) MaxPN = i;
			for (i=d; i<=4; i++) Tab->item(n, i)->setText(L[i]);
			for (i=5; i<=7; i++) Tab->item(n, i+1)->setText(L[i]);
			Comm = "";
            for (i=8; (i < s ? L[i] == "laser" || L[i] == "?" || L[i] == "weak" || L[i] == "strong"
						|| L[i] == "overlap" : false); i++)
			{
				if (i==8) Comm = L[i];
				else Comm += " " + L[i];
			}
			Tab->item(n, CC)->setText(Comm);
			Tab->item(n, CF)->setText("-1");
			if ((i < s ? L[i] != "0" : false)) Tab->item(n, CFile)->setText(L[i]);
			else Tab->item(n, CFile)->setText("");
			if (i == 8 && 9 < s) 
			{
				if (L[9] == "laser") 
				{
					Tab->item(n, CC)->setText(Comm = L[9]);
					Tab->item(n, CSNR)->setText("");
				}
				else if (L[9] == "0") Tab->item(n, CSNR)->setText("");
				else Tab->item(n, CSNR)->setText(L[9]);
				if (s > 10) 
				{
					if (L[10] == "laser") Tab->item(n, CC)->setText(Comm = L[10]);
					else if (L[10] == "0") Tab->item(n, CDev)->setText("");
					else Tab->item(n, CDev)->setText(L[10]);
					if (s > 11) Tab->item(n, CC)->setText(L[11]);
					else if (Comm.isEmpty()) Tab->item(n, CC)->setText("");
				}
				else 
				{
					Tab->item(n, CDev)->setText("");
					if (Comm.isEmpty()) Tab->item(n, CC)->setText("");
				}
			}
			else 
			{
				Tab->item(n, CSNR)->setText("");
				Tab->item(n, CDev)->setText("");
				if (Comm.isEmpty()) Tab->item(n, CC)->setText("");
			}
		}
		for (i=0; i < TableNormCols; i++) Tab->item(n, i)->setIcon(P);
		//printf("n=%d, N=%d\n", n, N);
        if (molecule != 0)
        {
            QString CurPath = Tab->item(n, CFile)->text(), MolPath = molecule->getFileName();
            Tab->item(n, CFile)->setText(getAbsolutePath(CurPath, MolPath));
        }
    }
    Tab->blockSignals(false);
	Saved();
	emit DataChanged();
	return Success;
}

bool LineTable::writeData(QString Filename)
{
    int i, j;
    QString Buffer, LFile, B, MolFile = molecule->getFileName();
	QPixmap P;
	QFile Datei(Filename);
	write(&Datei);
    QTextStream S(&Datei);
	S << "Source: " << getSource() << "\n";
	S << "Name: " << getName() << "\n";
	S << "   PN  vo  Jo  vu  Ju  FC           Eout       err  Iso  File  SNR  deviation  comment\n";
    Tab->blockSignals(true);
	for (i=0; i<Tab->rowCount(); i++)
    {
		S << ((B = Tab->item(i, CPN)->text()).toInt() > 0 ? ("     " + B).right(5) 
														  : "    0");
		for (j=Cvs; j<=CJss; j++) 
		{
			Buffer = "    " + Tab->item(i, j)->text();
			S << Buffer.right(4);
		}
		S << ' ' << ("   " + Tab->item(i, CF)->text()).right(3);
		Buffer = "               " + Tab->item(i, CWN)->text();
		if (Buffer.indexOf('.') == -1) Buffer += ".0000";
		Buffer = Buffer.right(15);
		if (Buffer[0] != ' ') S << " " << Buffer;
		else S << Buffer;
		if ((Buffer = "          " + Tab->item(i, Cerr)->text()).indexOf('.') == -1) Buffer += ".000";
		S << ' ' << Buffer.right(9);
		S << ("    " + Tab->item(i, CIso)->text()).right(4); 
		if ((Buffer = Tab->item(i, CFile)->text()).isEmpty()) Buffer = LFile;
		else if (Buffer == "Aktuelle Markierungen") 
		{
		    Buffer = InFile;
			Tab->item(i, 7)->setText(Buffer);
		}
		else LFile = Buffer;
        if (molecule != 0) Buffer = getRelativePath(Buffer, MolFile);
		S << " " << Buffer;
		if ((Buffer = Tab->item(i, CSNR)->text()).isEmpty()) Buffer = "0";
		S << ("      " + Buffer).right(7); 
		if ((Buffer = Tab->item(i, CDev)->text()).isEmpty()) Buffer = "0";
		S << ("          " + Buffer).right(11);
		S << " " << Tab->item(i, CC)->text() << char(10);
	}
	for (i = Tab->rowCount() - NpL; i < Tab->rowCount(); i++)
		for (j=0; j < TableNormCols; j++) Tab->item(i, j)->setIcon(P);
    Datei.close();
	Tab->blockSignals(false);
	MaxPN += NpProg;
	NpProg = 0;
    Saved();
	NpL = 0;
	if (MW != 0 && transition != 0) 
		MW->LineTableSaved(transition->getUpperState(), transition->getLowerState());
	emit DataChanged();
	return true;
}

bool LineTable::writeExcPotFitInput(QString Filename)
{
	if (Filename.isEmpty()) return false;
	if (molecule == 0) return false;
	if (ShowUpTermTable() == 0) return false;
	IsoTab *Iso = molecule->getIso();
	if (Iso == 0) return false;
	int i, n = 0, m, ND = 0, NFD = 0, nIso = Iso->numIso, uIA[3], uIB[3], use = 0;
	QString **Data = termTable->getData(ND, n);
	ElState *St = (transition != 0 ? transition->getUpperState() : 0);
	int lambda = (St != 0 ? St->getLambda() : -1), Maxv = 0, MaxJ = 0; 
	Atom *A1 = molecule->getAtom1(), *A2 = molecule->getAtom2();
	int nIsoA = A1->getnIso(), nIsoB = A2->getnIso(), iA, iB, p1, p2;
	int IsoCA[nIsoA], IsoCB[nIsoB], IsoTA[nIso], IsoTB[nIso], nIA = 0, nIB = 0;
	QString IsoSA[3], IsoSB[3], IsoSAF[3], IsoSBF[3];
	bool Diff;
	double W, WSum, DS, *OUnc = new double[ND];
	if (A1 == A2)
		for (n=i=0; n < nIsoA; n++) for (m=n; m < nIsoA; m++)
	{
		IsoTA[i] = n;
		IsoTB[i++] = m;
	}
	else for (n=i=0; n < nIsoA; n++) for (m=0; m < nIsoB; m++)
	{
		IsoTA[i] = n;
		IsoTB[i++] = m;
	}
	for (n=0; n < nIsoA; n++) IsoCA[n] = 0;
	for (n=0; n < nIsoB; n++) IsoCB[n] = 0;
	TermEnergy *D = new TermEnergy[ND], *FD = 0, TBuff;
	for (n=m=0; n < ND; n++)
	{
		if (Data[n][4] == "nan") continue;
		D[m].Iso = Data[n][0].toInt() / 10 - 1;
		if ((D[m].v = Data[n][1].toInt()) > Maxv) Maxv = D[n].v;
		if ((D[m].J = Data[n][2].toInt()) > MaxJ) MaxJ = D[n].J;
		D[m].ef = (Data[n][3].indexOf('e') != -1 ? true : false);
		D[m].E = Data[n][4].toDouble();
		if (D[m].E == 0.0) continue;
		D[m].err = 0.01;
		OUnc[m] = Data[n][5].toDouble();
		IsoCA[IsoTA[D[m].Iso]]++;
		IsoCB[IsoTB[D[m++].Iso]]++;
	}
	ND = m;
	for (n=0; n < ND; n++)
	{
		for (m=n+1; D[m].Iso == D[n].Iso && D[n].v == D[m].v && D[n].J == D[m].J
			        && D[n].ef == D[m].ef; m++) ;
		if (m==n+1) continue;
		for (i=n, WSum = DS = 0.0; i<m; i++)
		{
			WSum += (W = ((W = OUnc[m]) != 0.0 ? 1.0 / W : 1.0));
			DS += W * D[i].E;
		}
		D[n].E = DS / WSum;
		for (i=n+1; i<m; i++) D[i].E = 0.0;
	}
	Destroy(Data, ND);
	for (n=0; (n < ND ? D[n].E != 0.0 : false); n++) ;
	for (m=n+1; m < ND; m++) if (D[m].E != 0.0) D[n++] = D[m];
	ND = n;
	uIA[0] = IsoTA[Iso->refIso];
	uIB[0] = IsoTB[Iso->refIso];
	for (i=m=0; i < nIsoA; i++) if (IsoCA[i] > m && i != uIA[0])
	{
		m = IsoCA[i];
		uIA[1] = i;
	}
	if (m>0) 
	{
		for (i=m=0; i < nIsoA; i++) if (IsoCA[i] > m && i != uIA[0] && i != uIA[1])
		{
			m = IsoCA[i];
			uIA[2] = i;
		}
		nIA = (m>0 ? 3 : 2);
	}
	else nIA = 1;
	for (i=m=0; i < nIsoB; i++) if (IsoCB[i] > m && i != uIB[0])
	{
		m = IsoCB[i];
		uIB[1] = i;
	}
	if (m>0)
	{
		for (i=m=0; i < nIsoB; i++) if (IsoCB[i] > m && i != uIB[0] && i != uIB[1])
		{
			m = IsoCB[i];
			uIB[2] = i;
		}
		nIB = (m>0 ? 3 : 2);
	}
	else nIB = 1;
	for (n=0; n < nIA; n++) 
		IsoSA[n] = ("    " + QString::number(A1->getnNuc(uIA[n]))).right(5)
				 + ("          " 
				   + QString::number(A1->getIsoMass(uIA[n]), 'f', 7)).right(13)
				 + "           0.0";
	for (n=0; n < nIB; n++) 
		IsoSB[n] = ("    " + QString::number(A2->getnNuc(uIB[n]))).right(5)
				 + ("          " 
				   + QString::number(A2->getIsoMass(uIB[n]), 'f', 7)).right(13)
				 + "           0.0";
	QFile File(Filename);
	QString Buffer, L1, LP1, LP2, IsoStr[nIso];
	QuestionBox *QB;
	for (n=0; n < nIso; n++)
	{
		for (m=0; m < nIA; m++) if (Iso->mNumIso1[n] == IsoSA[m].left(5).toInt()) IsoStr[n] += IsoSA[m].left(5);
		for (m=0; m < nIB; m++) if (Iso->mNumIso2[n] == IsoSB[m].left(5).toInt()) IsoStr[n] += IsoSB[m].left(5);
	}
	if (File.exists())
	{
		QB = new QuestionBox("MolSpektAnalysis", 
								"The file '" + Filename 
							+ "' already exists, do you want do update or overwrite it?", 
							this);
		QB->addButton("Update");
		QB->addButton("Overwrite");
		QB->addButton("Cancel");
		QB->exec();
		int R = QB->getResult(), nIAF, nIBF = 0;
		delete QB;
		if (R == 0)
		{
			File.open(QIODevice::ReadOnly);
			QTextStream RS(&File);
			L1 = RS.readLine();
			NFD = L1.left(5).toInt();
			if ((n = L1.mid(20, 5).toInt()) != lambda && lambda >= 0)
			{
				QB = new QuestionBox("MolSpektAnalysis", 
						"The value for lambda in the file of " + QString::number(n)
						+ " is different to the value of the excited state of " 
						+ QString::number(lambda) 
						+ ". Do you want to keep or change it?", this);
				QB->addButton("Keep it");
				QB->addButton("Changed it");
				QB->addButton("Cancel");
				QB->exec();
				i = QB->getResult();
				delete QB;
				if (i==1)
					L1 = L1.left(20) + ("     " + QString::number(lambda)).right(5)
					     + L1.right(L1.length() - 25);
				else if (i==2)
				{
					delete Iso;
					delete[] D;
					return true;
				}
			}
			if (nIA != (nIAF = L1.mid(10, 5).toInt()) 
				|| nIB != (nIBF = L1.mid(15, 5).toInt())) Diff = true;
			else Diff = false;
			for (n=0; n < nIAF; n++)
			{	
				IsoSAF[n] = RS.readLine();
				if (!Diff) if (IsoSAF[n] != IsoSA[n]) Diff = true;
			}
			for (n=0; n < nIBF; n++)
			{
				IsoSBF[n] = RS.readLine();
				if (!Diff) if (IsoSBF[n] != IsoSB[n]) Diff = true;
			}
			if (Diff)
			{
				QB = new QuestionBox("MolSpektAnalysis", 
						"The isotope data in the file is different than expected.", this);
				QB->addButton("Keep it");
				QB->addButton("Changed it");
				QB->addButton("Cancel");
				QB->exec();
				i = QB->getResult();
				delete QB;
				switch (i) 
				{
					case 0:
						for (n=0, nIA = nIAF; n < nIA; n++) IsoSA[n] = IsoSAF[n];
						for (n=0, nIB = nIBF; n < nIB; n++) IsoSB[n] = IsoSBF[n];
						break;
					case 1:
						L1 = L1.left(10) + ("     " + QString::number(nIA)).right(5)
						   + ("     " + QString::number(nIB)).right(5)
						   + L1.right(L1.length() - 20);
						break;
					case 2:
						delete Iso;
						delete[] D;
						return true;
						break;
				}
			}
			LP1 = RS.readLine();
			LP2 = RS.readLine();
			FD = new TermEnergy[NFD];
			for (n=m=0; n < NFD; n++)
			{
				Buffer = RS.readLine();
				iA = Buffer.left(5).toInt();
				iB = Buffer.mid(5, 5).toInt();
				FD[m].v = Buffer.mid(10, 5).toInt();
				FD[m].J = Buffer.mid(15, 5).toInt();
				if (FD[m].J != Buffer.mid(25, 5).toInt()) FD[m].ef = true;
				else FD[m].ef = false;
				Buffer = RS.readLine();
				for (i = 0; (i < nIso ? iA != Iso->mNumIso1[i] || iB != Iso->mNumIso2[i] 
									  : false); i++) ;
				if (i < nIso)
				{
					FD[m].E = Buffer.left(15).toDouble();
					FD[m].err = Buffer.mid(15, 15).toDouble();
					FD[m].Iso = i;
					for (p1 = m, p2 = (m - 1) / 2; p1 > 0 && isnSPG(FD[p1], FD[p2]);
							p1 = p2, p2 = (p1 - 1) / 2)
					{
						TBuff = FD[p1];
						FD[p1] = FD[p2];
						FD[p2] = TBuff;
					}
					m++;
				}
			}
			File.close();
		}
		else if (R==2)
		{
			delete[] D;
			delete Iso;
			return true;
		}
	}
	if (L1.isEmpty())
	{
		L1 = ("     " + QString::number(ND)).right(5) + "    0"
		   + ("     " + QString::number(nIA)).right(5)
		   + ("     " + QString::number(nIB)).right(5)
		   + ("     " + QString::number(lambda)).right(5)
		   + "          delta v=10 und delta energy 0.0016";
		LP1 = "    0" + (St != 0 ? (St->getDunTable() != 0 ? ("         " 
				          + QString::number(St->getDunTable()->getwe(), 'f', 2)).right(9) : "    85.90") : "    85.90")
			+ "   -0.00       0" + ("     " + QString::number(Maxv)).right(5) + "    0"
			+ ("     " + QString::number(MaxJ)).right(5) + "       v Verschiebung, vibr. spacing, level limits";
		LP2 = " 0.0000            0               energy shift of levels to minimum, shift of quantum number position";
	}
	QString *LL = new QString[ND + NFD];
	for (n=m=0; n < ND || (NFD > 0 ? FD[0].Iso != -1 : false); m++)
	{
		if (n < ND && NFD > 0 ? D[n].Iso == FD[0].Iso && D[n].v == FD[0].v && D[n].J == FD[0].J && D[n].ef == FD[0].ef
			  : false)
		{
			TBuff = D[n++];
			TBuff.err = FD[0].err;
			FD[0].Iso = -2;
		}
		else if (n < ND ? (NFD > 0 && FD[0].Iso != -1 ? isnSPG(D[n], FD[0]) : true) : false) TBuff = D[n++];
		else
		{
			switch (use)
			{
				case -1:
					TBuff.Iso = -1;
					break;
				case 0:
					QB = new QuestionBox("MolSpektAnalysis", 
						(FD[0].ef ? "The e-level of the isotope=" : "The f-level of the isotope ") 
						 + QString::number(FD[0].Iso) + " with v'=" + QString::number(FD[0].v)
						 + " and J'=" + QString::number(FD[0].J) 
						 + " is inside the file but not inside the current data set. Nevertheless keep it?", this);
					QB->addButton("Yes");
					QB->addButton("Yes to all");
					QB->addButton("No");
					QB->addButton("No to all");
					QB->addButton("Cancel");
					QB->exec();
					i = QB->getResult();
					delete QB;
					switch (i)
					{
						case 1:
							use = 1;
                            /* Falls through. */
						case 0:
							TBuff = FD[0];
							break;
						case 3:
							use = -1;
                            /* Falls through. */
						case 2:
							TBuff.Iso = -1;
							break;
						case 4:
							delete[] LL;
							delete[] FD;
							delete[] D;
							delete Iso;
							return true;
							break;
					}
					break;
				case 1:
					TBuff = FD[0];
					break;
			}
			FD[0].Iso = -2;
		}
		if (NFD > 0 ? FD[0].Iso == -2 : false) 
		{
			for (p1 = 0, p2 = 2; p2 <= NFD && FD[p1].Iso != -1; p2 = 2 * ((p1 = p2) + 1))
			{
				if (p2 < NFD ? isnSPG(FD[p2-1], FD[p2]) : true) p2--;
				FD[p1] = FD[p2];
			}
			FD[p1].Iso = -1;
		}
		if (TBuff.Iso == -1 ? true : IsoStr[TBuff.Iso].length() < 10) m--;
		else 
		{
			LL[m] = IsoStr[TBuff.Iso] + ("     " + QString::number(TBuff.v)).right(5) 
			      + ("     " + QString::number(TBuff.J)).right(5) + "   -1"
			      + ("     " + QString::number(TBuff.ef ? TBuff.J + 1 : TBuff.J)).right(5)
			      + IsoStr[TBuff.Iso]
			      + "    0    0    0    0\n" + ("               " + QString::number(TBuff.E, 'f', 4)).right(15)
			      + ("               " + QString::number(TBuff.err, 'f', 3)).right(15) + "    2\n";
		}
	}
	if (L1.left(5).toInt() != m) L1 = ("     " + QString::number(m)).right(5) + L1.right(L1.length() - 5);
	File.open(QIODevice::WriteOnly);
	QTextStream WS(&File);
	WS << L1 << "\n";
	for (n=0; n < nIA; n++) WS << IsoSA[n] << "\n";
	for (n=0; n < nIB; n++) WS << IsoSB[n] << "\n";
	WS << LP1 << "\n" << LP2 << "\n";
	for (n=0; n<m; n++) WS << LL[n];
	File.close();
	delete[] LL;
	delete[] FD;
	delete[] D;
	delete[] OUnc;
	delete Iso;
	return true;
}

void LineTable::WriteTFGS(QString FileName, int vsO, vsOListElement *vsOList)
{
	if (Iso == 0)
	{
		QMessageBox::information(this, "MolSpektAnalysis", 
			"Error: The line table has to be assigned to a molecule first!", QMessageBox::Ok);
		return;
	}
	int i, j, N, n, vs, lvs = 0, P, lP = -1, Js, lJs = -1, Jss;
	int ivsO = vsO, PN, lPN = -1, aIso, lIso = -1;
	N = Tab->rowCount();
	if (FileName.isEmpty()) 
		FileName = QFileDialog::getSaveFileName(this, "Write TF input file", 
			(MW != 0 ? MW->getDir(getType()) : ""), "All files (*.*)");
	if (FileName.isEmpty()) return;
	QFile Datei(FileName);
	Datei.open(QIODevice::WriteOnly | QIODevice::Append);
	QTextStream S(&Datei);
	QString Buffer, IBuff, F, lF;
	int *SO = new int[N], *S1 = heapSort(sortIJvP);
	vsOListElement *CvsOElement = vsOList;
	if (CvsOElement == 0)
	{
		CvsOElement = vsOList = new vsOListElement;
		CvsOElement->Iso = CvsOElement->Js = CvsOElement->vs = -1;
		CvsOElement->next = 0;
	}
	for (n=0; n<N; n++) SO[S1[n]] = n;
	//printf("Nach Sortprog\n");
	for (i=0; i<N; i++) 
	{
		n = Tab->item(SO[i], CIso)->text().toInt();
		if (10 * (j = n / 10) != n) continue;
		if (--j < 0) continue;
		S << (IBuff = ("    " + QString::number(Iso->mNumIso1[j])).right(5) 
				+ ("    " + QString::number(Iso->mNumIso2[j])).right(5));
		//printf("Nach1\n");
		vs = Tab->item(SO[i], Cvs)->text().toInt();
		Js = Tab->item(SO[i], CJs)->text().toInt();
		Jss = Tab->item(SO[i], CJss)->text().toInt();
		P = fabs(Js - Jss);
		F = Tab->item(SO[i], CFile)->text();
		PN = Tab->item(SO[i], CPN)->text().toInt();
		aIso = Tab->item(SO[i], CIso)->text().toInt();
		if (F.left(5) == "laser") F = F.right(F.length() - 6);
		if (lvs != vs || lJs != Js || lIso != aIso) 
		{
			while (CvsOElement->next != 0 ? CvsOElement->next->Iso <= aIso && CvsOElement->next->Js <= Js 
					&& CvsOElement->vs <= vs : false) CvsOElement = CvsOElement->next;
			if (CvsOElement->Iso == aIso && CvsOElement->Js == Js && CvsOElement->vs == vs)
				CvsOElement->curMaxOffset = ivsO = (CvsOElement->curMaxOffset >= vsO ? CvsOElement->curMaxOffset + 100 : vsO + vs);
			else ivsO = vsO + vs;
			lvs = vs;
			lJs = Js;
			lF = F;
			lP = P;
			lPN = PN;
			lIso = aIso;
		}
		else if (F != lF || P != lP || PN != lPN) 
		{
			if (CvsOElement->Iso == aIso && CvsOElement->Js == Js && CvsOElement->vs == vs)
				CvsOElement->curMaxOffset = ivsO = CvsOElement->curMaxOffset + 100;
			else if ((ivsO += 100) >= vsO + 1000)
			{
				vsOListElement *vsOBuff = new vsOListElement;
				vsOBuff->Iso = aIso;
				vsOBuff->vs = vs;
				vsOBuff->Js = Js;
				vsOBuff->curMaxOffset = ivsO;
				vsOBuff->next = CvsOElement->next;
				CvsOElement = CvsOElement->next = vsOBuff;
			}
			lF = F;
			lP = P;
			lPN = PN;
			lIso = aIso;
		}
		Buffer = "    " + QString::number(ivsO);
		S << Buffer.right(5);
		Buffer = "    " + QString::number(Js);
		S << Buffer.right(5);
		S << ("     " + Tab->item(SO[i], Cvss)->text()).right(5);
		S << ("     " + Tab->item(SO[i], CJss)->text()).right(5);
		S << IBuff << "    0    0    0    0\n";
		Buffer = "               ";
		Buffer += Tab->item(SO[i], CWN)->text();
		if (Buffer.indexOf(".") == -1) Buffer += ".00";
		S << Buffer.right(15);
		S << ("               " + Tab->item(SO[i], Cerr)->text()).right(15) << "    2\n";
	}
	Datei.close();
	delete[] S1;
	delete[] SO;
	if (transition != 0 && molecule != 0 && vsO == 1000)
	{
		ElState *S = transition->getLowerState();
		Transition *T;
		LineTable *L;
		for (i=0, N = molecule->getNumTransitions(); i<N; i++) 
			if ((T = molecule->getTransitionP(i))->getLowerState() == S) 
				//for (j=0, n = T->getNumLineTables(); j<n; j++)
		{
			L = T->getLineTable();
			if (L != this && L != 0) L->WriteTFGS(FileName, vsO += 1000, vsOList);
		}
	}
	while (vsOList != 0)
	{
		CvsOElement = vsOList->next;
		delete vsOList;
		vsOList = CvsOElement;
	}
}

void LineTable::W2AI(QString sldc)
{	
	if (MW == 0) return;
	if (molecule == 0 || transition == 0)
	{
		QMessageBox::information(this, "MolSpektAnalysis", "Error: This line table has to be assigned to a molecule and an electronic transition first!");
		return;
	}
	if (!ShowUpTerm()) return;
	IsoTab *IsoT = molecule->getIso();
	int k, l, i, j, N, n, mvsldc = -1, v, JF, J, ZINC = 0, Mv = -1, mv = 1000, MJ = 0, mJ = 1000;
    double **sLD = 0;
	ElState *US = transition->getUpperState();
	k = (US != 0 ? US->getOmega() : 0);
	int OmegaSQ = k*k, MinDPR = 0, MinQL = 0, ZINI1 = 0, ZINI2 = 0;
	QString Buffer, IBuff, ZIL[4], *ZIC = 0, *ZII1 = 0, *ZII2 = 0, *Iso;
	QString ZInpFile = QFileDialog::getSaveFileName(this, "Write ZweiAt-Input.dat", MW->getDir(), 
									"Zweiat Inputdateien(*.dat)", 0, QFileDialog::DontConfirmOverwrite);
	if (ZInpFile.isEmpty()) return;
	QFile Datei(ZInpFile);
	if (Datei.exists())
		if (QMessageBox::question(this, "MolSpektAnalysis", 
			"The selected file already exists, do you want to update the frequency list only?",
   			QMessageBox::Yes | QMessageBox::No, QMessageBox::Yes) == QMessageBox::Yes)
	{
		if (!Datei.open(QIODevice::ReadOnly)) 
    	{
			QString Fehlermeldung = "Die Datei " + ZInpFile + " konnte nicht geöffnet werden.";
			QMessageBox::information( this, "MolSpektAnalysis", Fehlermeldung, QMessageBox::Ok);
			return;
		}
		printf("Einlesen der Isotopendaten\n");
		QTextStream SI(&Datei);
		for (i=0; i<4; i++) ZIL[i] = SI.readLine();
		ZINC = ZIL[3].mid(5, 5).toInt();
		ZINI1 = ZIL[3].mid(10, 5).toInt();
		ZINI2 = ZIL[3].mid(15, 5).toInt();
		OmegaSQ = ZIL[3].mid(20, 5).toInt();
		ZIC = new QString[ZINC];
		for (i=0; i < ZINC; i++) ZIC[i] = SI.readLine();
		ZII1 = new QString[ZINI1];
		ZII2 = new QString[ZINI2];
		for (i=0; i < ZINI1; i++) ZII1[i] = SI.readLine();
		for (i=0; i < ZINI2; i++) ZII2[i] = SI.readLine();
		Datei.close();
	}
	if (!Datei.open(QIODevice::WriteOnly)) 
    {
		QString Fehlermeldung = "Die Datei " + ZInpFile + " konnte nicht geöffnet werden.";
		QMessageBox::information( this, "MolSpektAnalysis", Fehlermeldung, QMessageBox::Ok);
		return;
    }
	//printf("Vor ShowUpTermTable\n");
	ShowUpTermTable();
	//printf("Nach ShowUpTermTable\n");
	if (termTable == 0) return;
	QString **Data = termTable->getData(N, i);
	if (!sldc.isEmpty())
	{
		QFile iDatei(sldc);
		if (!iDatei.open(QIODevice::ReadOnly))
		{
			QMessageBox::information( this, "MolSpektAnalysis", 
				"Fehler beim Lesen der Datei " + sldc + "!", QMessageBox::Ok);
			return;
		}
		QTextStream SR(&iDatei);
		Buffer = SR.readLine();
		if (Buffer.left(4) != "sldc") 
		{
			QMessageBox::information( this, "MolSpektAnalysis", 
				"Die Datei " + sldc + "hat nicht das richtige Dateiformat!", QMessageBox::Ok);
			return;
		}
		QStringList SL;
		for (i=0; i<N; i++) if ((j = Data[i][1].toInt()) > mvsldc) mvsldc = j; 
		sLD = Create(mvsldc + 1, 2);
		for (i=0; i <= mvsldc; i++) sLD[i][0] = sLD[i][1] = 0.0;
		while (!SR.atEnd())
		{
			SL = SR.readLine().split("	");
			if ((n = SL.count()) > 1) i = SL[0].toInt();
			else i=-1;
			if (i < 0 || i > mvsldc) continue;
			if (n >= 2) sLD[i][0] = SL[1].toDouble();
			if (n >= 4) sLD[i][1] = SL[3].toDouble();
		}
	}
	//printf("Nach Einlesen von sldc\n");
    j = N;
	for (n=0, N=0; (N<j ? Data[N][0] != "" : false); N++)
		if ((Data[N][6].toInt() >= MinQL && Data[N][3] == "f")
			|| (Data[N][7].toInt() >= MinDPR 
				   && Data[N][3] == "e"))
	{
		l = (k = Data[N][0].toInt()) / 10;
		if (l > 0 && k == 10 * l) n++;
		i = Data[N][1].toInt();
		if (i < mv) mv = i;
		if (i > Mv) Mv = i;
		i = Data[N][2].toInt();
		if (i < mJ) mJ = i;
		if (i > MJ) MJ = i;
	}
	//printf("Nach TermTable\n");
	if (mv < 0) mv = 0;
    if (ZINI1 == 0 || ZINI2 == 0)
	{
		Atom *atom = molecule->getAtom1();
		int RI1, RI2;
		molecule->getRefIso(RI1, RI2);
		printf("RI1=%d, RI2=%d\n", RI1, RI2);
		ZINI1 = atom->getnIso();
		ZII1 = new QString[ZINI1];
		for (i=0; i < ZINI1; i++)
			ZII1[i] = ("    " + QString::number(atom->getnNuc(i))).right(5) 
					+ "   " + QString::number(atom->getIsoMass(i), 'f', 8);
		for (i = RI1, Buffer = ZII1[RI1]; i>0; i--) ZII1[i] = ZII1[i-1];
		ZII1[0] = Buffer;
		atom = molecule->getAtom2();
		ZINI2 = atom->getnIso();
		ZII2 = new QString[ZINI2];
		for (i=0; i < ZINI2; i++)
			ZII2[i] = ("    " + QString::number(atom->getnNuc(i))).right(5)
					+ "   " + QString::number(atom->getIsoMass(i), 'f', 8);
		for (i = RI2, Buffer = ZII2[RI2]; i>0; i--) ZII2[i] = ZII2[i-1];
		ZII2[0] = Buffer;
	}
	//printf("Nach Iso, ZINC=%d, US=%d\n", ZINC, US);
	if (ZINC <= 0)
	{
		DunTable *DT = (US != 0 ? US->getDunTable() : 0);
		if (DT != 0)
		{
			int KMax, *LMax;
			double **Par, **Kor, **SpinR, **adCorr;
			DT->getData(KMax, LMax, Par, Kor, SpinR, adCorr);
			for (i = ZINC = 0; i <= KMax; i++) 
			{
				ZINC += LMax[i] + 1;
				if (Kor != 0) for (j=0; j <= LMax[i]; j++) if (Kor[i][j] != 0.0) ZINC++;
			}
			ZIC = new QString[ZINC];
			for (i=j=0; j <= KMax; j++) for (k=0; k <= LMax[j]; k++)
			{
				ZIC[i++] = ("    " + QString::number(k)).right(5) 
						 + ("    " + QString::number(j)).right(5)
						 + (Par[j][k] != 0.0 ? "    1    10.0000000000000000E+00    0    2"
										     : "    1    10.0000000000000000E+00    1    2");
				if (Kor != 0 ? Kor[j][k] != 0.0 : false) 
					ZIC[i++] = ("    " + QString::number(k)).right(5) 
							 + ("    " + QString::number(j)).right(5)
							 + "    6    10.0000000000000000E+00    0    2";
			}
			Destroy(Par, KMax);
			if (Kor != 0) Destroy(Kor, KMax);
			if (SpinR != 0) Destroy(SpinR, KMax);
			if (adCorr != 0) Destroy(adCorr, KMax);
			delete[] LMax;
		}
		//printf("Vor ZIC\n");
		if (ZINC <= 0)
		{
			ZINC = 10;
			ZIC = new QString[10];
			ZIC[0] = "    0    0    1    10.0000000000000000E+00    0    2";
			ZIC[1] = "    1    0    1    10.0000000000000000E+00    0    2";
			ZIC[2] = "    2    0    1    10.0000000000000000E+00    0    2";
			ZIC[3] = "    3    0    1    10.0000000000000000E+00    0    2";
			ZIC[4] = "    0    1    1    10.0000000000000000E+00    0    2";
			ZIC[5] = "    1    1    1    10.0000000000000000E+00    0    2";
			ZIC[6] = "    2    1    1    10.0000000000000000E+00    0    2";
			ZIC[7] = "    0    2    1    10.0000000000000000E+00    0    2";
			ZIC[8] = "    1    2    1    10.0000000000000000E+00    0    2";
			ZIC[9] = "    0    3    1    10.0000000000000000E+00    0    2";
		}
		//printf("Vor ZIL\n");
		ZIL[0] = "    0    0    2    2";
		ZIL[1] = ("    " + QString::number(mv)).right(5) + ("    " + QString::number(Mv)).right(5)
			   + ("    " + QString::number(mJ)).right(5) + ("    " + QString::number(MJ)).right(5)
			   + "   -1    0" + ("    " + QString::number(mJ + 1)).right(5) 
			   + ("    " + QString::number(MJ + 1)).right(5) + "   0.003  2";
		ZIL[2] = molecule->getName() + ", term energies of " 
			   + (US != 0 ? "the state " + US->getName() : QString("an unknown state"))
			   + ", calculated from lines assigned with the program MolSpektAnalysis";
		ZIL[3] = "     " + ("    " + QString::number(ZINC)).right(5) 
			   + ("    " + QString::number(ZINI1)).right(5) + ("    " + QString::number(ZINI2)).right(5)
			   + ("    " + QString::number(OmegaSQ)).right(5);
	}
	//printf("Vor IsoT\n");
	Iso = new QString[IsoT->numIso];
	for (i=0; i < IsoT->numIso; i++) 
		Iso[i] = ("    " + QString::number(IsoT->mNumIso1[i])).right(5) 
			   + ("    " + QString::number(IsoT->mNumIso2[i])).right(5);
	QTextStream S(&Datei);
	//printf("Vor write\n");
	ZIL[3] = ("     " + QString::number(n)).right(5) + ZIL[3].right(ZIL[3].length() - 5);
	for (i=0; i<4; i++) S << ZIL[i] << "\n";
	for (i=0; i < ZINC; i++) S << ZIC[i] << "\n";
	for (i=0; i < ZINI1; i++) S << ZII1[i] << "\n";
	for (i=0; i < ZINI2; i++) S << ZII2[i] << "\n";
	//printf("Vor Schleife\n");
	for (i=0; i<N; i++) 
		if ((Data[i][6].toInt() >= MinQL && Data[i][3] == "f")
		   || (Data[i][7].toInt() >= MinDPR && Data[i][3] == "e"))
    {
		l = (k = Data[i][0].toInt()) / 10;
		if (l <= 0 || k != 10 * l) continue;
		S << (IBuff = Iso[l-1]);
		S << ("     " + Data[i][1]).right(5) 
				<< ("    " + Data[i][2]).right(5) << "   -1";
		if (Data[i][3] == "f") 
			S << ("     " + Data[i][2]).right(5);
		else 
		{
		    Buffer = "     ";
	    	Buffer += QString::number(Data[i][2].toInt() + 1);
	    	S << Buffer.right(5);
		}
		S << IBuff << "    0    0    0    0\n";
		Buffer = "               ";
		if (mvsldc == -1 || Data[i][3] == "f") 
			Buffer += Data[i][4];
		else 
		{
			v = Data[i][1].toInt();
			J = Data[i][2].toInt();
			JF = J * (J + 1) - OmegaSQ;
			Buffer += QString::number(Data[i][4].toDouble() 
					+ JF * (sLD[v][0] + JF * sLD[v][1]), 'g', 11);
		}
		S << Buffer.right(15);
		Buffer = "               ";
		/*IBuff = TermTable->text(i, 5);
		if (IBuff == "" || IBuff.toDouble() == 0.0) Buffer += "0.02";
		else Buffer += TermTable->text(i, 5);*/
		Buffer += "0.010";
		S << Buffer.right(15) << "    2\n";
    }
	S << "0\n";
	Destroy(Data, N);
	delete IsoT;
	delete[] Iso;
	delete[] ZIC;
	delete[] ZII1;
	delete[] ZII2;
}

void LineTable::Assignvs()
{
	//printf("LineTable::Assignvs()\n");
	int i, nr = Tab->rowCount();
	ShowUpTerm();
	//printf("Nach ShowUpTerm\n");
	double TE = 0.0, dB;
	int vs = 0, I=0, J=0;
	if (molecule == 0 || transition == 0)
	{
		QMessageBox::information(this, tr("QT4MolSpektAn"), 
						tr("The line table has to be assigned to a molecule and a transition first!"));
		return;
	}
	ElState *US = transition->getUpperState();
	if (US == 0)
	{
		QMessageBox::information(this, tr("QT4MolSpektAn"), 
						tr("The transition has to be assigned to an upper electronic state first!"));
		return;
	}
	TermTable *TT = US->getTermTable();
	//printf("Nach getTermTable\n");
	if (TT == 0) 
	{
		DunTable *DT = US->getDunTable();
		if (DT != 0) DT->calcTermEnergies(TT);
	}
	if (TT == 0)
	{
		QMessageBox::information(this, tr("QT4MolSpektAn"), 
			tr("A term energy table or a Dunham coefficient set for the upper state has to be loaded first!"));
		return;
	}
	QDialog *D = new QDialog(this);
	D->setWindowTitle("Select tolerancy");
	QGridLayout *L = new QGridLayout(D);
	L->addWidget(new QLabel("Please select the tolerancy\n for the v' assignment:", D), 0, 0, 1, 2);
	L->addWidget(new QLabel("Tolerancy [cm^-1]:", D), 1, 0);
	QLineEdit *Tol = new QLineEdit("5", D);
	L->addWidget(Tol, 1, 1);
	L->setRowMinimumHeight(2, 20);
	QPushButton *BOK = new QPushButton("OK", D);
	L->addWidget(BOK, 3, 0);
	QPushButton *BCancel = new QPushButton("Cancel", D);
	L->addWidget(BCancel, 3, 1);
	connect(BOK, SIGNAL(clicked()), D, SLOT(accept()));
	connect(BCancel, SIGNAL(clicked()), D, SLOT(reject()));
	if (D->exec() == QDialog::Rejected) return;
	double ***UE = 0, ****UD = TT->getData(), AT = Tol->text().toDouble();
	int NumC = TT->getNumComp(), mvs = TT->getMaxv();
	//printf("mvs=%d, NumC=%d\n", mvs, NumC);
	Tab->blockSignals(true);
	for (i=0; i<nr; i++)
	{
		//printf("i=%d, nr=%d, CEav=%d, CCalc=%d, N=%d\n", i, nr, CEav, CCalc, Tab->columnCount());
		dB = Tab->item(SO[i], CEUp)->text().toDouble();
		//printf("Z-\n");
		if (TE != dB)
		{
			//printf("Z0\n");
			TE = dB;
			//printf("Z1\n");
			if ((J = Tab->item(SO[i], CJs)->text().toInt()) 
				!= Tab->item(SO[i], CJss)->text().toInt() 
						  || NumC == 1) 
				UE = UD[0];
			else UE = UD[1];
			//printf("Z2\n");
			I = (Tab->item(SO[i], CIso)->text().toInt() - 1) / 10;
			//printf("vs=%d, mvs=%d, I=%d, J=%d\n", vs, mvs, I, J);
			for (vs=0; (vs <= mvs ? UE[I][vs][J] < TE : false); vs++) ;
			//printf("Z4\n");
			if ((vs > 0 ? (vs <= mvs ? UE[I][vs][J] - TE > TE - UE[I][vs-1][J] : true)
				: false)) vs--;
			//printf("Z5\n");
			if (fabs(UE[I][vs][J] - TE) > AT) vs = -1;
		}
		//printf("Z6\n");
		Tab->item(SO[i], Cvs)->setText(QString::number(vs));
		//printf("Z7\n");
		if (vs >= 0) 
		{
			//printf("Z8, I=%d, vs=%d, J=%d\n", I, vs, J);
			Tab->setItem(SO[i], CCalc, new QTableWidgetItem(QString::number(UE[I][vs][J], 'g', 9)));
			//printf("Z9\n");
			Tab->setItem(SO[i], COmC, new QTableWidgetItem(
				QString::number(TE - UE[I][vs][J], 'g', 6)));
		}
		else 
		{
			//printf("ZA\n");
			Tab->setItem(SO[i], CCalc, new QTableWidgetItem(""));
			//printf("ZB\n");
			Tab->setItem(SO[i], COmC, new QTableWidgetItem(""));
		}
	}
	Tab->blockSignals(false);
	Changed();
	//printf("Ende Assign v'\n");
}

void LineTable::AssignFC()
{
	if (transition == 0) 
	{
		QMessageBox::information(this, "QT4MolSpektAn", 
			"The line table has to be assigned to a molecule in order to be able to assign the lines to fine structure components!");
		return;
	}
	ElState *XState = transition->getLowerState(), *EState = transition->getUpperState();
	if (XState == 0 || EState == 0)
	{
		QMessageBox::information(this, "QT4MolSpektAn", 
			"The line table has to be assigned to a known transition in order to be able to assign the lines to fine structure components!");
		return;
	}
	TermTable *XTT = XState->getTermTable(), *ETT = EState->getTermTable();
	if (XTT == 0 || ETT == 0)
	{
		QMessageBox::information(this, "QT4MolSpektAn", 
			"For both electronic states have to be term energies available in order to be able to assign the lines to fine structure components!");
		return;
	}
	int XNC = XTT->getNumComp(), ENC = ETT->getNumComp(), *XCT, *ECT, n, m;
	if (XNC < 2 || ENC < 2)
	{
		QMessageBox::information(this, "QT4MolSpektAn", 
			"For both term energy tables have to be different fine structure components available in order to be able to assign the lines to fine structure components!");
		return;
	}
	XTT->getCompT(XNC, XCT);
	ETT->getCompT(ENC, ECT);
	int cTX[XNC], cTE[XNC], cTL[XNC], NI = molecule->getNumIso();
	for (n=m=0; n <= XNC && n <= ENC; n++) if (XCT[n] >= 0 && ECT[n] >= 0)
	{
		cTX[m] = XCT[n];
		cTE[m] = ECT[n];
		cTL[m++] = n;
	}
	if ((XNC = m) < 2)
	{
		QMessageBox::information(this, "QT4MolSpektAn", 
			"The available fine structure components for the two electronic states have not enough overlap or incompatible labeling!");
		return;
	}
	int XNv = XTT->getMaxv() + 1, ENv = ETT->getMaxv() + 1;
	int XNJ = XTT->getMaxJ() + 1, ENJ = ETT->getMaxJ() + 1, NCol = Tab->columnCount();
	double ****XData = XTT->getData(), ****EData = ETT->getData(), RS, BS = 0, F[1000];
	int *LSO = heapSort(sortIJvP), NR = Tab->rowCount(), PN, nPN = 0, i, j, vs, Js, I, c, bc = 0;
	int *XIT = XTT->getIsoT(), *EIT = ETT->getIsoT(), vss[1000], Jss[1000], LI[1000];
	int *LO = new int[NR];
	for (n=0; n < NR; n++) LO[LSO[n]] = n;
	delete[] LSO;
	Tab->blockSignals(true);
	for (n=0, m=1, PN = Tab->item(LO[0], CPN)->text().toInt(); m <= NR; m++)
		if (m < NR ? (nPN = Tab->item(LO[m], CPN)->text().toInt()) != PN : true)
	{
		I = (Tab->item(LO[n], CIso)->text().toInt() - 1) / 10;
		vs = Tab->item(LO[n], Cvs)->text().toInt();
		Js = Tab->item(LO[n], CJs)->text().toInt();
		if (NCol > COmC) for (i=n; i<m; i++) for (c = CEUp; c <= COmC; c++) 
					Tab->item(LO[i], c)->setText("");
		if ((I >= 0 && I < NI ? XIT[I] >= 0 && EIT[I] >= 0 : false)
			&& vs >= 0 && vs < ENv && Js >= 0 && Js < ENJ 
			&& EData[0][EIT[I]][vs][Js] != 0.0)
		{
			for (i=n, j=0; i<m; i++)
			{
				vss[j] = Tab->item(LO[i], Cvss)->text().toInt();
				Jss[j] = Tab->item(LO[i], CJss)->text().toInt();
				LI[j] = LO[i];
				if (vss[j] >= 0 && vss[j] < XNv && Jss[j] >= 0 && Jss[j] < XNJ
						&& XData[0][XIT[I]][vss[j]][Jss[j]] != 0.0) 
					F[j++] = Tab->item(LO[i], CWN)->text().toDouble();
			}
			if (j>0)
			{
				for (c=0; c < XNC; c++)
				{
					for (i=0, RS = 0.0; i<j; i++) 
						RS += fabs(XData[cTX[c]][XIT[I]][vss[i]][Jss[i]] + F[i]
								   - EData[cTE[c]][EIT[I]][vs][Js]);
					if (RS < BS || c==0)
					{
						bc = c;
						BS = RS;
					}
				}
				for (i=n; i<m; i++) 
					Tab->item(LO[i], CF)->setText(' ' + QString::number(cTL[bc]));
				if (NCol > COmC)
				{
					for (i=0, RS = 0.0; i<j; i++)
						RS += (F[i] += XData[cTX[c]][XIT[I]][vss[i]][Jss[i]]);
					RS /= j;
					for (i=0; i<j; i++)
					{
						Tab->item(LI[i], CEUp)->setText(QString::number(F[i], 'f', 4));
						Tab->item(LI[i], CEUma)->setText(
										QString::number(F[i] - RS, 'g', 8));
						Tab->item(LI[i], COmC)->setText(QString::number(F[i]
								- EData[cTE[c]][EIT[I]][vs][Js], 'g', 8));
					}
					for (i=n; i<m; i++)
					{
						Tab->item(LO[i], CEav)->setText(QString::number(RS, 'f', 4));
						Tab->item(LO[i], CEdJ)->setText(
							QString::number(RS / (Js * (Js + 1)), 'f', 6));
						Tab->item(LO[i], CCalc)->setText(QString::number(
							EData[cTE[c]][EIT[I]][vs][Js], 'f', 4));
					}
				}
			}
			else for (i=n; i<m; i++) Tab->item(LO[i], CF)->setText("-1");
		}
		else for (i=n; i<m; i++) Tab->item(LO[i], CF)->setText("-1");
		n=m;
		PN = nPN;
	}
	Tab->blockSignals(false);
	delete[] LO;
	delete[] XCT;
	delete[] XIT;
	delete[] EIT;
	delete[] ECT;
	Changed();
}

void LineTable::MarkLines(int* Iso, int* vs, int* Js, double* WN, int N)
{
	int n, lJs, lIso = -1, lvs = -2, r, r1, R = Tab->rowCount(), MC = Tab->columnCount() - 1;
	double lWN = 0.0;
	for (r=0; r<R; r++)
	{
		lJs = Tab->item(r, CJs)->text().toInt();
		if (lvs != -2)
		{
			lvs = -2;
			if (lIso != -1)
			{
				lIso = -1;
				lWN = 0.0;
			}
		}
		for (n=0; n<N; n++) if (lJs == Js[n]) 
			if ((lvs == -2 ? (lvs = Tab->item(r, Cvs)->text().toInt()) : lvs) == vs[n])
				if ((lIso == -1 ? (lIso = (Tab->item(r, CIso)->text().toInt() - 1) / 10) : lIso) == Iso[n])
					if (fabs((lWN == 0.0 ? (lWN = Tab->item(r, CWN)->text().toDouble()) : lWN) - WN[n]) < 1e-4) break;
		if (n<N) break;
	}
	if (r==R) return;
	QList<QTableWidgetSelectionRange> SL = Tab->selectedRanges();
	for (n=0; n < SL.count(); n++) Tab->setRangeSelected(SL[n], false);
	Tab->scrollToItem(Tab->item(r, 0), QAbstractItemView::PositionAtTop);
	for (r1 = r++; r<R; r++)
	{
		lJs = Tab->item(r, CJs)->text().toInt();
		if (lvs != -2)
		{
			lvs = -2;
			if (lIso != -1)
			{
				lIso = -1;
				lWN = 0.0;
			}
		}
		for (n=0; n<N; n++) if (lJs == Js[n]) 
			if ((lvs == -2 ? (lvs = Tab->item(r, Cvs)->text().toInt()) : lvs) == vs[n])
				if ((lIso == -1 ? (lIso = (Tab->item(r, CIso)->text().toInt() - 1) / 10) : lIso) == Iso[n])
					if (fabs((lWN == 0.0 ? (lWN = Tab->item(r, CWN)->text().toDouble()) : lWN) - WN[n]) < 1e-4) break;
		if (n<N && r1 == -1) r1 = r;
		else if (n==N && r1 != -1) 
		{
			Tab->setRangeSelected(QTableWidgetSelectionRange(r1, 0, r-1, MC), true);
			r1 = -1;
		}
	}
	if (!isVisible()) show();
	activateWindow();
	Tab->setFocus();
}

void LineTable::MarkSelected()
{
    //printf("LineTable::MarkSelected\n");
    if (MW == 0) return;
	int k=0, SR, r, AnzahlMarker, s, n;
    double WN;
	ElState *lState, *uState;
	if (transition != 0)
	{
		lState = transition->getLowerState();
		uState = transition->getUpperState();
	}
	else lState = uState = 0;
	QString lSN = (lState != 0 ? lState->getName() : ""), uSN = (uState != 0 ? uState->getName() : ""), Comment;
    QStringList SL;
	QList<int> iL;
	NR = Tab->rowCount();
	SelR = Tab->selectedRanges();
	for (s=0; s < SelR.count(); s++) for (r = SelR[s].topRow(); r <= SelR[s].bottomRow(); r++)
	{
		for (mSpectrum = Tab->item(r, CFile)->text(), n = r - 1; mSpectrum.isEmpty() 
				&& n>=0; n--) 
			mSpectrum = Tab->item(n, CFile)->text();
		if (!mSpectrum.isEmpty())
		{
			for (n=0; (n < SL.count() ? SL[n] != mSpectrum : false); n++) ;
			if (n == SL.count())
			{
				SL << mSpectrum;
				iL << 1;
			}
			else iL[n]++;
		}
	}
	for (n=r=0; n < iL.count(); n++) if (iL[n] > r) r = iL[s=n];
	mSpectrum = SL[s];
	mIso = Tab->item(SR = Tab->currentRow(), CIso)->text().toInt();
	mJs = Tab->item(SR, CJs)->text().toInt();
	mvs = Tab->item(SR, Cvs)->text().toInt();
	if (mSpectrum.isEmpty()) return;
	Spektrum *spektrum = MW->getSpectrum(mSpectrum);
	if (spektrum == 0) return;
	if (!spektrum->isVisible()) spektrum->show();
    Marker *marker;
	spektrum->GetMarker(AnzahlMarker, marker);
	if (AnzahlMarker == 0)
	{
		spektrum->editFind();
		spektrum->GetMarker(AnzahlMarker, marker);
		if (AnzahlMarker == 0) return;
	}
    for (s=0; s < SelR.count(); s++) for (r = SelR[s].topRow(); r <= SelR[s].bottomRow(); r++)
	{
		if (Tab->item(r, CFile)->text() != mSpectrum) continue;
		WN = Tab->item(r, CWN)->text().toDouble();
		while ((k >= 0 ? marker[k].Line[0] > WN : false)) k--;
		k++;
		while ((k < AnzahlMarker ? marker[k].Line[0] < WN || k==0 : false)) k++;
		if ((k < AnzahlMarker ? marker[k].Line[0] - WN > WN - marker[k-1].Line[0] : true)) k--;
		//printf("marker[%d].Line[0]=%f, WN=%f\n", k, marker[k].Line[0], WN);
		if (fabs(marker[k].Line[0] - WN) < 0.001)
		{
		    //printf("Markiere marker[%d]\n", k);
		    marker[k].vs = Tab->item(r, Cvs)->text().toInt();
		    marker[k].Js = Tab->item(r, CJs)->text().toInt();
		    marker[k].vss = Tab->item(r, Cvss)->text().toInt();
		    marker[k].Jss = Tab->item(r, CJss)->text().toInt();
			marker[k].Iso = (Tab->item(r, CIso)->text().toInt() - 1) / 10;
			marker[k].FC = Tab->item(r, CF)->text().toInt();
		    marker[k].DD = Tab->item(r, CDev)->text().toDouble();
			marker[k].Mol = molecule;
		    marker[k].IsoName = (Iso != 0 ? Iso->texName[marker[k].Iso] : "");
			marker[k].lState = lSN;
			marker[k].uState = uSN;
			marker[k].LState = lState;
			marker[k].UState = uState;
			marker[k].DisplayData = true;
		    marker[k].Marked = true;
			marker[k].uncertainty = Tab->item(r, Cerr)->text().toDouble();
			Comment = Tab->item(r, CC)->text();
			marker[k].satellite = Comment.indexOf("satellite", 0, Qt::CaseInsensitive) >= 0;
			marker[k].overlap = Comment.indexOf("overlap", 0, Qt::CaseInsensitive) >= 0;
		}
	}
	spektrum->FocusSelected();
    activateWindow();
    setFocus();
  	spektrum->activateWindow();
	spektrum->setFocus();
	//printf("Ende MarkSelected, mSpectrum=%s\n", mSpectrum.ascii());
}

void LineTable::ShowCalcRelInt(int NumWFPoints)
{
	Potential *lPot = 0, *uPot = 0;
	ElState *lS, *uS;
	if (transition != 0)
	{
		if ((lS = transition->getLowerState()) != 0) lPot = lS->getPotential();
		if ((uS = transition->getUpperState()) != 0) uPot = uS->getPotential();
	}
	if (lPot == 0 || uPot == 0)
	{
		QMessageBox::information(this, tr("QT4MolSpektAn"), 
			tr("To use this function for both assigned states potentials have to be available first!"));
		return;
	}
	SortProg();
	RemoveDoubled();
	SetPN();
	if (Tab->columnCount() <= CFCF) Tab->setColumnCount(CFCF + 1);
	Tab->setHorizontalHeaderItem(CFCF, new QTableWidgetItem("calc. rel. int."));
	int NR = Tab->rowCount(), r, n, s, PN, I, Js, vs, MJ, Mv, v; 
	TermTable *lTerm = lS->getTermTable(), *uTerm = uS->getTermTable();
	int uNv = uTerm->getMaxv() + 1, uNJ = uTerm->getMaxJ() + 1, uNIso = uTerm->getNumIso();
	int lNv = lTerm->getMaxv() + 1, lNJ = lTerm->getMaxJ() + 1, lNIso = lTerm->getNumIso();
	double ****uData = uTerm->getData(), ****lData = lTerm->getData();
	int lJ[2 * lNv];
	double lE[2 * lNv], FCF[2 * lNv], IS, FS;
	if (uNIso < lNIso) lNIso = uNIso;
	for (r=0; r < NR; r=s)
	{
		I = Tab->item(r, CIso)->text().toInt() / 10 - 1;
		//printf("r=%d, s=%d, I=%d\n", r, s, I);
		for (IS = 0.0, n = Mv = MJ = 0, PN = Tab->item(r, CPN)->text().toInt(), s=r; 
				   (s < NR ? Tab->item(s, CPN)->text().toInt() == PN : false); s++, n++)
		{
			IS += Tab->item(s, CSNR)->text().toDouble();
			if ((lJ[n] = Tab->item(s, CJss)->text().toInt()) > MJ) MJ = lJ[n];
			if ((v = Tab->item(s, Cvss)->text().toInt()) > Mv) Mv = v;
			if (MJ < uNJ && Mv < uNv && I < lNIso) lE[n] = lData[0][I][v][lJ[n]];
		}
	//printf("PN=%d, MJ=%d, lNJ=%d, Js=%d, uNJ=%d, I=%d, lNIso=%d\n", PN, MJ, lNJ, Js, uNJ, I, lNIso);
		if (MJ >= lNJ || Mv >= lNv  || (Js = Tab->item(r, CJs)->text().toInt()) >= uNJ 
				  || (vs = Tab->item(r, Cvs)->text().toInt()) >= uNv || I >= lNIso)
			continue;
		//printf("r=%d, s=%d, Js=%d, n=%d\n", r, s, Js, n);
        lPot->getFastFCF(I, uPot, Js, 0, uData[0][I][vs][Js], n, lE, lJ, FCF, NumWFPoints);
		//printf("Nach getFastFCF, n=%d\n", n);
		for (FS = 0.0, v=0; v<n; v++) 
		{
			//printf("v=%d, n=%d\n", v, n);
			//printf("FS=%g, FCF[%d]=%g\n", FS, v, FCF[v]);
			FS += FCF[v];
		}
		//printf("FS=%g\n", FS);
		for (IS /= FS, v=0, n=r; n<s; v++, n++) 
			Tab->setItem(n, CFCF, 
					new QTableWidgetItem(("     " + QString::number(FCF[v] * IS, 'g', 5)).right(6)));
	}
}

void LineTable::ShowGSDeviations()
{
	ElState *S;
	if (transition == 0)
	{
		QMessageBox::information(this, tr("MolSpektAnalysis"), 
			tr("To use this function the linetable as to be assigned to a transition first!"));
		return;
	}
	if ((S = transition->getLowerState()) == 0)
	{
		QMessageBox::information(this, tr("MolSpektAnalysis"), 
	tr("To use this function the linetable as to be assigned to an electronic ground state first!"));
		return;
	}
	TermTable *TT;
	if ((TT = S->getTermTable()) == 0)
	{
		QMessageBox::information(this, tr("MolSpektAnalysis"), 
			tr("To use this function for the ground state term energies have to be available first!"));
		return;
	}
	SetPN();
	int n, m, pn = 0, apn = -1, N = Tab->rowCount(), nv = TT->getMaxv() + 1, nI = TT->getNumIso();
	int nJ = TT->getMaxJ(), v, J, I, ldr = 0;
	double Dev, US = 0.0, SS = 0.0, Sig, *E = new double[N], *UT = new double[N], ****Data = TT->getData();
	double MD = 0.0, AD;
	bool l = false;
	int *S1 = heapSort(isnSPG), *SA = new int[N];
	for (n=0; n<N; n++) SA[S1[n]] = n;
	for (n=m=0; n<N; m++)
	{
		if (m<N ? (pn = Tab->item(SA[m], CPN)->text().toInt()) != apn : true)
		{
			for (US /= SS; n<m; n++)
			{
				if (UT[n] == 0.0) continue;
				Dev = UT[n] - US;
				Tab->item(SA[n], CDev)->setText(QString::number(Dev, 'f', 4));
				if (fabs(Dev) > 4.0 * E[n]) 
				{
					l = true;
					if ((AD = fabs(Dev) / E[n]) > MD) 
					{
						MD = AD;
						ldr = n;
					}
				}
			}
			//printf("n=%d, m=%d, SS=%g, US=%g\n", n, m, SS, US);
			if (n==N || (l && Tab->item(SA[ldr], CC)->text().indexOf("pd") == -1)) break;
			apn = pn;
			SS = US = 0.0;
		}
		if ((I = Tab->item(SA[m], CIso)->text().toInt() / 10 - 1) >= nI 
				   || (J = Tab->item(SA[m], CJss)->text().toInt()) >= nJ 
				   || (v = Tab->item(SA[m], Cvss)->text().toInt()) >= nv)
		{
			//printf("I=%d, nI=%d, J=%d, nJ=%d, v=%d, nv=%d\n", I, nI, J, nJ, v, nv);
			UT[m] = 0.0;
			continue;
		}
		if (Data[0][I][v][J] == 0.0)
		{
			//printf("Data[0][%d][%d][%d]==0!!\n", I, v, J);
			UT[m] = 0.0;
			continue;
		}
		E[m] = Tab->item(SA[m], Cerr)->text().toDouble();
		UT[m] = Data[0][I][v][J] + Tab->item(SA[m], CWN)->text().toDouble();
		Sig = 1.0 / (E[m] * E[m]);
		SS += Sig;
		US += Sig * UT[m];
	}
	if (l) Tab->selectRow(ldr);
	delete[] E;
	delete[] UT;
	delete[] S1;
	delete[] SA;
}

void LineTable::TestProgressions(int NumWFPoints)
{
    if (MW == 0) return;
	SortProg();
	RemoveDoubled();
	SetPN();
	int k=0, l=0, n, m=0, r, AnzahlMarker, PN = 1, s, RC, i, j;
    double WN, MDev, D;
	int NR = Tab->rowCount(), NC = Tab->columnCount();
	bool Success, *R = new bool[NR], sat;
	int *P = new int[NR];
	QString Spekt, Buffer;
	ElState *lState, *uState;
	Marker *marker, *laser;
	Spektrum *spektrum = MW->CreateSpectrum();
	if (transition != 0)
	{
		lState = transition->getLowerState();
		uState = transition->getUpperState();
	}
	else lState = uState = 0;
	Tab->blockSignals(true);
	for (r=0; r < NR; PN++)
	{
		for (s=r; (s < NR ? Tab->item(s, CPN)->text().toInt() == PN : false); s++) R[s-r] = true;
		printf("PN=%d, r=%d, s=%d, l=%d, NR=%d\n", PN, r, s, l, NR);
		for (Spekt = "", n=r; mSpectrum.isEmpty() && n<s; n++) Spekt = Tab->item(r, CFile)->text();
		if (!spektrum->readData(Spekt))
		{
			if (l>r) for (n=r; n<s; n++, l++) for (m=0; m < NC; m++) 
						Tab->item(l, m)->setText(Tab->item(n, m)->text());
			r=s;
			continue;
		}
		spektrum->editFind();
		spektrum->GetMarker(AnzahlMarker, marker);
		laser = spektrum->GetLaserLine();
		sat = (Tab->item(r, CC)->text().indexOf("satellite") > -1 ? true : false);
		printf("AnzahlMarker=%d\n", AnzahlMarker);
		for (Success = false, k=0; !Success; )
		{
			spektrum->ClearMarked();
			printf("r=%d, s=%d\n", r, s);
			for (n = r; n < s; n++) if (R[n-r])
			{
				WN = Tab->item(n, CWN)->text().toDouble();
				while ((k >= 0 ? marker[k].Line[0] > WN : false)) k--;
				k++;
				while ((k < AnzahlMarker ? marker[k].Line[0] < WN || k==0 : false)) k++;
				if ((k < AnzahlMarker ? marker[k].Line[0] - WN > WN - marker[k-1].Line[0] : true)) k--;
				//printf("marker[%d].Line[0]=%f, WN=%f\n", k, marker[k].Line[0], WN);
				if (fabs(marker[k].Line[0] - WN) < 0.001)
				{
		    		//printf("Markiere marker[%d]\n", k);
		    		marker[k].vs = Tab->item(n, Cvs)->text().toInt();
		    		marker[k].Js = Tab->item(n, CJs)->text().toInt();
		    		marker[k].vss = Tab->item(n, Cvss)->text().toInt();
		    		marker[k].Jss = Tab->item(n, CJss)->text().toInt();
					marker[k].Iso = Tab->item(n, CIso)->text().toInt() / 10 - 1;
					marker[k].LState = lState;
					marker[k].UState = uState;
					marker[k].DisplayData = true;
		    		marker[k].Marked = true;
					P[n] = k;
					printf("n=%d, r=%d, s=%d, k=%d\n", n, r, s, k);
				}
				else P[n] = -1;
			}
			printf("Nach Schleife\n");
            spektrum->TTransition(NumWFPoints);
			for (n=0, MDev = 0.0; n < AnzahlMarker; n++) if (marker[n].Marked) 
					if ((D = fabs(marker[n].DD)) > MDev)
			{
				MDev = D;
				m=n;
			}
			printf("m=%d, MDev=%e, MaxDev=%e\n", m, MDev, CMaxSearchDev);
			if (MDev > CMaxSearchDev) 
			{
				for (n=r; P[n] != m; n++) ;
				printf("n=%d, m=%d\n", n, m);
				R[n-r] = false;
			}
			else Success = true;
		}
		r=s;
		for (n=m=0; n < AnzahlMarker; n++) if (marker[n].Marked) 
		{
			m++;
			k=n;
		}
		printf("m=%d, l=%d, s=%d\n", m, l, s);
		if (m<4 || (uState != 0 && marker[k].vs == -1)) continue;
		if (l+m < s) 
		{
			for (k=0; k < AnzahlMarker; k++) if (marker[k].Marked)
			{
				Tab->item(l, CPN)->setText(QString::number(PN));
				Tab->item(l, Cvs)->setText(("    " + QString::number(marker[k].vs)).right(4));
				Tab->item(l, CJs)->setText(("    " + QString::number(marker[k].Js)).right(4));
				Tab->item(l, Cvss)->setText(("    " + QString::number(marker[k].vss)).right(4));
				Tab->item(l, CJss)->setText(("    " + QString::number(marker[k].Jss)).right(4));
				Tab->item(l, CIso)->setText(("    " 
						+ QString::number(10 * (marker[k].Iso + 1))).right(4));
				i = (marker[k].Line[0] < 10000 ? 1 : 0);
				if ((j = (Buffer = "      " 
								 + QString::number(marker[k].Line[0], 'g', 9 - i)).length() + i) < 16)
					Buffer += QString(".0000").right(16 - j);
				Tab->item(l, CWN)->setText(Buffer.right(15));
				Tab->item(l, Cerr)->setText(marker[k].overlap ? "     0.020" : "     0.005");
    			Tab->item(l, CFile)->setText(Spekt);
				Tab->item(l, CSNR)->setText(("    " 
						+ QString::number(marker[k].SNR, 'g', 5)).right(6));
				Tab->item(l, CDev)->setText(("    " 
						+ QString::number(marker[k].DD, 'g', 5)).right(10));
				if (sat) Buffer = "satellite";
				else if (marker + k == laser) Buffer = "laser";
				else Buffer = ""; 
				if (marker[k].overlap) Buffer += "overlap";
				Tab->setItem(l++, CC, new QTableWidgetItem(Buffer));
			}
			continue;
		}
		RC = Tab->rowCount();
		Tab->setRowCount(RC + m);
		for (n = RC, k=0; k < AnzahlMarker; k++) if (marker[k].Marked)
		{
			Tab->setItem(n, CPN, new QTableWidgetItem(QString::number(PN)));
			Tab->setItem(n, Cvs, new QTableWidgetItem(("    " 
					+ QString::number(marker[k].vs)).right(4)));
			Tab->setItem(n, CJs, new QTableWidgetItem(("    " 
					+ QString::number(marker[k].Js)).right(4)));
			Tab->setItem(n, Cvss, new QTableWidgetItem(("    " 
					+ QString::number(marker[k].vss)).right(4)));
			Tab->setItem(n, CJss, new QTableWidgetItem(("    " 
					+ QString::number(marker[k].Jss)).right(4)));
			Tab->setItem(n, CIso, 
				new QTableWidgetItem(("    " + QString::number(10 * (marker[k].Iso + 1))).right(4)));
			i = (marker[k].Line[0] < 10000 ? 1 : 0);
			if ((j = (Buffer = "      " 
						  + QString::number(marker[k].Line[0], 'g', 9 - i)).length() + i) < 16)
				Buffer += QString(".0000").right(16 - j);
			Tab->setItem(n, CWN, new QTableWidgetItem(Buffer.right(15)));
			Tab->setItem(n, Cerr, new QTableWidgetItem(
						 marker[k].overlap ? "     0.020" : "     0.005"));
    		Tab->setItem(n, CFile, new QTableWidgetItem(Spekt));
			Tab->setItem(n, CSNR, new QTableWidgetItem(
						 ("    " + QString::number(marker[k].SNR, 'g', 5)).right(6)));
			Tab->setItem(n, CDev, new QTableWidgetItem(
							 ("    " + QString::number(marker[k].DD, 'g', 5)).right(10)));
			if (sat) Buffer = "satellite";
			else if (marker + k == laser) Buffer = "laser";
			else Buffer = ""; 
			if (marker[k].overlap) Buffer += "overlap";
			Tab->setItem(n++, CC, new QTableWidgetItem(Buffer));
		}
	}
	RC = Tab->rowCount();
	for (n = NR; n < RC; n++, l++) for (k=0; k < NC; k++) Tab->setItem(l, k, Tab->takeItem(n, k));
	Tab->setRowCount(l);
	Tab->blockSignals(false);
	delete[] R;
	delete[] P;
	Changed();
}

void LineTable::TakeOnChanges()
{
	//printf("mSpectrum=%s\n", mSpectrum.ascii());
	Spektrum *Spekt = (MW != 0 ? MW->getSpectrum(mSpectrum) : 0);
	if (Spekt == 0) return;
	int j, k, n = Tab->rowCount(), c, nc = Tab->columnCount();
	QString B, Buffer;
	Tab->blockSignals(true);
	for (j=0; j < SelR.size(); j++) for (k=SelR[j].topRow(); k <= SelR[j].bottomRow(); k++)
		if (Tab->item(k, CFile)->text() == mSpectrum) for (c=0; c < nc; c++) Tab->setItem(k, c, 0);
	for (j=k=0; j <= lRow; j++) if (Tab->item(j, 0) == 0) k++;
	lRow -= k;
	for (j=0; (j < n ? Tab->item(j, 0) != 0 : false); j++) ;
    for (k=j; k < n; k++) if (Tab->item(k, 0) != 0)
	{
		for (c=0; c < nc; c++) Tab->setItem(j, c, Tab->takeItem(k, c));
		j++;
	}	
    Tab->setRowCount(j);
	Tab->blockSignals(false);
	writeData();
}

void LineTable::addData(QString **Data, int NR, int NC)
{
	int r, n, c, C = Tab->columnCount(), N = Tab->rowCount();
	Tab->blockSignals(true);
	if (NC < C) Tab->setColumnCount(NC);
	else if (C < NC) NC = C;
	Tab->setRowCount(N + NR);
	for (r=N-1; r >= N - NpL; r--) for (c=0; c < NC; c++) 
			Tab->setItem(r + NR, c, Tab->takeItem(r, c));
	for (n=0, r = N - NpL; n < NR; n++, r++) 
	{
		for (c=0; c < NC; c++) Tab->setItem(r, c, new QTableWidgetItem(Data[n][c]));
		delete[] Data[n];
	}
	delete[] Data;
	Tab->blockSignals(false);
	Changed();
}

void LineTable::AddMarked(int AnzahlMarker, Marker *marker, Marker *LaserLine, QString SpektFile)
{
	//printf("LineTable::AddMarked Anzahl Zeilen=%d\n", Tab->rowCount());
	AcceptAssignments(SpektFile, false);
    if (marker == NULL) return;
	int i, j, n, N, k=0, l, c, ColC = Tab->columnCount(), RC = Tab->rowCount();
	ElState *LState, *UState;
	if (transition == 0) 
	{
		LState = 0;
		UState = 0;
	}
	else
	{
		LState = transition->getLowerState();
		UState = transition->getUpperState();
	}
	//printf("Nach States\n");
	i = SpektFile.lastIndexOf(QRegExp("[\\/]")) + 1;
	j = SpektFile.indexOf('.', i);
	j = (j > 0 ? j : SpektFile.length()) - i;
	QString Buffer, SFN = SpektFile.mid(i, j);
	printf("%s\n", SFN.toLatin1().data());
    QPixmap Pix(10, 10);
    QPainter P(&Pix);
    P.setPen(QColor(255, 0, 0));
    P.setFont(QFont("Arial", 10));
    P.drawText(0, 10, "N");
	//printf("Vor Schleife\n");
	for (N=0, i=0; i < AnzahlMarker; i++) 
		if (marker[i].Marked && marker[i].LState == LState && marker[i].UState == UState) N++;
		//if (i == 19765 || i == 19766) printf("new lines found %d!\n", i);
	//}
	//printf("N=%d, LState=%d, UState=%d\n", N, LState, UState);
	Marker *Lines[N], *LBuff = marker;
	Tab->clearSelection();
	for (i=j=0; i < AnzahlMarker; i++) 
		if (marker[i].Marked && marker[i].LState == LState && marker[i].UState == UState)
	{ 
		for (n=0; (n < RC ? Tab->item(n, CFile)->text().indexOf(SFN) == -1 || marker[i].Iso != (Tab->item(n, CIso)->text().toInt() - 1) / 10 
				|| marker[i].FC != Tab->item(n, CF)->text().toInt() || marker[i].vs != Tab->item(n, Cvs)->text().toInt()
				|| marker[i].Js != Tab->item(n, CJs)->text().toInt() || marker[i].vss != Tab->item(n, Cvss)->text().toInt() 
				|| marker[i].Jss != Tab->item(n, CJss)->text().toInt() 
				|| fabs(marker[i].Line[0] - Tab->item(n, CWN)->text().toDouble()) > 1e-3 : false); n++) ;
		if (n < RC) Tab->setRangeSelected(QTableWidgetSelectionRange(k=n, 0, n, ColC - 1), true);
		else Lines[j++] = marker + i;
		marker[i].added = true;
	}
	if (j==0)
	{
		Tab->scrollToItem(Tab->item(k, 0));
		return;
	}
	NpL += (N=j);
	while (LBuff != 0) for (i=1, LBuff = 0; i<N; i++)
			if (Lines[i-1]->Iso > Lines[i]->Iso || (Lines[i-1]->Iso == Lines[i]->Iso 
						 && (Lines[i-1]->vs > Lines[i]->vs 
						 || (Lines[i-1]->vs == Lines[i]->vs && Lines[i-1]->Js > Lines[i]->Js))))
	{
		LBuff = Lines[i-1];
		Lines[i-1] = Lines[i];
		Lines[i] = LBuff;
	}
	Tab->blockSignals(true);
 	Tab->setRowCount((i = Tab->rowCount()) + N);
	for (n=0; n<N; n++)
	{
	    //printf("i=%d, n=%d, N=%d\n", i, n, N);
		if (n>0 ? Lines[n]->vs != Lines[n-1]->vs || Lines[n]->Iso != Lines[n-1]->Iso 
				  || Lines[n]->Js != Lines[n-1]->Js : true) NpProg++;
	    Tab->setItem(i, CPN, new QTableWidgetItem(Pix, QString::number(MaxPN + NpProg)));

		if (Lines[n]->DisplayData)
		{
			Tab->setItem(i, Cvs, new QTableWidgetItem(Pix, 
						 ("    " + QString::number(Lines[n]->vs)).right(4)));
			Tab->setItem(i, CJs, new QTableWidgetItem(Pix,
						 ("    " + QString::number(Lines[n]->Js)).right(4)));
			Tab->setItem(i, Cvss, new QTableWidgetItem(Pix,
						 ("    " + QString::number(Lines[n]->vss)).right(4)));
			Tab->setItem(i, CJss, new QTableWidgetItem(Pix,
						 ("    " + QString::number(Lines[n]->Jss)).right(4)));
			Tab->setItem(i, CIso, new QTableWidgetItem(Pix, 
						 ("    " + QString::number(10 * (Lines[n]->Iso + 1))).right(4)));
			Tab->setItem(i, CF, new QTableWidgetItem(Pix, 
						 ("    " + QString::number(Lines[n]->FC)).right(4)));
		}
		else 
		{
			Tab->setItem(i, Cvs, new QTableWidgetItem(Pix, "  -1"));
			Tab->setItem(i, CJs, new QTableWidgetItem(Pix, "  -1"));
			Tab->setItem(i, Cvss, new QTableWidgetItem(Pix, "  -1"));
			Tab->setItem(i, CJss, new QTableWidgetItem(Pix, "  -1"));
			Tab->setItem(i, CF, new QTableWidgetItem(Pix, "  -1"));
			Tab->setItem(i, CIso, new QTableWidgetItem(Pix, "  10"));
		}
		k = (Lines[n]->Line[0] < 10000 ? 1 : 0);
		if ((l = (Buffer = "      " + QString::number(Lines[n]->Line[0], 'g', 9 - k)).length() + k) 
				   < 16)
			Buffer += QString(".0000").right(16 - l);
		Tab->setItem(i, CWN, new QTableWidgetItem(Pix, Buffer.right(15)));
		Tab->setItem(i, Cerr, new QTableWidgetItem(Pix, 
					 (Lines[n]->DisplayData && Lines[n]->uncertainty > 0.0 ? 
						("          " + QString::number(Lines[n]->uncertainty, 'g', 6)).right(10) 
						: "     0.005")));
    	Tab->setItem(i, CFile, new QTableWidgetItem(Pix, SpektFile));
		Tab->setItem(i, CSNR, new QTableWidgetItem(Pix, 
							 ("    " + QString::number(Lines[n]->SNR, 'g', 5)).right(6)));
		Tab->setItem(i, CDev, new QTableWidgetItem(Pix, 
							 ("    " + QString::number(Lines[n]->DD, 'g', 5)).right(10)));
		if (Lines[n]->satellite) Buffer = "satellite";
		else if (Lines[n] == LaserLine) Buffer = "laser";
		else Buffer = ""; 
		if (Lines[n]->overlap) Buffer += (Buffer.isEmpty() ? "overlap" : ",overlap");
		for (c = CC + 1; c < ColC; c++) Tab->setItem(i, c, new QTableWidgetItem());
		Tab->setItem(i++, CC, new QTableWidgetItem(Pix, Buffer));
    }
	Tab->scrollToItem(Tab->item(i-1, 0));
	Tab->blockSignals(false);
	Changed();
	//printf("Ende Addmarked, Anzahl Zeilen=%d\n", Tab->rowCount());
}

void LineTable::AcceptAssignments(QString SpektFile, bool accept)
{
	//printf("LineTable::AcceptAssignments\n");
	QPixmap P;
	int NR = Tab->rowCount(), NC = Tab->columnCount(), r1, r2, r3, rd, c, n, m, PN, lPN;
	/*if (mSpectrum == SpektFile)
	{
		if (accept)
		{
			for (r1 = 0; r1 < NR; r1++)
			{
				if (Tab->item(r1, CFile)->text() != mSpectrum) continue;
				if (Tab->item(r1, CIso)->text().toInt() != mIso 
						  || Tab->item(r1, CJs)->text().toInt() != mJs
						  || Tab->item(r1, Cvs)->text().toInt() != mvs) continue;
				Tab->selectRow(r1);
			}
			DeleteRows();
			NR = Tab->rowCount();
		}
		else
		{
			mSpectrum = "";
			return;
		}
	}*/
	if (accept)
	{
		for (r1 = r2 = NR - NpL; (r2 < NR ? Tab->item(r2, CFile)->text() != SpektFile : false); r2++) ;
		for (r3 = r2; (r3 < NR ? Tab->item(r3, CFile)->text() == SpektFile : false); r3++) ;
		if (r2 > r1 && r3 > r2)
		{
			Matrix<QTableWidgetItem*> B(rd = r3 - r2, NC);
			for (n=0, m = r2; n < rd; n++, m++) for (c=0; c < NC; c++) 
					B.R[n][c] = Tab->takeItem(m, c);
			for (n = r3 - 1, m = r2 - 1; m >= r1; n--, m--) for (c=0; c < NC; c++)
					Tab->setItem(n, c, Tab->takeItem(m, c));
			for (n = r1, m=0; m < rd; n++, m++) for (c=0; c < NC; c++)
					Tab->setItem(n, c, B.R[m][c]);
			r3 = (r2 = r1) + rd;
		}
		if (r2 < r3) for (n=r2, PN = Tab->item(n, CPN)->text().toInt(), MaxPN++, NpProg--; n < r3; n++) 
		{
			if ((lPN = Tab->item(n, CPN)->text().toInt()) != PN) 
			{
				MaxPN++;
				NpProg--;
			}
			Tab->item(n, CPN)->setText(QString::number(MaxPN));
			for (c=0; c < NC; c++) Tab->item(n, c)->setIcon(P);
			//if ((r1 = Tab->item(n, 0)->text().toInt()) > MaxPN) MaxPN = r1;
		}
		NpL += r2 - r3;
		emit DataChanged();
	}
	else
	{
		//printf("NR = %d, NpL = %d\n", NR, NpL);
		for (r1 = NR - NpL; (r1 < NR ? Tab->item(r1, CFile)->text() != SpektFile : false); r1++) ;
		for (r2 = r1; (r2 < NR ? Tab->item(r2, CFile)->text() == SpektFile : false); r2++) ;
		for (n = r1, m = r2; m < NR; n++, m++) for (c=0; c < NC; c++)
				Tab->setItem(n, c, Tab->takeItem(m, c));
		Tab->setRowCount(n);
		NpL += r1 - r2;
	}
	//printf("Ende AcceptAssignments\n");
}

void LineTable::TabSelChanged()
{
	QList<QTableWidgetSelectionRange> SR = Tab->selectedRanges();
	int r, n, N, *Js, MJ, Mv, NIso, Jss, I, vss;
	bool UTA = true;
	double *E, ****ELU = 0, uE;
	ElState *LS;
	TermTable *T;
	if (Tab->columnCount() <= CEUp) UTA = false;
	else for (r=N=0; r < SR.count(); r++) for (n = SR[r].topRow(); n <= SR[r].bottomRow(); n++)
	{
		if (Tab->item(n, CEUp)->text().toDouble() > 0.0) N++;
		else UTA = false;
	}
	if (!UTA) if (transition != 0) if ((LS = transition->getLowerState()) != 0) if ((T = LS->getTermTable()) != 0)
	{
		ELU = T->getData();
		NIso = T->getNumIso();
		Mv = T->getMaxv();
		MJ = T->getMaxJ();
		if (ELU != 0)
		{
			for (r=N=0; r < SR.count(); r++) for (n = SR[r].topRow(); n <= SR[r].bottomRow(); n++)
				if ((I = (Tab->item(n, CIso)->text().toInt() - 1) / 10) < NIso 
						&& (vss = Tab->item(n, Cvss)->text().toInt()) <= Mv
						&& (Jss = Tab->item(n, CJss)->text().toInt()) < MJ) 
					if (vss >= 0 && Jss >= 0 && I >= 0 ? ELU[0][I][vss][Jss] != 0.0 : false) N++;
		}
		else if (N>0) UTA = true;
	}
	if (N>0 && ELU != 0)
	{
		Js = new int[N];
		E = new double[N];
		if (UTA)
		{
			for (r=N=0; r < SR.count(); r++) for (n = SR[r].topRow(); n <= SR[r].bottomRow(); n++)
				if ((uE = Tab->item(n, CEUp)->text().toDouble()) > 0.0)
			{
				Js[N] = Tab->item(n, CJs)->text().toInt();
				E[N++] = Tab->item(n, CEUp)->text().toDouble();
			}
		}
		else for (r=N=0; r < SR.count(); r++) for (n = SR[r].topRow(); n <= SR[r].bottomRow(); n++)
			if ((I = (Tab->item(n, CIso)->text().toInt() - 1) / 10) < NIso && (vss = Tab->item(n, Cvss)->text().toInt()) <= Mv
					&& (Jss = Tab->item(n, CJss)->text().toInt()) < MJ)
				if (vss >= 0 && Jss >= 0 && I >= 0 ? ELU[0][I][vss][Jss] != 0.0 : false)
		{
			Js[N] = Tab->item(n, CJs)->text().toInt();
			E[N++] = Tab->item(n, CWN)->text().toDouble() + ELU[0][I][vss][Jss];
		}
	}
	else
	{
		Js = 0;
		E = 0;
		N=0;
	}
	if (N > 0 || SelE != 0)
	{
		if (SelE != 0)
		{
			delete[] SelE;
			delete[] SelJs;
		}
		SelE = E;
		SelJs = Js;
		NSel = N;
		emit SelChanged();
	}
}

void LineTable::getSelData(int *&Js, double *&E, int &N)
{
	Js = SelJs;
	E = SelE;
	N = NSel;
}

void LineTable::getViewnE(int*& Js, double*& E, int& N)
{
	int NR = Tab->rowCount();
	TermTable *TT = (transition != 0 ? (transition->getLowerState() != 0 ? transition->getLowerState()->getTermTable() : 0) : 0);
	if (TT == 0)
	{
		N=0;
		Js = 0;
		E = 0;
		return;
	}
	double ****ELU = TT->getData();
	int n, i, J, v, NIso = TT->getNumIso(), Mv = TT->getMaxv(), MJ = TT->getMaxJ();
	bool *VR = new bool[NR];
	getViewnRows(VR);
	for (n=N=0; n < NR; n++) if (VR[n]) N++;
	Js = new int[N];
	E = new double[N];
	for (n=N=0; n < NR; n++) if (VR[n])
	{
		i = Tab->item(CIso, n)->text().toInt();
		J = Tab->item(CJss, n)->text().toInt();
		v = Tab->item(Cvss, n)->text().toInt();
		if (i>=0 && i < NIso && J>=0 && J <= MJ && v >= 0 && v <= Mv)
		{
			Js[N] = Tab->item(CJs, n)->text().toInt();
			E[N++] = Tab->item(CWN, n)->text().toDouble() + ELU[0][i][v][J];
		}
	}
	delete[] VR;
}

bool LineTable::ShowUpTerm()
{
	//printf("LineTable::ShowUpTerm\n");
	if (molecule == 0 || transition == 0)
	{
		QMessageBox::information(this, tr("QT4MolSpektAn"), 
						tr("The line table has to be assigned to a molecule and a transition first!"));
		return false;
	}
	ElState *LS = transition->getLowerState();
	if (LS == 0)
	{
		QMessageBox::information(this, tr("QT4MolSpektAn"), 
						 tr("The transition has to be assigned to a lower electronic state first!"));
		return false;
	}
	TermTable *TT = LS->getTermTable();
	DunTable *DT;
	if (TT == 0) 
	{
		DT = LS->getDunTable();
		if (DT != 0) DT->calcTermEnergies(TT);
	}
	if (TT == 0)
	{
		QMessageBox::information(this, tr("QT4MolSpektAn"), 
	tr("A term energy table or a Dunham coefficient set for the lower state has to be loaded first!"));
		return false;
	}
	double ****UT = 0, ****ELU = TT->getData(), mvs = 0, Jeo = 0;
	int MJ = TT->getMaxJ(), Mv = TT->getMaxv(), *CT, MCT, *IsoT = TT->getIsoT();
	int *UIsoT = 0, *UCT = 0, UMCT = 0, UNC = 0;
	float S = 0.0;
	TT->getCompT(MCT, CT);
	ElState *US = transition->getUpperState();
	if (US != 0)
	{
		if (US->getTermTable() == 0 && (DT = US->getDunTable()) != 0) 
			DT->calcTermEnergies(TT = 0);
		if ((TT = US->getTermTable()) != 0) 
		{
			UT = TT->getData();
			mvs = TT->getMaxv();
			Jeo = TT->getMaxJ();
			UIsoT = TT->getIsoT();
			TT->getCompT(UMCT, UCT);
			S = US->getS();
			UNC = TT->getNumComp();
		}
	}
	int n, i, c, I, j, v, N = Tab->rowCount(), NI = (Iso != 0 ? Iso->numIso : 0);
    int **Z = CreateInt(N, 6), li[2] = {-1, -1};
	QString SpektFile, Buffer;
//    bool swaped = true;
    double E, T[N], dE, Calc, ES, EF;
	if (SO != 0) delete[] SO;
	int *LSO = heapSort(isnSPG);
	SO = new int[N];
	for (n=0; n<N; n++) SO[LSO[n]] = n;
	delete[] LSO;
	NSO = N;
	Tab->blockSignals(true);
	if (UT != 0) Tab->setColumnCount(COmC+1);
	else Tab->setColumnCount(CEdJ+1);
	HeaderLabels << "calc. rel. int." << "E_UpTerm" << "E_avarage" << "E_Up - E_av" << "E/(J*(J+1))";
	if (UT != 0) HeaderLabels << "Calc" << "obs-calc";
	Tab->setHorizontalHeaderLabels(HeaderLabels);
	//printf("Vor Schleife1\n");
    for (i=0; i<N; i++)
    {
		if ((Buffer = Tab->item(SO[i], CFile)->text()).isEmpty()) 
			Tab->item(SO[i], CFile)->setText(SpektFile);
		else if (Buffer.left(6) == " laser")
		{
		    for (j=6; Buffer[j] == ' '; j++) ;
	    	SpektFile = Buffer.right(Buffer.length() - j + 1);
		}
		else SpektFile = Buffer;
		j = Tab->item(SO[i], CJss)->text().toInt();
		v = Tab->item(SO[i], Cvss)->text().toInt();
		Z[i][0] = Tab->item(SO[i], CIso)->text().toInt();
		I = (Z[i][0] - 1) / 10;
		Z[i][1] = Tab->item(SO[i], Cvs)->text().toInt();
		Z[i][2] = Tab->item(SO[i], CJs)->text().toInt();
		Z[i][3] = (int)fabs(Z[i][2] - j);
		Z[i][4] = j;
		Z[i][5] = Tab->item(SO[i], CPN)->text().toInt();
		c = Tab->item(SO[i], CF)->text().toInt();
		c = (c >= 0 && c <= MCT ? (CT[c] >= 0 ? CT[c] : 0) : 0);
		//printf("i=%d, I=%d, v=%d, J=%d\n", i, I, v, j);
		E = (((I >= 0 && I < NI ? IsoT[I] >= 0 : false) && v >= 0 && v <= Mv && j >= 0 
				&& j <= MJ ? ELU[c][IsoT[I]][v][j] != 0.0 : false)? 
					Tab->item(SO[i], CWN)->text().toDouble() + ELU[c][IsoT[I]][v][j] : 0.0);
		Tab->setItem(SO[i], CEUp, new QTableWidgetItem(QString::number(E, 'g', 9)));
		T[i] = E;
    }
	//printf("Vor Schleife2\n");
    for (E=ES=0.0, i=n=0; n<=N; n++)
    {
		//printf("E=%f, n=%d, N=%d, i=%d\n", E, n, N, i);
		if (n < N ? (n > i && (Z[i][0] != Z[n][0] || Z[i][1] != Z[n][1] || Z[i][2] != Z[n][2] 
				  || Z[i][3] != Z[n][3] || Z[i][5] != Z[n][5])) : true)
		{
	    	E /= ES;
	    	j = li[Z[n-1][3]];
	    	if ((j>0 ? Z[j][0] == Z[n-1][0] && Z[j][1] == Z[n-1][1] && 
					   Z[j][2] == Z[n-1][2] - 1 : false)) dE = (E - T[j]) / (2 * Z[n-1][2]);
	    	else dE = 0.0;
	    	for (j=i; j<n; j++) 
	    	{
				/*if (n==1482) 
					printf("Z[%d][0]=%d, Z[%d][1]=%d, Z[%d][2]=%d, Z[%d][4]=%d\n", 
						   j, Z[j][0], j, Z[j][1], j, Z[j][2], j, Z[j][4]);*/
				Tab->setItem(SO[j], CEav, 
							 new QTableWidgetItem(QString::number(E, 'g', 9)));
				Tab->setItem(SO[j], CEUma, 
							 new QTableWidgetItem(QString::number(T[j] - E, 'g', 6)));
				if (dE > 0.0) 
					Tab->setItem(SO[j], CEdJ, 
								 new QTableWidgetItem(QString::number(dE, 'g', 6)));
				I = (Z[j][0] - 1) / 10;
				if (T[j] != 0.0 && UT != 0 && Z[j][1] >= 0 && Z[j][1] <= mvs 
					&& Z[j][2] <= Jeo && UIsoT[I] >= 0)
				{
					//printf("Beginn Neu, j=%d\n", j);
					v = Tab->item(SO[j], Cvss)->text().toInt();
					
					//printf("I=%d, Z[j][0]=%d, Z[j][1]=%d, Z[j][2]=%d, Z[j][4]=%d, v=%d\n", 
						//   I, Z[j][0], Z[j][1], Z[j][2], Z[j][4], v);
					//printf("ELU=%f, eT=%f\n", ELU[I][v][Z[j][4]], eT[I][Z[j][1]][Z[j][2]]);
					if (S==0.0) c = (Z[j][2] != Z[j][4] || UNC == 1 ? 0 : 1); 
					else
					{
						c = Tab->item(SO[i], CF)->text().toInt();
						c = (c >= 0 && c <= UMCT ? (UCT[c] >= 0 ? UCT[c] : 0) : 0);
					}
					Calc = UT[c][UIsoT[I]][Z[j][1]][Z[j][2]];
					//printf("E=%f, Calc=%f\n", E, Calc);
					Tab->setItem(SO[j], CCalc, 
								 new QTableWidgetItem(QString::number(Calc, 'g', 9)));
					Tab->setItem(SO[j], COmC, 
								 new QTableWidgetItem(QString::number(E - Calc, 'g', 6)));
					//printf("Ende Neu\n");
				}
	    	}
	    	//if (n==1482) printf("n=%d, i=%d\n", n, i);
	    	if (n < N)
	    	{
				li[Z[n-1][3]] = n - 1;
				i = n;
				T[n-1] = E;
				E = ES = 0.0;
	    	}
		}
		if (n<N) 
		{
			EF = Tab->item(SO[n], Cerr)->text().toDouble();
			EF = 1.0 / (EF * EF);
			E += T[n] * EF;
			ES += EF;
		}
    }
	//printf("Vor Destroy\n");
	Tab->blockSignals(false);
    Destroy(Z, N);
	//printf("Ende ShowUpTerm\n");
	Changed();
	return true;
}

void LineTable::FindBigDiff()
{
	bool ChDD = false;
	QList<QTableWidgetSelectionRange> Sel = Tab->selectedRanges();
	int i, n = Sel.size(), N = Tab->rowCount(), nc = Tab->columnCount();
	//printf("nc=%d, COmC=%d\n", nc, COmC);
	if (NSO != N) ShowUpTerm();
	QString Buffer, File, FBuffer;
	//printf("nc=%d, COmC=%d\n", nc, COmC);
    if (++lRow >= N - 1) lRow = 0;
    for (i=0; i<n; i++) if (lRow >= Sel[i].topRow() && lRow <= Sel[i].bottomRow())
    {
		Buffer = Tab->item(lRow, CEav)->text();
		File = Tab->item(lRow, CFile)->text();
    }
	//printf("Buffer=%s, File=%s\n", Buffer.ascii(), File.ascii());
    while ((lRow<N ? (fabs(Tab->item(lRow, CEUma)->text().toDouble()) < 0.03 
			&& (nc >= 14 && ChDD ? fabs(Tab->item(lRow, COmC)->text().toDouble()) < 0.1 : true)) 
			|| (Tab->item(lRow, CEav)->text() == Buffer && Tab->item(lRow, CFile)->text() == File)
			: false))
		lRow++;
    if (lRow < N) 
    {
		for (i=0; i<n; i++) Tab->setRangeSelected(Sel[i], false);
		File = Tab->item(lRow, CFile)->text();
		if (File.left(6) == " laser") 
		{	
	    	for (i=6; File[i] == ' '; i++) ;
	    	File = File.right(File.length() - i);
		}
		else File = File.right(File.length() - 1);
		/*printf("File=%s\n", File.ascii());
		if (File != InFile) 
		{
		    emit OpenRequest(File);
	    	raise();
		    setActiveWindow();
	    	setFocus();
	    	printf("InFile=%s\n", InFile.ascii());
		}*/
		for (n=lRow, Buffer=Tab->item(lRow, CEav)->text(); Buffer == Tab->item(n, CEav)->text(); n--) ;
		n++;
		//printf("Buffer=%s, Tab->text(n, 9)=%s\n", Buffer.ascii(), Tab->text(n, CEav).ascii());
		//Tab->ensureCellVisible(lRow, CEUma);
		while (Tab->item(n, CEav)->text() == Buffer)
		{
		    if ((FBuffer = Tab->item(n, CFile)->text()).left(6) == " laser") 
				for (i=6; FBuffer[i] == ' '; i++) ;
	    	else i = 1;
		    FBuffer = FBuffer.right(FBuffer.length() - i);
		    if (FBuffer == File) Tab->selectRow(n);
		    //printf("FBuffer=%s, File=%s\n", FBuffer.ascii(), File.ascii());
	    	n++;
		}
		//MarkSelected();
    }
    //printf("lRow = %d\n", lRow);
}

void LineTable::ShowWeakProgressions()
{
	int NRow = Tab->rowCount(), fRow = 0, aRow, av = -1, aJ = -1, lv = 0, lJ = 0;
	QString aFile, lFile;
	bool Weak = false;
	if (lRow >= NRow) aRow = 0;
	else aRow = lRow + 1;
	for (;aRow != lRow; aRow++)
	{
		if (aRow==NRow) aRow = 0;
		aFile = Tab->item(aRow, CFile)->text();
		aJ = Tab->item(aRow, CJs)->text().toInt();
		av = Tab->item(aRow, Cvs)->text().toInt();
		if (av != lv || aJ != lJ || aFile != lFile)
		{
			if (Weak)
			{
				for (lRow = fRow; lRow < aRow; lRow++) Tab->selectRow(lRow);
				//Tab->ensureCellVisible(lRow, CSNR);
				return;
			}
			else Weak = true;
			lJ = aJ;
			lv = av;
			lFile = aFile;
			fRow = aRow;
		}
		else if ((Tab->item(aRow, CSNR)->text().toDouble() > 4 && Tab->item(aRow, CC)->text().isEmpty())
				|| Tab->item(aRow, CC)->text() == "satellite") Weak = false;
	}
	QMessageBox::information( this, "MolSpektAnalysis", 
			"There are no really weak progressions found inside the list.", QMessageBox::Ok);
}

TableWindow* LineTable::ShowUpTermTable()
{
    if (MW == 0) return 0;
	int i, j=0, k=0, n, iJs = 0, AD = 0, N = 0, iJss;
    double Diff, SqSum = 0.0, TV = 0.0, TE, MinE = 0.0, MaxE = 0.0;
    QString Buffer, Js, Jss, vs, vss, Iso, File, J, **Data;
    if (Tab->columnCount() == TableNormCols || NSO != Tab->rowCount()) 
		if (!ShowUpTerm()) return 0;
    if (termTable == NULL) 
    {
		termTable = new TableWindow(TextTable1, MW, molecule);
		termTable->setWindowTitle("UpTermTable to " + getName());
		termTable->setHorizontalHeader(QStringList() << "Isotop" << "v" << "J" << "Par" 
				<< "Termenergie" << "Standardabweichung" << "Anzahl Uebergaenge" 
				<< "Anz. vollst. Doubl." << "Min E" << "Max E");
    }
    for (i=0; i < NSO; i++) if (Tab->item(SO[i], CEav)->text() != Buffer)
    {
		Buffer = Tab->item(SO[i], CEav)->text();
		N++;
    }
    Data = CreateQString(NSO, 10);
    for (n=i=0; i <= NSO; i++) 
    {
		if ((i < NSO ? Tab->item(SO[i], CEav)->text() != Buffer : true))
		{
	    	if (i > 0)
	    	{
				Data[n][0] = Iso;
				Data[n][1] = vs;
				Data[n][2] = Js;
				if (Jss == Js) Data[n][3] = "f";
				else Data[n][3] = "e";
				Data[n][4] = Buffer;
				if (j > 1) 
					Data[n][5] = QString::number(sqrt(SqSum / (j * (j - 1))), 'g', 6);
				Data[n][6] = QString::number(j);
				if (Jss != Js) Data[n][7] = QString::number(AD);
				Data[n][8] = QString::number(MinE, 'g', 9);
				Data[n++][9] = QString::number(MaxE, 'g', 9);
	    	}
	    	if (i < NSO)
	    	{
				TV = (Buffer = Tab->item(SO[i], CEav)->text()).toDouble();
				iJs = (Js = Tab->item(SO[i], CJs)->text()).toInt();
				Jss = Tab->item(SO[i], CJss)->text();
				Iso = Tab->item(SO[i], CIso)->text();
				vs = Tab->item(SO[i], Cvs)->text();
				j = 0;
				SqSum = 0.0;
				MinE = MaxE = TV;
				AD = 0;
	    	}
		}
		if (i < NSO)
		{
			Diff = TV - (TE = Tab->item(SO[i], CEUp)->text().toDouble());
			if (TE < MinE) MinE = TE;
			else if (TE > MaxE) MaxE = TE;
			SqSum += Diff * Diff;
	    	j++;
			if (Jss != Js)
			{
				iJss = 2 * iJs - Jss.toInt();
				vss = Tab->item(SO[i], Cvss)->text();
				File = Tab->item(SO[i], CFile)->text();
				//printf("Jss=%s, vss=%s, vs=%s, Js=%s, File=%s\n", 
					 //  Jss.ascii(), vss.ascii(), vs.ascii(), Js.ascii(), File.ascii());
				for (k=i-1; (k>=0 ? Tab->item(SO[k], CJs)->text() == Js 
								 && Tab->item(SO[k], Cvs)->text() == vs 
								 && Tab->item(SO[k], CIso)->text() == Iso : false); k--) 
					if (Tab->item(SO[k], Cvss)->text().toInt() == vss.toInt() 
						&& Tab->item(SO[k], CJss)->text().toInt() 
								== 2 * Js.toInt() - iJss 
						&& Tab->item(SO[k], CFile)->text() == File) AD++;
			}
			//printf("k=%d, i=%d, AD=%d\n", k, i, AD);
		}
    }
    termTable->setData(Data, n, 10);
	Destroy(Data, NSO);
    return termTable;
}

void LineTable::Updatevs(int *nvs)
{
	int i, nr = Tab->rowCount();
	QString Buffer;
	Tab->blockSignals(true);
	for (i=0; i<nr; i++) 
	{
		Buffer = "   " + QString::number(nvs[i]);
		Tab->item(i, Cvs)->setText(Buffer.right(4));
	}
	Tab->blockSignals(false);
	Changed();
}

void LineTable::RemoveDoubled()
{
	int i, j, nr = Tab->rowCount(), c, C = Tab->columnCount(), *S1 = heapSort(sortfRemDoubl);
	int n1, n2, m1, m2, S[nr];
	QString F1, F2, Comment;
	for (i=0; i < nr; i++) S[S1[i]] = i;
	Tab->blockSignals(true);
	for (i=1, F2 = Tab->item(S[0], CFile)->text(); i < nr; i++)
	{
		F1 = F2;
		F2 = Tab->item(S[i], CFile)->text();
		//printf("I=%d, vs=%d, Js=%d, vss=%d, Jss=%d, WN=%f\n", Tab->item(S[i], CIso)->text().toInt(), Tab->item(S[i], Cvs)->text().toInt(),
			//   Tab->item(S[i], CJs)->text().toInt(), Tab->item(S[i], Cvss)->text().toInt(), Tab->item(S[i], CJss)->text().toInt(),
			  // Tab->item(S[i], CWN)->text().toDouble());
		if (Tab->item(S[i], Cvss)->text().toInt() == Tab->item(S[i-1], Cvss)->text().toInt()) 
			if (Tab->item(S[i], CJss)->text().toInt() == Tab->item(S[i-1], CJss)->text().toInt())
			if (Tab->item(S[i], CJs)->text().toInt() == Tab->item(S[i-1], CJs)->text().toInt())
		{
			n1 = F1.lastIndexOf(QRegExp("[\\/]"));
			n2 = F2.lastIndexOf(QRegExp("[\\/]"));
			m1 = ((m1 = F1.lastIndexOf('.')) != -1 ? m1 : F1.length());
			m2 = ((m2 = F2.lastIndexOf('.')) != -1 ? m2 : F2.length());
			if (F1.mid(n1 + 1, m1 - n1 - 1) == F2.mid(n2 + 1, m2 - n2 - 1)
				&& Tab->item(S[i], CIso)->text().toInt() == Tab->item(S[i-1], CIso)->text().toInt()
				&& fabs(Tab->item(S[i], CWN)->text().toDouble() - Tab->item(S[i-1], CWN)->text().toDouble()) < 1e-4)
			{
				Tab->setItem(S[i], 0, 0);
				if (F1 != F2)
				{
					QFile File1(F1), File2(F2);
					if (!File1.exists() && File2.exists()) Tab->item(S[i-1], CFile)->setText(F2);
				}
				Comment = Tab->item(S[i], CC)->text();
				if (Tab->item(S[i-1], CC)->text().length() < Comment.length()) Tab->item(S[i-1], CC)->setText(Comment);
				if (Tab->item(S[i-1], Cerr)->text().toDouble() < Tab->item(S[i], Cerr)->text().toDouble())
					Tab->item(S[i-1], Cerr)->setText(Tab->item(S[i], Cerr)->text());
			}
		}
	}
	for (j=0; (j < nr ? Tab->item(j, 0) != 0 : false); j++) ;
	for (i=j+1; i < nr; i++) if (Tab->item(i, 0) != 0) 
	{
		for (c=0; c<C; c++) Tab->setItem(j, c, Tab->takeItem(i, c));
		j++;
	}
	Tab->setRowCount(j);
	Tab->blockSignals(false);
	Changed();
	delete[] S1;
} 

void LineTable::setData(int nR, int nC, QStringList &H, QTableWidgetItem ***D)
{
	int r, c;
	Tab->blockSignals(true);
	Tab->setRowCount(nR);
	Tab->setColumnCount(nC);
	Tab->setHorizontalHeaderLabels(H);
	for (r=0; r < nR; r++) for (c=0; c < nC; c++) Tab->setItem(r, c, D[r][c]);
	Tab->blockSignals(false);
	Changed();
}

void LineTable::splitTable()
{
	if (MW == 0) return;
	LineTable *NT = MW->CreateLineTable();
	if (NT == 0) return;
	int i, j, k, c, R = Tab->rowCount(), C = Tab->columnCount();
	QTableWidgetItem ***IS = new QTableWidgetItem**[R];
	QStringList L;
	ElState *lS, *uS;
	if (transition != 0) 
	{
		lS = transition->getLowerState();
		uS = transition->getUpperState();
	}
	else lS = uS = 0;
	Tab->blockSignals(true);
	for (k=0; (k<R ? Tab->item(k, Cvs)->text().toInt() != -1 : false); k++) ;
	for (j=k, i=0; k<R; k++)
	{
		if (Tab->item(k, Cvs)->text().toInt() == -1)
		{
			IS[i] = new QTableWidgetItem*[C];
			for (c=0; c<C; c++) IS[i][c] = Tab->takeItem(k, c);
			i++;
		}
		else 
		{
			for (c=0; c<C; c++) Tab->setItem(j, c, Tab->takeItem(k, c));
			j++;
		}
	}
	Tab->setRowCount(j);
	Tab->blockSignals(false);
	for (c=0; c<C; c++) L << Tab->horizontalHeaderItem(c)->text();
	NT->setData(i, C, L, IS);
	NT->setSource(getSource());
	NT->setName(getName());
	if (molecule != 0) molecule->addLineTable(NT, lS);
	NT->show();
	for (j=0; j<i; j++) delete[] IS[j];
	delete[] IS;
	Changed();
}

void LineTable::sortUpTermIvJ()
{
	if (Tab->columnCount() < CEav)
	{
		QMessageBox::information( this, "MolSpektAnalysis", 
			"There is no upper term energy table calculated for this line table!", QMessageBox::Ok);
		return;
	}
	/*bool S = true;
	int n, I, Iv, v, vv, J, Jv, N = termTable->rowCount(), c, C = termTable->columnCount();
	QTableWidgetItem *Item;
	while (S)
	{
		S = false;
		Iv = termTable->item(0, 0)->text().toInt();
		vv = termTable->item(0, 1)->text().toInt();
		Jv = termTable->item(0, 2)->text().toInt();
		for (n = 1; n  < N; n++) 
		{
			I = termTable->item(n, 0)->text().toInt();
			v = termTable->item(n, 1)->text().toInt();
			J = termTable->item(n, 2)->text().toInt();
			if (Iv > I || (Iv == I && (vv > v || (vv == v && Jv > J))))
			{
				for (c=0; c<C; c++)
				{
					Item = termTable->takeItem(n-1, c);
					termTable->setItem(n-1, c, termTable->takeItem(n, c));
					termTable->setItem(n, c, Item);
				}
				S = true;
			}
			else 
			{
				Iv = I;
				vv = v;
				Jv = J;
			}
		}
	}*/
	int *SArray = heapSort(sortUtIvJ);
	sortTab(SArray);
}

void LineTable::SortfRemDoubled()
{
	sortTab(heapSort(sortfRemDoubl));
}

void LineTable::SortProg()
{
	SortIJvP();
}

void LineTable::SortIJvP()
{
	int *SArray = heapSort(sortIJvP);
	sortTab(SArray);
}

void LineTable::SortIvPJ()
{
	int *SArray = heapSort(sortIvPJ);
	sortTab(SArray);
}

void LineTable::SortFPInt()
{
	int *SArray = heapSort(sortFPInt);
	sortTab(SArray);
}

void LineTable::SortSpectrum()
{
	int *SArray = heapSort(sortBySpectrum);
	sortTab(SArray);
}

void LineTable::sortTab(int *S2)
{
	Tab->blockSignals(true);
	QString F1, F2, FN;
	int l, i, n, N = Tab->rowCount();
	for (n=0; n<N; n++)
	{
		l = (F2 = Tab->item(n, CFile)->text()).length();
		if (l==0) F2 = F1;
		else
		{
			for (i=0; (i < l ? F2[i] == QChar(' ') : false); i++) ;
			F2 = F2.right(l-=i);
			if (F2.left(5) == "laser")
			{
				F2 = F2.right(l-=5);
				for (i=0; (i < l ? F2[i] == QChar(' ') : false); i++) ;
				F2 = F2.right(l-i);
				if (F2.isEmpty()) F2 = F1;
				Tab->item(n, CC)->setText("laser");
			}
			F1 = F2;
		}
		Tab->item(n, CFile)->setText(F2);
	}
	Tab->blockSignals(false);
	TableWindow::sortTab(S2);
}

void LineTable::SetError()
{
	QDialog SDialog(this);
	SDialog.setWindowTitle("Set error");
	QGridLayout *L = new QGridLayout(&SDialog);
	QCheckBox *MaxMin = new QCheckBox("Only max/min values", &SDialog);
	QLabel L1("Please specify the error to\nbe set for the selected lines.", &SDialog);
	L->addWidget(&L1, 0, 0, 1, 2);
	L->addWidget(MaxMin, 1, 0, 1, 2);
	QLabel L2("Good lines:", &SDialog);
	L->addWidget(&L2, 2, 0);
	QLineEdit NE(QString::number(Error, 'g', 5), &SDialog);
	L->addWidget(&NE, 2, 1);
	NE.setAlignment(Qt::AlignRight);
	QLabel L3("Overlapping lines:", &SDialog);
	L->addWidget(&L3, 3, 0);
	QLineEdit OE(QString::number(OvError, 'g', 5), &SDialog);
	L->addWidget(&OE, 3, 1);
	OE.setAlignment(Qt::AlignRight);
	QPushButton OK("OK", &SDialog), Cancel("Cancel", &SDialog);
	L->setRowMinimumHeight(4, 20);
	L->addWidget(&OK, 5, 0);
	L->addWidget(&Cancel, 5, 1);
	connect(&OK, SIGNAL(clicked()), &SDialog, SLOT(accept()));
	connect(&Cancel, SIGNAL(clicked()), &SDialog, SLOT(reject()));
	if (SDialog.exec() == QDialog::Rejected) return;
	Error = NE.text().toDouble();
	OvError = OE.text().toDouble();
	QString SErr = ("          " + QString::number(Error, 'g', 5)).right(10), Buffer;
	QString OErr = ("          " + QString::number(OvError, 'g', 5)).right(10);
	int i, n = Tab->rowCount(), r, rc, m=0;
	QList<QTableWidgetSelectionRange> SR = Tab->selectedRanges();
	Tab->blockSignals(true);
	if (MaxMin->isChecked())
	{
		if ((rc = SR.size()) > 0) for (r=0; r < rc; r++) for (i = SR[r].topRow(), n = SR[r].bottomRow(); i<=n; i++)
		{
			Buffer = Tab->item(i, CC)->text();
			if (Buffer.indexOf("overlap") > -1 || Buffer.indexOf("laser") > -1 || Buffer.indexOf("weak") > -1)
			{
				if (Tab->item(i, Cerr)->text().toDouble() < OvError) Tab->item(i, Cerr)->setText(OErr); 
			}
			else if (Tab->item(i, Cerr)->text().toDouble() > Error) Tab->item(i, Cerr)->setText(SErr);
		}
		else for (i=0, n = Tab->rowCount(); i<n; i++)
		{
			Buffer = Tab->item(i, CC)->text();
			if (Buffer.indexOf("overlap") > -1 || Buffer.indexOf("laser") > -1 || Buffer.indexOf("weak") > -1)
			{
				if (Tab->item(i, Cerr)->text().toDouble() < OvError) Tab->item(i, Cerr)->setText(OErr); 
			}
			else if (Tab->item(i, Cerr)->text().toDouble() > Error) Tab->item(i, Cerr)->setText(SErr);
		}
	}
	else for (r=0, rc = SR.size(); r < rc; r++)
	{
		m = SR[r].topRow();
		n = SR[r].bottomRow();
		for (i=m; i<=n; i++)
		{
			Buffer = Tab->item(i, CC)->text();
			if (Buffer.indexOf("overlap") > -1 || Buffer.indexOf("laser") > -1 
						 || Tab->item(i, CSNR)->text().toDouble() < 5.0) 
				Tab->item(i, Cerr)->setText(OErr) ;
			else Tab->item(i, Cerr)->setText(SErr);
		}
	}
	Tab->blockSignals(false);
	Changed();
}

void LineTable::Shiftvup()
{
	int i, m, n, r, rc;
	QList<QTableWidgetSelectionRange> SR = Tab->selectedRanges();
	Tab->blockSignals(true);
	for (r=0, rc = SR.size(); r < rc; r++)
	{
		m = SR[r].topRow();
		n = SR[r].bottomRow();
		for (i=m; i<=n; i++) 
			Tab->item(i, Cvss)->setText(QString::number(Tab->item(i, Cvss)->text().toInt() + 1));
	}
	Tab->blockSignals(false);
	Changed();
}

void LineTable::Shiftvdown()
{
	int i, m, n, r, rc;
	QList<QTableWidgetSelectionRange> SR = Tab->selectedRanges();
	Tab->blockSignals(true);
	for (r=0, rc = SR.size(); r < rc; r++)
	{
		m = SR[r].topRow();
		n = SR[r].bottomRow();
		for (i=m; i<=n; i++) 
			Tab->item(i, Cvss)->setText(QString::number(Tab->item(i, Cvss)->text().toInt() - 1));
	}
	Tab->blockSignals(false);
	Changed();
}

void LineTable::ShiftJup()
{
	int i, m, n, r, rc;
	QList<QTableWidgetSelectionRange> SR = Tab->selectedRanges();
	Tab->blockSignals(true);
	for (r=0, rc = SR.size(); r < rc; r++)
	{
		m = SR[r].topRow();
		n = SR[r].bottomRow();
		for (i=m; i<=n; i++) 
		{
			Tab->item(i, CJs)->setText(QString::number(Tab->item(i, CJs)->text().toInt() + 1));
			Tab->item(i, CJss)->setText(QString::number(Tab->item(i, CJss)->text().toInt() + 1));
		}
	}
	Tab->blockSignals(false);
	Changed();
}

void LineTable::ShiftJdown()
{
	int i, m, n, r, rc;
	QList<QTableWidgetSelectionRange> SR = Tab->selectedRanges();
	Tab->blockSignals(true);
	for (r=0, rc = SR.size(); r < rc; r++)
	{
		m = SR[r].topRow();
		n = SR[r].bottomRow();
		for (i=m; i<=n; i++) 
		{
			Tab->item(i, CJs)->setText(QString::number(Tab->item(i, CJs)->text().toInt() - 1));
			Tab->item(i, CJss)->setText(QString::number(Tab->item(i, CJss)->text().toInt() - 1));
		}
	}
	Tab->blockSignals(false);
	Changed();
}

void LineTable::Shiftvsdown()
{
	int i, m, n, r, rc;
	QList<QTableWidgetSelectionRange> SR = Tab->selectedRanges();
	Tab->blockSignals(true);
	for (r=0, rc = SR.size(); r < rc; r++)
	{
		m = SR[r].topRow();
		n = SR[r].bottomRow();
		for (i=m; i<=n; i++) 
			Tab->item(i, Cvs)->setText(QString::number(Tab->item(i, Cvs)->text().toInt() - 1));
	}
	Tab->blockSignals(false);
	Changed();
}

void LineTable::Shiftvsup()
{
	int i, m, n, r, rc;
	QList<QTableWidgetSelectionRange> SR = Tab->selectedRanges();
	Tab->blockSignals(true);
	for (r=0, rc = SR.size(); r < rc; r++)
	{
		m = SR[r].topRow();
		n = SR[r].bottomRow();
		for (i=m; i<=n; i++) 
			Tab->item(i, Cvs)->setText(QString::number(Tab->item(i, Cvs)->text().toInt() + 1));
	}
	Tab->blockSignals(false);
	Changed();
}

void LineTable::ShiftIso()
{
	IsoTab *Iso = (molecule != 0 ? molecule->getIso() : 0);
	int M = ((Iso != 0 ? Iso->numIso : 0) - 1) * 10;
	if (M == -10)
	{
		QMessageBox::information(this, "MolSpektAnalysis", 
								 "This line table has to be assigned to a molecule first!",
								 QMessageBox::Ok);
		return;
	}
	int i, n, m, r, rc, v;
	QList<QTableWidgetSelectionRange> SR = Tab->selectedRanges();
	Tab->blockSignals(true);
	for (r=0, m=0, rc = SR.size(); r < rc; r++)
	{
		m = SR[r].topRow();
		n = SR[r].bottomRow();
		for (i=m; i<=n; i++) 
		{
			if ((v = Tab->item(i, CIso)->text().toInt()) <= M) 
				Tab->item(i, CIso)->setText(QString::number(v + 10));
			else Tab->item(i, CIso)->setText(QString::number(v - M));
		}
	}
	Tab->blockSignals(false);
	Changed();
}

void LineTable::SetvssAscending()
{
	QList<QTableWidgetSelectionRange> SR = Tab->selectedRanges();
	int i, r = SR.size(), l=-1, n, m;
	if (r==0)
	{
		QMessageBox::information(this, "MolSpektAnalysis", "You have to select some rows first!",
								 QMessageBox::Ok);
		return;
	}
	Tab->blockSignals(true);
	for (i=0; i<r; i++)
	{
		n = SR[i].bottomRow();
		for (m = SR[i].topRow(); m<=n; m++)
		{
			if (l==-1) l = Tab->item(i, Cvss)->text().toInt() + 1;
			else Tab->item(i, Cvss)->setText(QString::number(l++));
		}
	}
	Tab->blockSignals(false);
	Changed();
}

void LineTable::Delete()
{
	QList<QTableWidgetItem*> SI = Tab->selectedItems();
	int i, n = SI.size();
	Tab->blockSignals(true);
	for (i=0; i<n; i++) SI[i]->setText("");
	Tab->blockSignals(false);
	Changed();
}

void LineTable::cutRows(int& numRows, int& numColumns, QString**& Data)
{
    if (NpL > 0)
	{
		int NR = Tab->rowCount();
		int SR = NR - NpL, n, r;
		QList<QTableWidgetSelectionRange> SelR = Tab->selectedRanges();
		bool Rem[NpL];
		for (n=0; n < NpL; n++) Rem[n] = false;
		for (n=0; n < SelR.size(); n++) for (r = SelR[n].topRow(); r <= SelR[n].bottomRow(); r++) if (r >= SR) Rem[r - SR] = true;
		for (n=r=0; n < NpL; n++) if (Rem[n]) r++;
		NpL -= r;
	}
	TableWindow::cutRows(numRows, numColumns, Data);
}

void LineTable::DeleteRows()
{
    if (Tab == 0) return;
	QList<QTableWidgetSelectionRange> SR = Tab->selectedRanges();
	int j, k, n = Tab->rowCount(), c, nc = Tab->columnCount();
	int pr = n - NpL;
	QString B, Buffer;
	Tab->blockSignals(true);
	for (j=0; j < SR.size(); j++) 
	{
		for (k=SR[j].topRow(); k <= SR[j].bottomRow(); k++)
		{
			if (Tab->item(k, CFile) != 0 ? 
				!(B = Tab->item(k, CFile)->text()).isEmpty() : false) Buffer = B;
			for (c=0; c < nc; c++) Tab->setItem(k, c, 0);
			if (k >= pr) NpL--;
		}
		if (k < n ? (Tab->item(k, CFile) != 0 ? Tab->item(k, CFile)->text().isEmpty() 
				: false) : false) 
			Tab->item(k, CFile)->setText(Buffer);
	}
	for (j=k=0; j <= lRow; j++) if (Tab->item(j, 0) == 0) k++;
	lRow -= k;
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

void LineTable::deleteRows(int* rows, int N)
{
	if (Tab == 0) return;
	int j, k, n = Tab->rowCount(), c, nc = Tab->columnCount();
	int pr = n - NpL;
	QString B, Buffer;
	Tab->blockSignals(true);
	for (j=0; j<N; j++)
	{
		if (Tab->item(rows[j], CFile) != 0 ? !(B = Tab->item(rows[j], CFile)->text()).isEmpty() : false) Buffer = B;
		for (c=0; c < nc; c++) Tab->setItem(rows[j], c, 0);
		if (rows[j] > pr) NpL--;
		if ((k = rows[j] + 1) < n ? (Tab->item(k, CFile) != 0 ? Tab->item(k, CFile)->text().isEmpty() 
				: false) : false)
			Tab->item(k, CFile)->setText(Buffer);
	}
	for (j=k=0; j <= lRow; j++) if (Tab->item(j, 0) == 0) k++;
	lRow -= k;
	for (j=0; (j<n ? Tab->item(j, 0) != 0 : false); j++) ;
	for (k=j; k<n; k++) if (Tab->item(k, 0) != 0)
	{
		for (c=0; c < nc; c++) Tab->setItem(j, c, Tab->takeItem(k, c));
		j++;
	}
	Tab->setRowCount(j);
	Tab->blockSignals(false);
	Changed();
}

void LineTable::DeleteRows(int NL, int* PN, int* vss, int* Jss)
{
	int R = Tab->rowCount(), r, C = Tab->columnCount(), c, n, sr = R - NpL;
	Tab->blockSignals(true);
	for (r=0; r<R; r++)
	{
		for (n=0; (n < NL ? Tab->item(r, CPN)->text().toInt() != PN[n] 
					|| (vss != 0 && Jss != 0 ? Tab->item(r, Cvss)->text().toInt() != vss[n] 
							|| Tab->item(r, CJss)->text().toInt() != Jss[n] : false) : false); n++) ;
		if (n < NL)
		{
			Tab->setItem(r, 0, 0);
			if (r >= sr) NpL--;
		}
	}
	for (r=0; (r<R ? Tab->item(r, 0) != 0 : false); r++) ;
	for (n=r+1; n<R; n++) if (Tab->item(n, 0) != 0)
	{
		for (c=0; c<C; c++) Tab->setItem(r, c, Tab->takeItem(n, c));
		r++;
	}
	Tab->blockSignals(false);
	Changed();
}

void LineTable::SetError(int NL, int* PN, int* vss, int* Jss, double* Err)
{
	int R = Tab->rowCount(), r, n;
	Tab->blockSignals(true);
	for (r=0; r<R; r++)
	{
		for (n=0; (n < NL ? Tab->item(r, CPN)->text().toInt() != PN[n] || Tab->item(r, Cvss)->text().toInt() != vss[n]
							|| Tab->item(r, CJss)->text().toInt() != Jss[n] : false); n++) ;
		if (n < NL) Tab->item(r, Cerr)->setText(QString::number(Err[n], 'f', 3));
	}
	Tab->blockSignals(false);
	Changed();
}

void LineTable::SetError(int NL, int* PN, int* vss, int* Jss, QString* Err)
{
	int R = Tab->rowCount(), r, n;
	Tab->blockSignals(true);
	for (r=0; r<R; r++)
	{
		for (n=0; (n < NL ? Tab->item(r, CPN)->text().toInt() != PN[n] || Tab->item(r, Cvss)->text().toInt() != vss[n]
							|| Tab->item(r, CJss)->text().toInt() != Jss[n] : false); n++) ;
		if (n < NL) Tab->item(r, Cerr)->setText(Err[n]);
	}
	Tab->blockSignals(false);
	Changed();
}

Progression LineTable::getSelectedProgression()
{
	Progression P;
	int L[1000], r = Tab->currentRow(), R = Tab->rowCount(), n;
	if (r == -1)
	{
		P.N = 0;
		return P;
	}
	int I = P.Iso = Tab->item(r, CIso)->text().toInt();
	int J = P.Js = Tab->item(r, CJs)->text().toInt();
	int v = P.vs = Tab->item(r, Cvs)->text().toInt(), PN = Tab->item(r, CPN)->text().toInt();
	QString F = Tab->item(r, CFile)->text();
	for (r=n=0; r<R; r++) if (Tab->item(r, CPN)->text().toInt() == PN)
		if (Tab->item(r, CFile)->text() == F) 
			if (Tab->item(r, CIso)->text().toInt() == I 
					&& Tab->item(r, Cvs)->text().toInt() == v
					&& Tab->item(r, CJs)->text().toInt() == J) L[n++] = r;
	P.L = new Line[P.N = n];
	for (n=0; n < P.N; n++)
	{
		P.L[n].E = Tab->item(L[n], CWN)->text().toDouble();
		P.L[n].Jss = Tab->item(L[n], CJss)->text().toInt();
		P.L[n].err = Tab->item(L[n], Cerr)->text().toDouble();
		P.L[n].vss = Tab->item(L[n], Cvss)->text().toInt();
	}
	return P;
}

void LineTable::findSimilarProgression(Progression P)
{
	if (P.N == 0) 
	{
		printf("LineTable::findSimilarProgression error: no Progression selected\n");
		return;
	}
	int n, mv, Mv = P.L[0].vss + 1, N, LM, M = (P.Js == P.L[0].Jss ? P.N - 1 : P.N / 2 - 1);
	int v[1000], vs, Js, I, PN, RC = Tab->rowCount(), cr, m, nwL;
	double S[1000], B;
	QString F;
	for (n=0, mv = Mv - 2; n < P.N; n++)
	{
		if (P.L[n].vss <= mv) mv = P.L[n].vss - 1;
		else if (P.L[n].vss >= Mv) Mv = P.L[n].vss + 1;
	}
	n = cr = Tab->currentRow();
	F = Tab->item(n, CFile)->text();
	vs = Tab->item(n, Cvs)->text().toInt();
	Js = Tab->item(n, CJs)->text().toInt();
	I = Tab->item(n, CIso)->text().toInt();
	PN = Tab->item(n, CPN)->text().toInt();
	for (n++; (n < RC ? Tab->item(n, CFile)->text() == F 
			&& Tab->item(n, Cvs)->text().toInt() == vs
			&& Tab->item(n, CJs)->text().toInt() == Js
			&& Tab->item(n, CIso)->text().toInt() == I
			&& Tab->item(n, CPN)->text().toInt() == PN : false); n++) ;
	for (F = Tab->item(n, CFile)->text(), 
		 vs = Tab->item(n, Cvs)->text().toInt(),
		 Js = Tab->item(n, CJs)->text().toInt(),
		 I = Tab->item(n, CIso)->text().toInt(),
		 PN = Tab->item(n, CPN)->text().toInt(), N=0; n != cr; n++)
	{
		if (n == RC) 
		{
			n=0;
			if (cr == 0) break;
		}
		if (Tab->item(n, CFile)->text() != F 
			|| Tab->item(n, Cvs)->text().toInt() != vs
			|| Tab->item(n, CJs)->text().toInt() != Js
			|| Tab->item(n, CIso)->text().toInt() != I
			|| Tab->item(n, CPN)->text().toInt() != PN)
		{
			LM = (Tab->item((n!=0 ? n : RC) - 1, CJss)->text().toInt() == Js ? M : 2 * M); 
			for (vs = 0; vs != -1; ) for (vs = -1, m = 1; m < N; m++) if (S[m-1] < S[m])
			{
				B = S[m-1];
				S[m-1] = S[m];
				S[m] = B;
				vs = v[m-1];
				v[m-1] = v[m];
				v[m] = vs;
			}
			for (m = nwL = 0; m < LM && m < N; m++)
			{
				vs = Tab->item(v[m], Cvss)->text().toInt();
				F = Tab->item(v[m], CC)->text();
				if ((vs < mv || vs > Mv) && F.indexOf("overlap") == -1
											&& F.indexOf("laser") == -1) nwL++;
			}
			if (nwL <= 1)
			{
				vs = (n!=0 ? n : RC);
				for (m = vs - N; m < vs; m++) Tab->selectRow(m);
				return;
			}
			vs = Tab->item(n, Cvs)->text().toInt();
			Js = Tab->item(n, CJs)->text().toInt();
			F = Tab->item(n, CFile)->text();
			I = Tab->item(n, CIso)->text().toInt();
			PN = Tab->item(n, CPN)->text().toInt();
			N=0;
		}
		S[N] = Tab->item(n, CSNR)->text().toDouble();
		v[N++] = n;
	}
	QMessageBox::information(this, "MolSpektAnalysis", 
			"No similar progression found!", QMessageBox::Ok);
}

void LineTable::setvs()
{
	QDialog SDialog(this);
	SDialog.setWindowTitle("Set vs");
	SDialog.setMaximumSize(200, 130);
	SDialog.setMinimumSize(200, 130);
	QLabel L1("Please specify the v' value to\nbe set for the selected lines.", &SDialog);
	L1.setGeometry(10, 10, 180, 40);
	QLabel L2("New v':", &SDialog);
	L2.setGeometry(10, 60, 65, 20);
	QLineEdit NE("  -1", &SDialog);
	NE.setGeometry(100, 60, 90, 20);
	NE.setAlignment(Qt::AlignRight);
	QPushButton OK("OK", &SDialog), Cancel("Cancel", &SDialog);
	OK.setGeometry(10, 100, 70, 20);
	Cancel.setGeometry(120, 100, 70, 20);
	connect(&OK, SIGNAL(clicked()), &SDialog, SLOT(accept()));
	connect(&Cancel, SIGNAL(clicked()), &SDialog, SLOT(reject()));
	if (SDialog.exec() == QDialog::Rejected) return;
	QString vs = NE.text();
	int i, m=0, n = Tab->rowCount(), r, rc;
	QList<QTableWidgetSelectionRange> SR = Tab->selectedRanges();
	Tab->blockSignals(true);
	for (r=0, rc = SR.size(); r < rc; r++)
	{
		m = SR[r].topRow();
		n = SR[r].bottomRow();
		for (i=m; i<=n; i++) Tab->item(i, Cvs)->setText(vs);
	}
	Tab->blockSignals(false);
	Changed();
}

void LineTable::setFC(QString FC)
{
	Tab->blockSignals(true);
	QList<QTableWidgetSelectionRange> SR = Tab->selectedRanges();
	int r, n;
	for (n=0; n < SR.count(); n++) for (r = SR[n].topRow(); r <= SR[n].bottomRow(); r++)
		Tab->item(r, CF)->setText(FC);
	Tab->blockSignals(false);
	Changed();
}

void LineTable::sortbyvs()
{
	int *SArray = heapSort(sortByvs);
	sortTab(SArray);
}

void LineTable::SetPN()
{
	int n, aJ = -1, nJ, ndJ, adJ = -1, aI=-1, nI, av=-1, nv, aF = -1, nF, N = Tab->rowCount(), i, j;
	int *SA = new int[N];
	double nWn, aWn = 0.0;
	QString aFi, nFi;
	bool UTA = ShowUpTerm();
	int *tSO = heapSort(sortForSPN);
	for (n=0; n<N; n++) SA[tSO[n]] = n;
	delete[] tSO;
	MaxPN = 0;
	Tab->blockSignals(true);
	for (n=0; n < N; n++)
	{
		nJ = Tab->item(SA[n], CJs)->text().toInt();
		ndJ = nJ - Tab->item(SA[n], CJss)->text().toInt();
		ndJ *= ndJ;
		nI = Tab->item(SA[n], CIso)->text().toInt();
		nv = Tab->item(SA[n], Cvs)->text().toInt();
		nWn = (UTA ? Tab->item(SA[n], CEUp)->text().toDouble() : 0.0);
		i = (nFi = Tab->item(SA[n], CFile)->text()).lastIndexOf(QRegExp("[\\/]")) + 1;
		j = nFi.indexOf('.', i);
		nFi = nFi.mid(i, (j>=0 ? j : nFi.length()) - i);
		nF = Tab->item(SA[n], CF)->text().toInt();
		if (nFi.isEmpty()) Tab->item(SA[n], CFile)->setText(nFi = aFi);
		if (nJ != aJ || ndJ != adJ || nv != av || fabs(nI - aI) > 9.0 || nF != aF || nFi != aFi 
			|| fabs(nWn - aWn) > 1e2 * Tab->item(SA[n], Cerr)->text().toDouble())
		{
			aJ = nJ;
			av = nv;
			aI = nI;
			aF = nF;
			aFi = nFi;
			aWn = nWn;
			adJ = ndJ;
			MaxPN++;
		}
		Tab->item(SA[n], CPN)->setText(QString::number(MaxPN));
	}
	Tab->blockSignals(false);
	delete[] SA;
	Changed();
}

void LineTable::getKnownLevels(int NI, int &mvs, int &mvss, int &mJs, int &mJss, bool ***&uL, 
							   bool ***&lL)
{
	//printf("LineTab::getKnownLevels\n");
	int n, m, i, o, I, NR = Tab->rowCount();
	for (n=0, mvs = 0, mvss = 0, mJs = 0, mJss = 0; n < NR; n++)
	{
		if ((m = Tab->item(n, Cvs)->text().toInt()) > mvs) mvs = m;
		if ((m = Tab->item(n, Cvss)->text().toInt()) > mvss) mvss = m;
		if ((m = Tab->item(n, CJs)->text().toInt()) > mJs) mJs = m;
		if ((m = Tab->item(n, CJss)->text().toInt()) > mJss) mJss = m;
	}
	uL = CreateBool(NI, mvs + 1, mJs + 1);
	lL = CreateBool(NI, mvss + 1, mJss + 1);
	for (i=0; i < NI; i++) 
	{
		for (n=0; n <= mvs; n++) for (m=0; m <= mJs; m++) uL[i][n][m] = false;
		for (n=0; n <= mvss; n++) for (m=0; m <= mJss; m++) lL[i][n][m] = false;
	}
	for (i=0; i < NR; i++)
	{
		I = (o = Tab->item(i, CIso)->text().toInt()) / 10;
		if (o == 10 * I && I > 0 && I <= NI)
		{
			I--; 
			if ((m = Tab->item(i, Cvs)->text().toInt()) >= 0 
						  && (n = Tab->item(i, CJs)->text().toInt()) >= 0)
				uL[I][m][n] = true;
			if ((m = Tab->item(i, Cvss)->text().toInt()) >= 0 
						  && (n = Tab->item(i, CJss)->text().toInt()) >= 0)
				lL[I][m][n] = true;
		}
	}
	//printf("NR=%d, mJs=%d, mJss=%d\n", NR, mJs, mJss);
}

void LineTable::findErrors()
{
	int n, N = Tab->columnCount();
	for (n=0; (n<N ? Tab->item(n, Cerr)->text().toDouble() != 0.0 : false); n++) ;
	if (n<N) Tab->selectRow(n);
}
