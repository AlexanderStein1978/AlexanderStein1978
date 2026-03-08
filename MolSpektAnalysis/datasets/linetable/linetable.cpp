//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
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
#include "linetablecore.h"
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
#include <QHeaderView>


LineTable::LineTable(MainWindow *MW, Molecule *M, Transition *T) : TableViewWindow(new LineTableCore, LineTab, MW, M)
{
	reinterpret_cast<LineTableCore*>(mCore)->setLineTable(this);
	setFilter("Line table (*.lines)");
	setFileExt(".lines");
	table->setModel(mCore);
	table->verticalHeader()->setVisible(true);
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
	table->setColumnWidth(LineTableCore::Cvs, 50);
	table->setColumnWidth(LineTableCore::CJs, 50);
	table->setColumnWidth(LineTableCore::Cvss, 50);
	table->setColumnWidth(LineTableCore::CJss, 50);
	table->setColumnWidth(LineTableCore::CF, 30);
	table->setColumnWidth(LineTableCore::CWN, 150);
    table->setColumnWidth(LineTableCore::CFile, 250);
	resize(600, 600);
	//connect(Tab, SIGNAL(itemChanged(QTableWidgetItem*)), this, SLOT(Changed()));
	connect(table, SIGNAL(SelChanged()), this, SLOT(TabSelChanged()));
	Saved();	
}

LineTable::~LineTable()
{
	//printf("LineTable::~LineTable\n");
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
	if (Mol != 0) Iso = Mol->getIso();
	else Iso = 0;
	mCore->setMolecule(Mol);
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
	LineTableCore* lineTableCore = reinterpret_cast<LineTableCore*>(mCore);
	int i, j, nr = lineTableCore->rowCount(), s, *n, N, iB = 0, **p, *lines, I, AM;
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
		T = lineTableCore->getSourceFile(lines[i]);
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
		T = lineTableCore->getSourceFile(lines[i]);
		if (T != "" || !L) 
		{
			for (s=0, L = false; (s<N ? T.indexOf(SpektFile[s]) == -1 : false); s++) ;
			if (s<N) L = true;
		}
		if (L)
		{
			//printf("i=%d, n[0]=%d\n", lines[i], n[0]);
			F[s][n[s]] = lineTableCore->getWaveNumber(lines[i]);
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
					I = (lineTableCore->getIso(p[s][i]) - 1)/10;
					if (I < 0 || I >= Iso->numIso) I = 0;
					if (marker[j].IsoName == Iso->texName[I]) marker[j].DisplayData = false;
				}
				else
				{
					//printf("i=%d, j=%d\n", i, j);
					//printf("p[%d][%d]=%d, Iso->numIso=%d\n", s, i, p[s][i], Iso->numIso);
					marker[j].vs = lineTableCore->get_vs(p[s][i]);
					marker[j].Js = lineTableCore->getJs(p[s][i]);
					marker[j].vss = lineTableCore->get_vss(p[s][i]);
					marker[j].Jss = lineTableCore->getJss(p[s][i]);
					marker[j].Iso = (lineTableCore->getIso(p[s][i]) - 1)/10;
					marker[j].DD = lineTableCore->getObsCalc(p[s][i]);
					marker[j].FC = lineTableCore->getFineStructureQN(p[s][i]);
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
    return mCore->getNSources();
}

void LineTable::getLines(int **Zuordnung, double *Energien, double *Unc)
{
	//printf("LineTable::getLines()\n");
	int N;
	LineTableCore* ltc = reinterpret_cast<LineTableCore*>(mCore);
    for (N=0; N < mCore->rowCount(); ++N)
    {
		Zuordnung[N][0] = (ltc->getIso(N) - 1) / 10;
		Zuordnung[N][1] = ltc->get_vs(N);
		Zuordnung[N][2] = ltc->getJs(N);
		Zuordnung[N][3] = ltc->get_vss(N);
		Zuordnung[N][4] = ltc->getJss(N);
		Zuordnung[N][5] = ltc->getFineStructureQN(N);
		Energien[N] = ltc->getWaveNumber(N);
		Unc[N] = ltc->getUncertainty(N);
    }
	//printf("Ende von getLines\n");
}

/*void LineTable::getLines(const QString &Filename, double **Lines, int *numLines)
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
}*/

void LineTable::getLines(TableLine*& L, int& N)
{
	int r;
	LineTableCore* ltc = reinterpret_cast<LineTableCore*>(mCore);
	if ((N = ltc->rowCount()) > 0) L = new TableLine[N];
	for (r=0; r<N; ++r)
	{
		L[r].dev = ltc->getObsCalc(r);
		L[r].err = ltc->getUncertainty(r);
		L[r].DevR = L[r].dev / L[r].err;
		L[r].FC = ltc->getFineStructureQN(r);
		L[r].File = ltc->getSourceFile(r);
		L[r].Iso = (ltc->getIso(r) - 1) / 10;
		L[r].isTE = false;
		L[r].Js = ltc->getJs(r);
		L[r].Jss = ltc->getJss(r);
		L[r].LTab = this;
		L[r].PN = ltc->getProgression(r);
		L[r].Row = r;
		L[r].SourceN = 0;
		L[r].vs = ltc->get_vs(r);
		L[r].vss = ltc->get_vss(r);
		L[r].WN = ltc->getWaveNumber(r);
	}
}

void LineTable::getSortedLines(TableLine*& L, int& N, int SortOrder)
{
	LineTableCore* ltc = reinterpret_cast<LineTableCore*>(mCore);
	if (SortOrder != 0 || (N = mCore->rowCount()) == 0)
	{
		N=0;
		L=0;
		return;
	}
	int r, *SO = heapSort(sortByFrequency);
	L = new TableLine[N];
	for (r=0; r<N; r++)
	{
		L[SO[r]].dev = ltc->getObsCalc(r);
		L[SO[r]].err = ltc->getUncertainty(r);
		L[SO[r]].DevR = L[r].dev / L[r].err;
		L[SO[r]].FC = ltc->getFineStructureQN(r);
		L[SO[r]].File = ltc->getSourceFile(r);
		L[SO[r]].Iso = (ltc->getIso(r) - 1) / 10;
		L[SO[r]].isTE = false;
		L[SO[r]].Js = ltc->getJs(r);
		L[SO[r]].Jss = ltc->getJss(r);
		L[SO[r]].LTab = this;
		L[SO[r]].PN = ltc->getProgression(r);
		L[SO[r]].Row = r;
		L[SO[r]].SourceN = 0;
		L[SO[r]].vs = ltc->get_vs(r);
		L[SO[r]].vss = ltc->get_vss(r);
		L[SO[r]].WN = ltc->getWaveNumber(r);
		L[SO[r]].SNR = ltc->getSNR(r);
	}
	delete[] SO;
}

int LineTable::getNgTE(int *mv, int mJ)
{
	LineTableCore* ltc = reinterpret_cast<LineTableCore*>(mCore);
	int n, I, N, P, p, NR = ltc->rowCount(), v, J;
	if (NSO != NR) ShowUpTerm();
	for (n=N=0, P=-1; n < NR; n++) if ((p = ltc->getProgression(SO[n])) != P)
	{
		if (ltc->getUpperEnergy(SO[n]) == 0.0) continue;
		v = ltc->get_vs(SO[n]);
		J = ltc->getJs(SO[n]);
		if ((v >= 0 && J >= 0) || (mv != 0 && (J <= mJ || v > mv[J]))) continue;
		I = ltc->getIso(SO[n]);
		while (10 * (I / 10) != I && n < NR) I = ltc->getIso(SO[n++]);
		if (n == NR) break;
		P = p;
		N++;
	}
	return N;
}

void LineTable::getgoodTE(int &N, TermEnergy *&E, int *mv, int mJ)
{
	//printf("LineTable::getgoodTE\n");
	LineTableCore* ltc = reinterpret_cast<LineTableCore*>(mCore);
	int n, I, P = -1, p, NR = ltc->rowCount(), J, v;
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
	for (n=N=0; n < NR; ++n) if ((p = ltc->getProgression(SO[n])) != P)
	{
		//printf("N=%d, n=%d\n", N, n);
		if (ltc->getUpperEnergy(SO[n]) == 0.0) continue;
		v = ltc->get_vs(SO[n]);
		J = ltc->getJs(SO[n]);
		if (mv != 0 ? ((v >= 0 && J >= 0) || (J <= mJ || v > mv[J])) : J < 0) continue;
		I = ltc->getIso(SO[n]);
		while (10 * (gli[N][0] = I / 10) != I && n < NR) I = ltc->getIso(SO[n++]);
		if (n == NR) break;
		P = p;
		gli[N][1] = SO[n];
		--(gli[N++][0]);
	}
	//printf("Vor Schleife2\n");
	E = new TermEnergy[N];
	for (n=0; n<N; ++n)
	{
		E[n].v = ltc->get_vs(gli[n][1]);
		E[n].J = ltc->getJs(gli[n][1]);
		E[n].Iso = gli[n][0];
		E[n].err = 0.01;
		E[n].E = ltc->getAverageUpperEnergy(gli[n][1]);
		E[n].PN = ltc->getProgression(gli[n][1]);
		E[n].File = ltc->getSourceFile(gli[n][1]);
		E[n].ef = ltc->getJss(gli[n][1]) != E[n].J;
		E[n].dev = E[n].DevR = 0.0;
		E[n].FC = ltc->getFineStructureQN(gli[n][1]);
	}
	//printf("Ende getgoodTE");
}

int LineTable::getNgL(int *mv, int mJ)
{
	LineTableCore* ltc = reinterpret_cast<LineTableCore*>(mCore);
	int n, I, N, NR = ltc->rowCount(), J;
	for (n=N=0; n < NR; ++n)
	{
		if (mv != 0 && (J = ltc->getJss(n)) <= mJ && ltc->get_vss(n) > mv[J]) continue;
		I = ltc->getIso(n);
		if (10 * (I / 10) == I) ++N;
	}
	return N;
}

void LineTable::getgoodLines(int &N, TableLine *&L, int *mv, int mJ, bool SortFunction(const TableViewWindowCore *const, const int, const int))
{
	//printf("Beginn LineTable::getgoodLines\n");
	LineTableCore* ltc = reinterpret_cast<LineTableCore*>(mCore);
	int n, i, NR = ltc->rowCount(), J;
	//printf("NR=%d\n", NR);
    int *SO = new int[NR], *S1 = heapSort(SortFunction);
	int gli[NR][2];
	for (n=0; n < NR; n++) SO[S1[n]] = n;
	for (n = N = 0; n < NR; ++n)
	{
		if (mv != 0 && (J = ltc->getJss(SO[n])) <= mJ && ltc->get_vss(SO[n]) > mv[J]) continue;
		//printf("n=%d, NR=%d\n", n, NR);
		i = ltc->getIso(SO[n]);
		if (10 * (gli[N][0] = i / 10) == i) 
		{
			gli[N][1] = SO[n];
			--(gli[N++][0]);
		}
	}
	L = new TableLine[N];
	//printf("n=%d, N=%d\n", n, N);
	for (n=0; n<N; ++n)
	{
		L[n].FC = ltc->getFineStructureQN(gli[n][1]);
		L[n].PN = ltc->getProgression(gli[n][1]);
		L[n].vs = ltc->get_vs(gli[n][1]);
		L[n].Js = ltc->getJs(gli[n][1]);
		L[n].vss = ltc->get_vss(gli[n][1]);
		L[n].Jss = ltc->getJss(gli[n][1]);
		L[n].Iso = gli[n][0];
		L[n].WN = ltc->getWaveNumber(gli[n][1]);
		L[n].err = ltc->getUncertainty(gli[n][1]);
		L[n].File = ltc->getSourceFile(gli[n][1]);
		L[n].Row = gli[n][1];
	}
	delete[] SO;
	delete[] S1;
}

int LineTable::getMaxvs()
{
	int i, v, n=0, nr = mCore->rowCount();
	for (i=0; i<nr; ++i) if ((v = reinterpret_cast<LineTableCore*>(mCore)->get_vs(i)) > n) n = v;
	return n;
}

int LineTable::getMaxvss()
{
	int i, v, n=0, nr = mCore->rowCount();
	for (i=0; i < nr; ++i) if ((v = reinterpret_cast<LineTableCore*>(mCore)->get_vss(i)) > n) n = v;
	return n;
}

int LineTable::getMaxJs()
{
	int i, v, n=0, nr = mCore->rowCount();
	for (i=0; i < nr; ++i) if ((v = reinterpret_cast<LineTableCore*>(mCore)->getJs(i)) > n) n = v;
	return n;
}

int LineTable::getMaxJss()
{
	int i, v, n=0, nr = mCore->rowCount();
	for (i=0; i < nr; ++i) if ((v = reinterpret_cast<LineTableCore*>(mCore)->getJss(i)) > n) n = v;
	return n;
}

void LineTable::getObsIso(bool *O, int N)
{
	int I, r, nr = mCore->rowCount();
	for (r=0; r<N; ++r) O[r] = false;
	for (r=0; r < nr; ++r) if ((I = reinterpret_cast<LineTableCore*>(mCore)->getIso(r) / 10 - 1) < N && I >= 0) O[I] = true;
}

void LineTable::getProgressions(int &N, Progression *&P, bool SortFunction(const TableViewWindowCore *const, const int, const int))
{
	//printf("LineTable::getProgressions\n");
	TableLine *L;
	QString S, aS;
	int vs = -1, Js = -1, Iso = -1, n, m, NP=0, NR, PN = -1;
    getgoodLines(NR, L, 0, 0, SortFunction);
	int lpb[NR];
	for (n=0; n < NR; ++n) if (L[n].vs != vs || L[n].Js != Js || L[n].Iso != Iso || L[n].File != S || L[n].PN != PN)
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
	for (n=0; n<N; ++n)
	{
		P[n].N = lpb[n+1] - lpb[n];
		P[n].vs = L[lpb[n]].vs;
		P[n].Js = L[lpb[n]].Js;
		P[n].Iso = L[lpb[n]].Iso;
		P[n].L = new Line[P[n].N];
        P[n].PNum = L[lpb[n]].PN;
		//printf("N=%d, P[%d].N=%d\n", N, n, P[n].N);
		for (m=0; m < P[n].N; ++m)
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
	table->blockSignals(true);
	mCore->blockSignals(true);
	for (n=0; n < NLines; ++n) mCore->setUncertainty(RowNumbers[n], NewUncertainty[n]);
	mCore->blockSignals(false);
	table->blockSignals(false);
	Changed();
}

bool LineTable::readData(QString Filename)
{
    //printf("LineTable::readData\n");
	int n, s, i, N = 0, d=0;
	bool FCA = false, Success = true;
    QFile Datei(Filename);
    if (!read(&Datei)) return false;
    QTextStream S1(&Datei);
    QString B, Buffer, B2 = S1.readLine();
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
    LineTableCore* ltc = reinterpret_cast<LineTableCore*>(mCore);
	ltc->setSowNotColumnCount(true);
	ltc->readData(N, Buffer, FCA, MaxPN, d);
//	ltc->removeNotNeededColumns();
	ltc->setSowNotColumnCount(false);
	Saved();
	emit DataChanged();
	return Success;
}

bool LineTable::writeData(QString Filename)
{
    int i, j;
    QString Buffer, B;
	QPixmap P;
	QFile Datei(Filename);
	write(&Datei);
    QTextStream S(&Datei);
	LineTableCore* ltc = reinterpret_cast<LineTableCore*>(mCore);
	S << "Source: " << getSource() << "\n";
	S << "Name: " << getName() << "\n";
	S << "   PN  vo  Jo  vu  Ju  FC           Eout       err  Iso  File  SNR  deviation  comment\n";
    table->blockSignals(true);
	ltc->blockSignals(true);
	ltc->writeData(S);
    Datei.close();
	ltc->blockSignals(false);
	table->blockSignals(false);
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
	if (A1 == A2) for (n=i=0; n < nIsoA; ++n) for (m=n; m < nIsoA; ++m)
	{
		IsoTA[i] = n;
		IsoTB[i++] = m;
	}
	else for (n=i=0; n < nIsoA; ++n) for (m=0; m < nIsoB; ++m)
	{
		IsoTA[i] = n;
		IsoTB[i++] = m;
	}
	for (n=0; n < nIsoA; ++n) IsoCA[n] = 0;
	for (n=0; n < nIsoB; ++n) IsoCB[n] = 0;
	TermEnergy *D = new TermEnergy[ND], *FD = 0, TBuff;
	for (n=m=0; n < ND; ++n)
	{
		if (Data[n][4] == "nan") continue;
		D[m].Iso = Data[n][0].toInt() / 10 - 1;
		if ((D[m].v = Data[n][1].toInt()) > Maxv) Maxv = D[n].v;
		if ((D[m].J = Data[n][2].toInt()) > MaxJ) MaxJ = D[n].J;
		D[m].ef = (Data[n][3].indexOf('e') != -1);
		D[m].E = Data[n][4].toDouble();
		if (D[m].E == 0.0) continue;
		D[m].err = 0.01;
		OUnc[m] = Data[n][5].toDouble();
		++(IsoCA[IsoTA[D[m].Iso]]);
		++(IsoCB[IsoTB[D[m++].Iso]]);
	}
	ND = m;
	for (n=0; n < ND; ++n)
	{
		for (m=n+1; D[m].Iso == D[n].Iso && D[n].v == D[m].v && D[n].J == D[m].J && D[n].ef == D[m].ef; ++m) ;
		if (m==n+1) continue;
		for (i=n, WSum = DS = 0.0; i<m; ++i)
		{
			WSum += (W = ((W = OUnc[m]) != 0.0 ? 1.0 / W : 1.0));
			DS += W * D[i].E;
		}
		D[n].E = DS / WSum;
		for (i=n+1; i<m; ++i) D[i].E = 0.0;
	}
	Destroy(Data, ND);
	for (n=0; (n < ND && D[n].E != 0.0); ++n) ;
	for (m=n+1; m < ND; ++m) if (D[m].E != 0.0) D[n++] = D[m];
	ND = n;
	uIA[0] = IsoTA[Iso->refIso];
	uIB[0] = IsoTB[Iso->refIso];
	for (i=m=0; i < nIsoA; ++i) if (IsoCA[i] > m && i != uIA[0])
	{
		m = IsoCA[i];
		uIA[1] = i;
	}
	if (m>0) 
	{
		for (i=m=0; i < nIsoA; ++i) if (IsoCA[i] > m && i != uIA[0] && i != uIA[1])
		{
			m = IsoCA[i];
			uIA[2] = i;
		}
		nIA = (m>0 ? 3 : 2);
	}
	else nIA = 1;
	for (i=m=0; i < nIsoB; ++i) if (IsoCB[i] > m && i != uIB[0])
	{
		m = IsoCB[i];
		uIB[1] = i;
	}
	if (m>0)
	{
		for (i=m=0; i < nIsoB; ++i) if (IsoCB[i] > m && i != uIB[0] && i != uIB[1])
		{
			m = IsoCB[i];
			uIB[2] = i;
		}
		nIB = (m>0 ? 3 : 2);
	}
	else nIB = 1;
	for (n=0; n < nIA; ++n)
		IsoSA[n] = ("    " + QString::number(A1->getnNuc(uIA[n]))).right(5)
				 + ("          " 
				   + QString::number(A1->getIsoMass(uIA[n]), 'f', 7)).right(13)
				 + "           0.0";
	for (n=0; n < nIB; ++n)
		IsoSB[n] = ("    " + QString::number(A2->getnNuc(uIB[n]))).right(5)
				 + ("          " 
				   + QString::number(A2->getIsoMass(uIB[n]), 'f', 7)).right(13)
				 + "           0.0";
	QFile File(Filename);
	QString Buffer, L1, LP1, LP2, IsoStr[nIso];
	QuestionBox *QB;
	for (n=0; n < nIso; ++n)
	{
		for (m=0; m < nIA; ++m) if (Iso->mNumIso1[n] == IsoSA[m].left(5).toInt()) IsoStr[n] += IsoSA[m].left(5);
		for (m=0; m < nIB; ++m) if (Iso->mNumIso2[n] == IsoSB[m].left(5).toInt()) IsoStr[n] += IsoSB[m].left(5);
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
						delete[] D;
						return true;
						break;
				}
			}
			LP1 = RS.readLine();
			LP2 = RS.readLine();
			FD = new TermEnergy[NFD];
			for (n=m=0; n < NFD; ++n)
			{
				Buffer = RS.readLine();
				iA = Buffer.left(5).toInt();
				iB = Buffer.mid(5, 5).toInt();
				FD[m].v = Buffer.mid(10, 5).toInt();
				FD[m].J = Buffer.mid(15, 5).toInt();
				if (FD[m].J != Buffer.mid(25, 5).toInt()) FD[m].ef = true;
				else FD[m].ef = false;
				Buffer = RS.readLine();
				for (i = 0; i < nIso && (iA != Iso->mNumIso1[i] || iB != Iso->mNumIso2[i]); ++i) ;
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
					++m;
				}
			}
			File.close();
		}
		else if (R==2)
		{
			delete[] D;
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
	for (n=m=0; n < ND || (NFD > 0 && FD[0].Iso != -1); ++m)
	{
		if (n < ND && NFD > 0 ? D[n].Iso == FD[0].Iso && D[n].v == FD[0].v && D[n].J == FD[0].J && D[n].ef == FD[0].ef
			  : false)
		{
			TBuff = D[n++];
			TBuff.err = FD[0].err;
			FD[0].Iso = -2;
		}
		else if (n < ND && ((NFD > 0 && FD[0].Iso != -1) || isnSPG(D[n], FD[0]))) TBuff = D[n++];
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
		if (NFD > 0 && FD[0].Iso == -2)
		{
			for (p1 = 0, p2 = 2; p2 <= NFD && FD[p1].Iso != -1; p2 = 2 * ((p1 = p2) + 1))
			{
				if (p2 < NFD && isnSPG(FD[p2-1], FD[p2])) p2--;
				FD[p1] = FD[p2];
			}
			FD[p1].Iso = -1;
		}
		if (TBuff.Iso == -1 || IsoStr[TBuff.Iso].length() < 10) --m;
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
	for (n=0; n<m; ++n) WS << LL[n];
	File.close();
	delete[] LL;
	delete[] FD;
	delete[] D;
	delete[] OUnc;
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
	int i, j, N;
	N = mCore->rowCount();
	if (FileName.isEmpty()) 
		FileName = QFileDialog::getSaveFileName(this, "Write TF input file", 
			(MW != nullptr ? MW->getDir(getType()) : ""), "All files (*.*)");
	if (FileName.isEmpty()) return;
	QFile Datei(FileName);
	Datei.open(QIODevice::WriteOnly | QIODevice::Append);
	QTextStream S(&Datei);
	QString Buffer, IBuff, F, lF;
	vsOListElement *CvsOElement = vsOList;
	if (CvsOElement == nullptr)
	{
		CvsOElement = vsOList = new vsOListElement;
		CvsOElement->Iso = CvsOElement->Js = CvsOElement->vs = -1;
		CvsOElement->next = nullptr;
	}
	reinterpret_cast<LineTableCore*>(mCore)->WriteTFGS(S, vsOList, vsO);
	Datei.close();
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
		for (i=0; i <= mvsldc; ++i) sLD[i][0] = sLD[i][1] = 0.0;
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
	for (n=0, N=0; (N<j && Data[N][0] != ""); ++N) if ((Data[N][6].toInt() >= MinQL && Data[N][3] == "f") || (Data[N][7].toInt() >= MinDPR && Data[N][3] == "e"))
	{
		l = (k = Data[N][0].toInt()) / 10;
		if (l > 0 && k == 10 * l) ++n;
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
		for (i=0; i < ZINI1; ++i)
			ZII1[i] = ("    " + QString::number(atom->getnNuc(i))).right(5) 
					+ "   " + QString::number(atom->getIsoMass(i), 'f', 8);
		for (i = RI1, Buffer = ZII1[RI1]; i>0; i--) ZII1[i] = ZII1[i-1];
		ZII1[0] = Buffer;
		atom = molecule->getAtom2();
		ZINI2 = atom->getnIso();
		ZII2 = new QString[ZINI2];
		for (i=0; i < ZINI2; ++i)
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
			for (i = ZINC = 0; i <= KMax; ++i)
			{
				ZINC += LMax[i] + 1;
				if (Kor != 0) for (j=0; j <= LMax[i]; j++) if (Kor[i][j] != 0.0) ++ZINC;
			}
			ZIC = new QString[ZINC];
			for (i=j=0; j <= KMax; j++) for (k=0; k <= LMax[j]; ++k)
			{
				ZIC[++i] = ("    " + QString::number(k)).right(5)
						 + ("    " + QString::number(j)).right(5)
						 + (Par[j][k] != 0.0 ? "    1    10.0000000000000000E+00    0    2"
										     : "    1    10.0000000000000000E+00    1    2");
				if (Kor != 0 ? Kor[j][k] != 0.0 : false) 
					ZIC[++i] = ("    " + QString::number(k)).right(5)
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
	for (i=0; i < IsoT->numIso; ++i)
		Iso[i] = ("    " + QString::number(IsoT->mNumIso1[i])).right(5) 
			   + ("    " + QString::number(IsoT->mNumIso2[i])).right(5);
	QTextStream S(&Datei);
	//printf("Vor write\n");
	ZIL[3] = ("     " + QString::number(n)).right(5) + ZIL[3].right(ZIL[3].length() - 5);
	for (i=0; i<4; ++i) S << ZIL[i] << "\n";
	for (i=0; i < ZINC; ++i) S << ZIC[i] << "\n";
	for (i=0; i < ZINI1; ++i) S << ZII1[i] << "\n";
	for (i=0; i < ZINI2; ++i) S << ZII2[i] << "\n";
	//printf("Vor Schleife\n");
	for (i=0; i<N; ++i)
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
	delete[] Iso;
	delete[] ZIC;
	delete[] ZII1;
	delete[] ZII2;
}

void LineTable::Assignvs()
{
	//printf("LineTable::Assignvs()\n");
	int i, nr = mCore->rowCount();
	ShowUpTerm();
	//printf("Nach ShowUpTerm\n");
	if (molecule == 0 || transition == 0)
	{
		QMessageBox::information(this, tr("QT4MolSpektAn"), tr("The line table has to be assigned to a molecule and a transition first!"));
		return;
	}
	ElState *US = transition->getUpperState();
	if (US == 0)
	{
		QMessageBox::information(this, tr("QT4MolSpektAn"), tr("The transition has to be assigned to an upper electronic state first!"));
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
	double ****UD = TT->getData(), AT = Tol->text().toDouble();
	int NumC = TT->getNumComp(), mvs = TT->getMaxv();
	//printf("mvs=%d, NumC=%d\n", mvs, NumC);
	table->blockSignals(true);
	mCore->blockSignals(true);
	reinterpret_cast<LineTableCore*>(mCore)->Assign_vs(UD, AT, NumC, mvs, SO);
	mCore->blockSignals(false);
	table->blockSignals(false);
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
	int cTX[XNC], cTE[XNC], cTL[XNC];
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
	int XNJ = XTT->getMaxJ() + 1, ENJ = ETT->getMaxJ() + 1;
	double ****XData = XTT->getData(), ****EData = ETT->getData();
	int *LSO = heapSort(sortIJvP), NR = mCore->rowCount();
	int *XIT = XTT->getIsoT(), *EIT = ETT->getIsoT();
	int *LO = new int[NR];
	for (n=0; n < NR; n++) LO[LSO[n]] = n;
	delete[] LSO;
	table->blockSignals(true);
	mCore->blockSignals(true);
	reinterpret_cast<LineTableCore*>(mCore)->AssignFC(LO, XIT, EIT, ENv, ENJ, EData, XNv, XNJ, XData, XNC, cTX, cTE, cTL);
	mCore->blockSignals(false);
	table->blockSignals(false);
	delete[] LO;
	delete[] XCT;
	delete[] XIT;
	delete[] EIT;
	delete[] ECT;
	Changed();
}

void LineTable::MarkLines(int* Iso, int* vs, int* Js, double* WN, int N)
{
	LineTableCore* ltb = reinterpret_cast<LineTableCore*>(mCore);
	int n, m, lJs, lIso = -1, lvs = -2, r, r1, R = ltb->rowCount(), MC = ltb->columnCount() - 1, *rows = new int[N];
	double lWN = 0.0;
	for (r=0; r<R; ++r)
	{
		lJs = ltb->getJs(r);
		if (lvs != -2)
		{
			lvs = -2;
			if (lIso != -1)
			{
				lIso = -1;
				lWN = 0.0;
			}
		}
		for (n=0; n<N; ++n) if (lJs == Js[n] && ((lvs == -2 ? (lvs = ltb->get_vs(r)) : lvs) == vs[n]) && ((lIso == -1 ? (lIso = (ltb->getIso(r) - 1) / 10) : lIso) == Iso[n])
			&& (fabs((lWN == 0.0 ? (lWN = ltb->getWaveNumber(r)) : lWN) - WN[n]) < 1e-4)) break;
		if (n<N) break;
	}
	if (r==R) return;
	DeselectEverything();
	QModelIndex topIndex = mCore->getIndex(r, 0);
	table->scrollTo(topIndex, QAbstractItemView::PositionAtTop);
	for (r1 = ++r, m=-1; r<R; ++r)
	{
		lJs = ltb->getJs(r);
		if (lvs != -2)
		{
			lvs = -2;
			if (lIso != -1)
			{
				lIso = -1;
				lWN = 0.0;
			}
		}
		for (n=0; n<N; ++n) if (lJs == Js[n] && ((lvs == -2 ? (lvs = ltb->get_vs(r)) : lvs) == vs[n]) && ((lIso == -1 ? (lIso = (ltb->getIso(r) - 1) / 10) : lIso) == Iso[n])
			&& fabs((lWN == 0.0 ? (lWN = ltb->getWaveNumber(r)) : lWN) - WN[n]) < 1e-4) break;
		if (n<N && r1 == -1) r1 = r;
		else if (n==N && r1 != -1) 
		{
			rows[++m] = r1;
			r1 = -1;
		}
	}
	TableViewWindow::MarkLines(rows, ++m);
	delete[] rows;
	if (!isVisible()) show();
	activateWindow();
	table->setFocus();
}

void LineTable::MarkSelected()
{
    //printf("LineTable::MarkSelected\n");
    if (MW == 0) return;
	int k=0, SR, r, s, AnzahlMarker, n, N, *rows;
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
	LineTableCore* ltb = reinterpret_cast<LineTableCore*>(mCore);
	NR = ltb->rowCount();
	table->getSelectedRows(rows, N);
	if (nullptr == rows) return;
	for (r = 0; r<=N; ++r)
	{
		for (mSpectrum = ltb->getSourceFile(rows[r]), n = rows[r] - 1; mSpectrum.isEmpty() && n>=0; --n) mSpectrum = ltb->getSourceFile(n);
		if (!mSpectrum.isEmpty())
		{
			for (n=0; (n < SL.count() ? SL[n] != mSpectrum : false); ++n) ;
			if (n == SL.count())
			{
				SL << mSpectrum;
				iL << 1;
			}
			else ++(iL[n]);
		}
	}
	for (n=r=0; n < iL.count(); ++n) if (iL[n] > r) r = iL[s=n];
	mSpectrum = SL[s];
	mIso = ltb->getIso(SR = rows[0]);
	mJs = ltb->getJs(SR);
	mvs = ltb->get_vs(SR);
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
    for (r=0; r<N; ++r)
	{
		if (ltb->getSourceFile(rows[r]) != mSpectrum) continue;
		WN = ltb->getWaveNumber(rows[r]);
		while ((k >= 0 ? marker[k].Line[0] > WN : false)) --k;
		++k;
		while ((k < AnzahlMarker ? marker[k].Line[0] < WN || k==0 : false)) ++k;
		if ((k < AnzahlMarker ? marker[k].Line[0] - WN > WN - marker[k-1].Line[0] : true)) --k;
		//printf("marker[%d].Line[0]=%f, WN=%f\n", k, marker[k].Line[0], WN);
		if (fabs(marker[k].Line[0] - WN) < 0.001)
		{
		    //printf("Markiere marker[%d]\n", k);
		    marker[k].vs = ltb->get_vs(rows[r]);
		    marker[k].Js = ltb->getJs(rows[r]);
		    marker[k].vss = ltb->get_vss(rows[r]);
		    marker[k].Jss = ltb->getJss(rows[r]);
			marker[k].Iso = (ltb->getIso(rows[r]) - 1) / 10;
			marker[k].FC = ltb->getFineStructureQN(rows[r]);
		    marker[k].DD = ltb->getObsCalc(rows[r]);
			marker[k].Mol = molecule;
		    marker[k].IsoName = (Iso != 0 ? Iso->texName[marker[k].Iso] : "");
			marker[k].lState = lSN;
			marker[k].uState = uSN;
			marker[k].LState = lState;
			marker[k].UState = uState;
			marker[k].DisplayData = true;
		    marker[k].Marked = true;
			marker[k].uncertainty = ltb->getUncertainty(rows[r]);
			Comment = ltb->getComment(rows[r]);
			marker[k].satellite = Comment.indexOf("satellite", 0, Qt::CaseInsensitive) >= 0;
			marker[k].overlap = Comment.indexOf("overlap", 0, Qt::CaseInsensitive) >= 0;
		}
	}
	spektrum->FocusSelected();
    activateWindow();
    setFocus();
  	spektrum->activateWindow();
	spektrum->setFocus();
	delete[] rows;
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
	LineTableCore* ltb = reinterpret_cast<LineTableCore*>(mCore);
	if (ltb->columnCount() <= LineTableCore::CFCF) ltb->setColumnCount(LineTableCore::CFCF + 1);
	int NR = ltb->rowCount(), r, n, s, PN, I, Js, vs, MJ, Mv, v;
	TermTable *lTerm = lS->getTermTable(), *uTerm = uS->getTermTable();
	int uNv = uTerm->getMaxv() + 1, uNJ = uTerm->getMaxJ() + 1, uNIso = uTerm->getNumIso();
	int lNv = lTerm->getMaxv() + 1, lNJ = lTerm->getMaxJ() + 1, lNIso = lTerm->getNumIso();
	double ****uData = uTerm->getData(), ****lData = lTerm->getData();
	int lJ[2 * lNv];
	double lE[2 * lNv], FCF[2 * lNv], IS, FS;
	if (uNIso < lNIso) lNIso = uNIso;
	for (r=0; r < NR; r=s)
	{
		I = ltb->getIso(r) / 10 - 1;
		//printf("r=%d, s=%d, I=%d\n", r, s, I);
		for (IS = 0.0, n = Mv = MJ = 0, PN = ltb->getProgression(r), s=r; (s < NR ? ltb->getProgression(s) == PN : false); ++s, ++n)
		{
			IS += ltb->getSNR(s);
			if ((lJ[n] = ltb->getJss(s)) > MJ) MJ = lJ[n];
			if ((v = ltb->get_vss(s)) > Mv) Mv = v;
			if (MJ < uNJ && Mv < uNv && I < lNIso) lE[n] = lData[0][I][v][lJ[n]];
		}
		//printf("PN=%d, MJ=%d, lNJ=%d, Js=%d, uNJ=%d, I=%d, lNIso=%d\n", PN, MJ, lNJ, Js, uNJ, I, lNIso);
		if (MJ >= lNJ || Mv >= lNv  || ltb->getJs(r) >= uNJ || (vs = ltb->get_vs(r)) >= uNv || I >= lNIso) continue;
		//printf("r=%d, s=%d, Js=%d, n=%d\n", r, s, Js, n);
        lPot->getFastFCF(I, uPot, Js, 0, uData[0][I][vs][Js], n, lE, lJ, FCF, NumWFPoints);
		//printf("Nach getFastFCF, n=%d\n", n);
		for (FS = 0.0, v=0; v<n; ++v)
		{
			//printf("v=%d, n=%d\n", v, n);
			//printf("FS=%g, FCF[%d]=%g\n", FS, v, FCF[v]);
			FS += FCF[v];
		}
		//printf("FS=%g\n", FS);
		for (IS /= FS, v=0, n=r; n<s; ++v, ++n) ltb->setFCF(n, FCF[v] * IS);
	}
}

void LineTable::ShowGSDeviations()
{
	ElState *S;
	if (transition == 0)
	{
		QMessageBox::information(this, "MolSpektAnalysis", tr("To use this function the linetable as to be assigned to a transition first!"));
		return;
	}
	if ((S = transition->getLowerState()) == 0)
	{
		QMessageBox::information(this, "MolSpektAnalysis", tr("To use this function the linetable as to be assigned to an electronic ground state first!"));
		return;
	}
	TermTable *TT;
	if ((TT = S->getTermTable()) == 0)
	{
		QMessageBox::information(this, "MolSpektAnalysis", tr("To use this function for the ground state term energies have to be available first!"));
		return;
	}
	SetPN();
	LineTableCore* ltb = reinterpret_cast<LineTableCore*>(mCore);
	int n, m, pn = 0, apn = -1, N = ltb->rowCount(), nv = TT->getMaxv() + 1, nI = TT->getNumIso();
	int nJ = TT->getMaxJ(), v, J, I, ldr = 0;
	double Dev, US = 0.0, SS = 0.0, Sig, *E = new double[N], *UT = new double[N], ****Data = TT->getData();
	double MD = 0.0, AD;
	bool l = false;
	int *S1 = heapSort(isnSPG), *SA = new int[N];
	for (n=0; n<N; ++n) SA[S1[n]] = n;
	for (n=m=0; n<N; ++m)
	{
		if (m<N || (pn = ltb->getProgression(SA[m]) != apn))
		{
			for (US /= SS; n<m; ++n)
			{
				if (UT[n] == 0.0) continue;
				Dev = UT[n] - US;
				ltb->setObsCalc(SA[n], Dev);
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
			if (n==N || (l && ltb->getComment(SA[ldr]).indexOf("pd") == -1)) break;
			apn = pn;
			SS = US = 0.0;
		}
		if ((I = ltb->getIso(SA[m]) / 10 - 1) >= nI || (J = ltb->getJss(SA[m])) >= nJ || (v = ltb->get_vss(SA[m])) >= nv)
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
		E[m] = ltb->getUncertainty(SA[m]);
		UT[m] = Data[0][I][v][J] + ltb->getWaveNumber(SA[m]);
		Sig = 1.0 / (E[m] * E[m]);
		SS += Sig;
		US += Sig * UT[m];
	}
	if (l) TableViewWindow::MarkLines(&ldr, 1);
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
	int NR = mCore->rowCount(), NC = mCore->columnCount();
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
	LineTableCore* ltb = reinterpret_cast<LineTableCore*>(mCore);
	table->blockSignals(true);
	mCore->blockSignals(false);
	for (r=0; r < ltb->rowCount(); ++PN)
	{
		for (s=r; s < NR && mCore->getProgression(s) == PN; ++s) R[s-r] = true;
		printf("PN=%d, r=%d, s=%d, l=%d, NR=%d\n", PN, r, s, l, NR);
		for (Spekt = "", n=r; mSpectrum.isEmpty() && n<s; ++n) Spekt = ltb->getSourceFile(r);
		if (!spektrum->readData(Spekt))
		{
			if (l>r) for (n=r; n<s; n++, l++) ltb->deleteRow(n);
			continue;
		}
		spektrum->editFind();
		spektrum->GetMarker(AnzahlMarker, marker);
		laser = spektrum->GetLaserLine();
		sat = ltb->getComment(r).indexOf("satellite") > -1;
		printf("AnzahlMarker=%d\n", AnzahlMarker);
		for (Success = false, k=0; !Success; )
		{
			spektrum->ClearMarked();
			printf("r=%d, s=%d\n", r, s);
			for (n = r; n < s; ++n) if (R[n-r])
			{
				WN = ltb->getWaveNumber(n);
				while (k >= 0 && marker[k].Line[0] > WN) --k;
				++k;
				while (k < AnzahlMarker && (marker[k].Line[0] < WN || k==0)) ++k;
				if (k < AnzahlMarker || marker[k].Line[0] - WN > WN - marker[k-1].Line[0]) --k;
				//printf("marker[%d].Line[0]=%f, WN=%f\n", k, marker[k].Line[0], WN);
				if (fabs(marker[k].Line[0] - WN) < 0.001)
				{
		    		//printf("Markiere marker[%d]\n", k);
		    		marker[k].vs = ltb->get_vs(n);
		    		marker[k].Js = ltb->getJs(n);
		    		marker[k].vss = ltb->get_vss(n);
		    		marker[k].Jss = ltb->getJss(n);
					marker[k].Iso = ltb->getIso(n) / 10 - 1;
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
			for (n=0, MDev = 0.0; n < AnzahlMarker; ++n) if (marker[n].Marked && (D = fabs(marker[n].DD) > MDev))
			{
				MDev = D;
				m=n;
			}
			printf("m=%d, MDev=%e, MaxDev=%e\n", m, MDev, CMaxSearchDev);
			if (MDev > CMaxSearchDev) 
			{
				for (n=r; P[n] != m; ++n) ;
				printf("n=%d, m=%d\n", n, m);
				R[n-r] = false;
			}
			else Success = true;
		}
		r=s;
		for (n=m=0; n < AnzahlMarker; ++n) if (marker[n].Marked)
		{
			++m;
			k=n;
		}
		printf("m=%d, l=%d, s=%d\n", m, l, s);
		if (m<4 || (uState != 0 && marker[k].vs == -1)) continue;
		if (l+m < s) 
		{
			for (k=0; k < AnzahlMarker; ++k) if (marker[k].Marked)
			{
				LineTableBaseData* data = new LineTableBaseData;
				data->progressionNumber = PN;
				data->vs = marker[k].vs;
				data->Js = marker[k].Js;
				data->vss = marker[k].vss;
				data->Jss = marker[k].Jss;
				data->isotope = 10 * (marker[k].Iso + 1);
				data->waveNumber = marker[k].Line[0];
				data->uncertainty = (marker[k].overlap ? 0.020 : 0.005);
    			data->file = Spekt;
				data->SNR = marker[k].SNR;
				data->obsMinusCalc = marker[k].DD;
				if (sat) Buffer = "satellite";
				else if (marker + k == laser) Buffer = "laser";
				else Buffer = ""; 
				if (marker[k].overlap) Buffer += "overlap";
				data->Comment = Buffer;
				ltb->addRow(data);
			}
			continue;
		}
	}
	mCore->blockSignals(false);
	table->blockSignals(false);
	delete[] R;
	delete[] P;
	Changed();
}

void LineTable::TakeOnChanges()
{
	//printf("mSpectrum=%s\n", mSpectrum.ascii());
	Spektrum *Spekt = (MW != 0 ? MW->getSpectrum(mSpectrum) : 0);
	if (Spekt == 0) return;
	int j, k, n = mCore->rowCount(), c, nc = mCore->columnCount(), *rows, N;
	QString B, Buffer;
	table->getSelectedRows(rows, N);
	if (nullptr == rows) return;
	table->blockSignals(true);
	mCore->blockSignals(true);
	for (k=0; k<N; ++k)	if (mCore->getSourceFile(rows[k]) == mSpectrum) mCore->deleteRow(k);
	mCore->blockSignals(false);
	table->blockSignals(false);
	delete[] rows;
	writeData();
}

void LineTable::addData(QString **Data, int NR, int NC)
{
	int r, n, C = mCore->columnCount(), N = mCore->rowCount();
	table->blockSignals(true);
	mCore->blockSignals(true);
	LineTableCore* ltc = reinterpret_cast<LineTableCore*>(mCore);
	if (C < NC) ltc->setColumnCount(NC);
	std::vector<BaseData*> data;
	for (n=0; n < NR; ++n) data.push_back(ltc->convertStringArrayToLineTableBaseData(Data[n], NC));
	ltc->insertRows(N - NpL, data);
	Destroy(Data, NR);
	mCore->blockSignals(false);
	table->blockSignals(false);
	Changed();
}

void LineTable::AddMarked(int AnzahlMarker, Marker *marker, Marker *LaserLine, QString SpektFile)
{
	//printf("LineTable::AddMarked Anzahl Zeilen=%d\n", Tab->rowCount());
	AcceptAssignments(SpektFile, false);
    if (marker == NULL) return;
	int i, j, n, N, k=0, l, c, RC = mCore->rowCount(), ColC = mCore->columnCount();
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
	//printf("Vor Schleife\n");
	for (N=0, i=0; i < AnzahlMarker; ++i)
		if (marker[i].Marked && marker[i].LState == LState && marker[i].UState == UState) ++N;
		//if (i == 19765 || i == 19766) printf("new lines found %d!\n", i);
	//}
	//printf("N=%d, LState=%d, UState=%d\n", N, LState, UState);
	Marker *Lines[N], *LBuff = marker;
	DeselectEverything();
	LineTableCore* ltc = reinterpret_cast<LineTableCore*>(mCore);
	for (i=j=0; i < AnzahlMarker; ++i)
		if (marker[i].Marked && marker[i].LState == LState && marker[i].UState == UState)
	{ 
		ltc->getMarkedLine(marker[i], SFN);
		if (n < RC) MarkLines(&(k=n), 1);
		else Lines[j++] = marker + i;
		marker[i].added = true;
	}
	if (j==0)
	{
		scrollTo(k);
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
	table->blockSignals(true);
	mCore->blockSignals(true);
 	mCore->setRowCount((i = RC) + N);
	ltc->AddMarked(Lines, LaserLine, AnzahlMarker, i, NpProg, MaxPN, SpektFile);
	scrollTo(i-1);
	mCore->blockSignals(false);
	table->blockSignals(false);
	Changed();
	//printf("Ende Addmarked, Anzahl Zeilen=%d\n", Tab->rowCount());
}

void LineTable::AcceptAssignments(QString SpektFile, bool accept)
{
	//printf("LineTable::AcceptAssignments\n");
	int NR = mCore->rowCount(), NC = mCore->columnCount(), r1, r2, r3, rd, c, n, m, PN, lPN;
	LineTableCore* ltc = reinterpret_cast<LineTableCore*>(mCore);
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
		for (r1 = r2 = NR - NpL; r2 < NR && ltc->getSourceFile(r2) != SpektFile; ++r2) ;
		for (r3 = r2; r3 < NR && ltc->getSourceFile(r3) == SpektFile; ++r3) ;
		if (r2 > r1 && r3 > r2)
		{
			BaseData** B = new BaseData*[rd = r3 - r2];
			for (n=0, m = r2; n < rd; n++, m++) B[n] = ltc->getData(m);
			for (n = r3 - 1, m = r2 - 1; m >= r1; n--, m--) ltc->setRow(ltc->getData(m), n);
			for (n = r1, m=0; m < rd; n++, m++) ltc->setRow(B[m], n);
			r3 = (r2 = r1) + rd;
			delete[] B;
		}
		if (r2 < r3) for (n=r2, PN = ltc->getProgression(n), ++MaxPN, --NpProg; n < r3; ++n)
		{
			if ((lPN = ltc->getProgression(n)) != PN)
			{
				++MaxPN;
				--NpProg;
			}
			ltc->setProgression(n, MaxPN);
		}
		NpL += r2 - r3;
		emit DataChanged();
	}
	else
	{
		//printf("NR = %d, NpL = %d\n", NR, NpL);
		for (r1 = NR - NpL; r1 < NR && ltc->getSourceFile(r1) != SpektFile; ++r1) ;
		for (r2 = r1; r2 < NR && ltc->getSourceFile(r2) == SpektFile; ++r2) ;
		for (n = r1, m = r2; m < NR; ++n, ++m) ltc->removeRow(n);
		NpL += r1 - r2;
	}
	//printf("Ende AcceptAssignments\n");
}

void LineTable::TabSelChanged()
{
	LineTableCore* ltc = reinterpret_cast<LineTableCore*>(mCore);
	int r, n, N, *Js, MJ, Mv, NIso, Jss, I, vss, *Rows, NR;
	bool UTA = true;
	double *E, ****ELU = 0, uE;
	ElState *LS;
	TermTable *T;
	table->getSelectedRows(Rows, NR);
	if (ltc->columnCount() <= LineTableCore::CEUp) UTA = false;
	else for (r=N=0; r < NR; ++r)
	{
		if (ltc->getUpperEnergy(Rows[r]) > 0.0) ++N;
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
			for (r=N=0; r < NR; ++r)
				if ((I = (ltc->getIso(Rows[r]) - 1) / 10) < NIso
						&& (vss = ltc->get_vss(Rows[r])) <= Mv
						&& (Jss = ltc->getJss(Rows[r])) < MJ)
					if (vss >= 0 && Jss >= 0 && I >= 0 && ELU[0][I][vss][Jss] != 0.0) ++N;
		}
		else if (N>0) UTA = true;
	}
	if (N>0 && ELU != 0)
	{
		Js = new int[N];
		E = new double[N];
		if (UTA)
		{
			for (r=N=0; r < NR; ++r)
				if ((uE = ltc->getUpperEnergy(Rows[r])) > 0.0)
			{
				Js[N] = ltc->getJs(Rows[r]);
				E[N++] = ltc->getUpperEnergy(Rows[r]);
			}
		}
		else for (r=N=0; r < NR; ++r)
			if ((I = (ltc->getIso(Rows[r]) - 1) / 10) < NIso && (vss = ltc->get_vss(Rows[r])) <= Mv
					&& (Jss = ltc->getJss(Rows[r])) < MJ)
				if (vss >= 0 && Jss >= 0 && I >= 0 ? ELU[0][I][vss][Jss] != 0.0 : false)
		{
			Js[N] = ltc->getJs(Rows[r]);
			E[N++] = ltc->getWaveNumber(Rows[r]) + ELU[0][I][vss][Jss];
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
	delete[] Rows;
}

void LineTable::getSelData(int *&Js, double *&E, int &N)
{
	Js = SelJs;
	E = SelE;
	N = NSel;
}

void LineTable::getViewnE(int*& Js, double*& E, int& N)
{
	LineTableCore* ltc = reinterpret_cast<LineTableCore*>(mCore);
	int NR = ltc->rowCount();
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
	for (n=N=0; n < NR; ++n) if (VR[n]) ++N;
	Js = new int[N];
	E = new double[N];
	for (n=N=0; n < NR; ++n) if (VR[n])
	{
		i = ltc->getIso(n);
		J = ltc->getJss(n);
		v = ltc->get_vss(n);
		if (i>=0 && i < NIso && J>=0 && J <= MJ && v >= 0 && v <= Mv)
		{
			Js[N] = ltc->getJs(n);
			E[N++] = ltc->getWaveNumber(n) + ELU[0][i][v][J];
		}
	}
	delete[] VR;
}

bool LineTable::ShowUpTerm()
{
	//printf("LineTable::ShowUpTerm\n");
	if (molecule == nullptr || transition == nullptr)
	{
		QMessageBox::information(this, "QT4MolSpektAn", tr("The line table has to be assigned to a molecule and a transition first!"));
		return false;
	}
	ElState *LS = transition->getLowerState();
	if (LS == nullptr)
	{
		QMessageBox::information(this, "QT4MolSpektAn", tr("The transition has to be assigned to a lower electronic state first!"));
		return false;
	}
	TermTable *TT = LS->getTermTable();
	DunTable *DT;
	if (TT == nullptr)
	{
		DT = LS->getDunTable();
		if (DT != nullptr) DT->calcTermEnergies(TT);
	}
	if (TT == nullptr)
	{
		QMessageBox::information(this, "QT4MolSpektAn", tr("A term energy table or a Dunham coefficient set for the lower state has to be loaded first!"));
		return false;
	}
	double ****UT = nullptr, ****ELU = TT->getData(), mvs = 0, Jeo = 0;
	int MJ = TT->getMaxJ(), Mv = TT->getMaxv(), *CT, MCT, *IsoT = TT->getIsoT();
	int *UIsoT = nullptr, *UCT = nullptr, UMCT = 0, UNC = 0;
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
	int n, j, N = mCore->rowCount(), NI = (Iso != 0 ? Iso->numIso : 0);
    int **Z = CreateInt(N, 6);
	QString SpektFile, Buffer;
	if (SO != 0) delete[] SO;
	int *LSO = heapSort(isnSPG);
	SO = new int[N];
	for (n=0; n<N; n++) SO[LSO[n]] = n;
	delete[] LSO;
	NSO = N;
	table->blockSignals(true);
	mCore->blockSignals(true);
	LineTableCore* ltc = reinterpret_cast<LineTableCore*>(mCore);
	if (UT != 0) ltc->setColumnCount(LineTableCore::COmC+1);
	else ltc->setColumnCount(LineTableCore::CEdJ+1);
	//printf("Vor Schleife1\n");
	ltc->ShowUpTerm(SO, MCT, CT, IsoT, ELU, Mv, MJ, UT, mvs, Jeo, UIsoT, S, UNC, UMCT, UCT);
	//printf("Vor Destroy\n");
	mCore->blockSignals(false);
	table->blockSignals(false);
    Destroy(Z, N);
	//printf("Ende ShowUpTerm\n");
	Changed();
	return true;
}

void LineTable::FindBigDiff()
{
	bool ChDD = false;
	LineTableCore* ltc = reinterpret_cast<LineTableCore*>(mCore);
	int i, n, N = ltc->rowCount(), nc = ltc->columnCount(), NR, *Rows;
	table->getSelectedRows(Rows, NR);
	//printf("nc=%d, COmC=%d\n", nc, COmC);
	if (NSO != N) ShowUpTerm();
	QString File, FBuffer;
	double EAV;
	//printf("nc=%d, COmC=%d\n", nc, COmC);
    if (++lRow >= N - 1) lRow = 0;
    for (i=0; i < NR; ++i) if (lRow == Rows[i])
    {
		EAV = ltc->getAverageUpperEnergy(lRow);
		File = ltc->getSourceFile(lRow);
    }
    while ((lRow<N && (fabs(ltc->getDiffToAverageUpperEnergy(lRow))) < 0.03 && (nc >= 14 && ChDD && fabs(ltc->getDeviationToCalculatedUpperEnergy(lRow)) < 0.1))
			|| (ltc->getAverageUpperEnergy(lRow) == EAV && ltc->getSourceFile(lRow) == File)) ++lRow;
    if (lRow < N) 
    {
		MarkLines(Rows, NR);
		File = ltc->getSourceFile(lRow);
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
		for (n=lRow, EAV = ltc->getAverageUpperEnergy(lRow); EAV == ltc->getAverageUpperEnergy(n); --n) ;
		++n;
		//Tab->ensureCellVisible(lRow, CEUma);
		while (ltc->getAverageUpperEnergy(n) == EAV)
		{
		    if ((FBuffer = ltc->getSourceFile(n)).left(6) == " laser")
				for (i=6; FBuffer[i] == ' '; ++i) ;
	    	else i = 1;
		    FBuffer = FBuffer.right(FBuffer.length() - i);
		    if (FBuffer == File) MarkLines(&n, 1);
		    //printf("FBuffer=%s, File=%s\n", FBuffer.ascii(), File.ascii());
	    	++n;
		}
		//MarkSelected();
    }
    //printf("lRow = %d\n", lRow);
    delete[] Rows;
}

void LineTable::ShowWeakProgressions()
{
	LineTableCore* ltc = reinterpret_cast<LineTableCore*>(mCore);
	int NRow = ltc->rowCount(), fRow = 0, aRow, av = -1, aJ = -1, lv = 0, lJ = 0;
	QString aFile, lFile;
	bool Weak = false;
	if (lRow >= NRow) aRow = 0;
	else aRow = lRow + 1;
	for (;aRow != lRow; ++aRow)
	{
		if (aRow == NRow) aRow = 0;
		aFile = ltc->getSourceFile(aRow);
		aJ = ltc->getJs(aRow);
		av = ltc->get_vs(aRow);
		if (av != lv || aJ != lJ || aFile != lFile)
		{
			if (Weak)
			{
				for (lRow = fRow; lRow < aRow; ++lRow) MarkLines(&lRow, 1);
				//Tab->ensureCellVisible(lRow, CSNR);
				return;
			}
			else Weak = true;
			lJ = aJ;
			lv = av;
			lFile = aFile;
			fRow = aRow;
		}
		else if ((ltc->getSNR(aRow) > 4 && ltc->getComment(aRow).isEmpty())
				|| ltc->getComment(aRow) == "satellite") Weak = false;
	}
	QMessageBox::information( this, "MolSpektAnalysis", "There are no really weak progressions found inside the list.", QMessageBox::Ok);
}

TableWindow* LineTable::ShowUpTermTable()
{
    if (MW == 0) return 0;
	int n;
    QString **Data;
    if ((mCore->columnCount() == LineTableCore::TableNormCols || NSO != mCore->rowCount()) && !ShowUpTerm()) return nullptr;
    if (termTable == nullptr)
    {
		termTable = new TableWidgetWindow(TextTable1, MW, molecule);
		termTable->setWindowTitle("UpTermTable to " + getName());
		termTable->setHorizontalHeader(QStringList() << "Isotop" << "v" << "J" << "Par" 
				<< "Termenergie" << "Standardabweichung" << "Anzahl Uebergaenge" 
				<< "Anz. vollst. Doubl." << "Min E" << "Max E");
    }
	reinterpret_cast<LineTableCore*>(mCore)->ShowUpTermtable(SO, n, Data);
    termTable->setData(Data, n, 10);
	Destroy(Data, NSO);
    return termTable;
}

void LineTable::Updatevs(int *nvs)
{
	int i, nr = mCore->rowCount();
	QString Buffer;
	table->blockSignals(true);
	mCore->blockSignals(true);
	for (i=0; i < nr; ++i) reinterpret_cast<LineTableCore*>(mCore)->set_vs(i, nvs[i]);
	mCore->blockSignals(false);
	table->blockSignals(false);
	Changed();
}

void LineTable::RemoveDoubled()
{
	int i, j, nr = mCore->rowCount(), *S1 = heapSort(sortfRemDoubl);
	int *S = new int[nr];
	QString F1, F2, Comment;
	for (i=0; i < nr; ++i) S[S1[i]] = i;
	table->blockSignals(true);
	mCore->blockSignals(true);
	reinterpret_cast<LineTableCore*>(mCore)->RemoveDoubled(S);
	mCore->blockSignals(false);
	table->blockSignals(false);
	Changed();
	delete[] S1;
	delete[] S;
} 

void LineTable::setData(const std::vector<LineTableCore::TableCols>& mColumns, const std::vector<BaseData*> data)
{
	table->blockSignals(true);
	mCore->blockSignals(true);
	reinterpret_cast<LineTableCore*>(mCore)->setData(mColumns, data);
	mCore->blockSignals(false);
	table->blockSignals(false);
	Changed();
}

void LineTable::splitTable()
{
	if (MW == 0) return;
	LineTable *NT = MW->CreateLineTable();
	if (NT == 0) return;
	LineTableCore* ltc = reinterpret_cast<LineTableCore*>(mCore);
	int i, j, k, c, R = ltc->rowCount(), C = ltc->columnCount();
	std::vector<BaseData*> data;
	ElState *lS, *uS;
	if (transition != 0) 
	{
		lS = transition->getLowerState();
		uS = transition->getUpperState();
	}
	else lS = uS = 0;
	table->blockSignals(true);
	ltc->blockSignals(true);
	for (k=0; k<R && ltc->get_vs(k) != -1; ++k) ;
	for ( ; k<R; ++k)
	{
		if (ltc->get_vs(k) == -1)
		{
			data.push_back(ltc->getRow(k));
			mCore->setData(k, nullptr);
		}
	}
	ltc->RemoveEmptyRows();
	ltc->blockSignals(false);
	table->blockSignals(false);
	NT->setData(ltc->getColumnVector(), data);
	NT->setSource(getSource());
	NT->setName(getName());
	if (molecule != 0) molecule->addLineTable(NT, lS);
	NT->show();
	Changed();
}

void LineTable::sortUpTermIvJ()
{
	if (mCore->columnCount() <  LineTableCore::CEav)
	{
		QMessageBox::information( this, "MolSpektAnalysis", 
			"There is no upper term energy table calculated for this line table!", QMessageBox::Ok);
		return;
	}
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
	LineTableCore* ltc = reinterpret_cast<LineTableCore*>(mCore);
	table->blockSignals(true);
	ltc->blockSignals(true);
	QString F1, F2, FN;
	int l, i, n, N = ltc->rowCount();
	for (n=0; n<N; ++n)
	{
		l = (F2 = ltc->getSourceFile(n)).length();
		if (l==0) F2 = F1;
		else
		{
			for (i=0; i<l && F2[i] == QChar(' '); ++i) ;
			F2 = F2.right(l-=i);
			if (F2.left(5) == "laser")
			{
				F2 = F2.right(l-=5);
				for (i=0; i<l && F2[i] == QChar(' '); ++i) ;
				F2 = F2.right(l-i);
				if (F2.isEmpty()) F2 = F1;
				ltc->setComment(n, "laser");
			}
			F1 = F2;
		}
		ltc->setSourceFile(n, F2);
	}
	ltc->sortTab(S2);
	ltc->blockSignals(false);
	table->blockSignals(false);
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
	int i, n = mCore->rowCount(), r, m=0, *rows, NR;
	LineTableCore* ltc = reinterpret_cast<LineTableCore*>(mCore);
	table->getSelectedRows(rows, NR);
	table->blockSignals(true);
	mCore->blockSignals(true);
	if (MaxMin->isChecked())
	{
		if (NR > 0) for (r=0; r < NR; ++r) ltc->setUncertainty(rows[r], Error, OvError);
		else for (i=0; i<n; ++i) ltc->setUncertainty(i, Error, OvError);
	}
	else for (r=0; r < NR; ++r)
	{
		QString Buffer = ltc->getComment(rows[r]);
		if (Buffer.indexOf("overlap") > -1 || Buffer.indexOf("laser") > -1 || ltc->getSNR(rows[r]) < 5.0) ltc->setUncertainty(rows[r], OvError);
		else ltc->setUncertainty(rows[r], Error);
	}
	mCore->blockSignals(false);
	table->blockSignals(false);
	Changed();
	delete[] rows;
}

void LineTable::Shiftvup()
{
	int r, *rows, NR;
	LineTableCore* ltc = reinterpret_cast<LineTableCore*>(mCore);
	table->blockSignals(true);
	mCore->blockSignals(true);
	for (r=0; r < NR; ++r) ltc->set_vss(rows[r], ltc->get_vss(rows[r]) + 1);
	mCore->blockSignals(false);
	table->blockSignals(false);
	Changed();
	delete[] rows;
}

void LineTable::Shiftvdown()
{
	int r, *rows, NR;
	LineTableCore* ltc = reinterpret_cast<LineTableCore*>(mCore);
	table->blockSignals(true);
	mCore->blockSignals(true);
	for (r=0; r < NR; ++r) ltc->set_vss(rows[r], ltc->get_vss(rows[r]) - 1);
	mCore->blockSignals(false);
	table->blockSignals(false);
	Changed();
	delete[] rows;
}

void LineTable::ShiftJup()
{
	int r, *rows, NR;
	LineTableCore* ltc = reinterpret_cast<LineTableCore*>(mCore);
	table->blockSignals(true);
	mCore->blockSignals(true);
	for (r=0; r < NR; ++r)
	{
		ltc->setJs(rows[r], ltc->getJs(rows[r]) + 1);
		ltc->setJss(rows[r], ltc->getJss(rows[r]) + 1);
	}
	mCore->blockSignals(false);
	table->blockSignals(false);
	Changed();
	delete[] rows;
}

void LineTable::ShiftJdown()
{
	int r, *rows, NR;
	LineTableCore* ltc = reinterpret_cast<LineTableCore*>(mCore);
	table->blockSignals(true);
	mCore->blockSignals(true);
	for (r=0; r < NR; ++r)
	{
		ltc->setJs(rows[r], ltc->getJs(rows[r]) - 1);
		ltc->setJss(rows[r], ltc->getJss(rows[r]) - 1);
	}
	mCore->blockSignals(false);
	table->blockSignals(false);
	Changed();
	delete[] rows;
}

void LineTable::Shiftvsdown()
{
	int r, *rows, NR;
	LineTableCore* ltc = reinterpret_cast<LineTableCore*>(mCore);
	table->blockSignals(true);
	mCore->blockSignals(true);
	for (r=0; r < NR; ++r) ltc->set_vs(rows[r], ltc->get_vs(rows[r]) - 1);
	mCore->blockSignals(false);
	table->blockSignals(false);
	Changed();
	delete[] rows;
}

void LineTable::Shiftvsup()
{
	int r, *rows, NR;
	LineTableCore* ltc = reinterpret_cast<LineTableCore*>(mCore);
	table->blockSignals(true);
	mCore->blockSignals(true);
	for (r=0; r < NR; ++r) ltc->set_vs(rows[r], ltc->get_vs(rows[r]) + 1);
	mCore->blockSignals(false);
	table->blockSignals(false);
	Changed();
	delete[] rows;
}

void LineTable::ShiftIso()
{
	int M = ((Iso != 0 ? Iso->numIso : 0) - 1) * 10;
	if (M == -10)
	{
		QMessageBox::information(this, "MolSpektAnalysis", 
								 "This line table has to be assigned to a molecule first!",
								 QMessageBox::Ok);
		return;
	}
	int r, v, *rows, NR;
	LineTableCore* ltc = reinterpret_cast<LineTableCore*>(mCore);
	table->blockSignals(true);
	mCore->blockSignals(true);
	for (r=0; r < NR; ++r)
	{
		if ((v = ltc->getIso(rows[r])) <= M) ltc->setIso(rows[r], v + 10);
		else ltc->setIso(rows[r], v - M);
	}
	mCore->blockSignals(false);
	table->blockSignals(false);
	Changed();
	delete[] rows;
}

void LineTable::SetvssAscending()
{
	LineTableCore* ltc = reinterpret_cast<LineTableCore*>(mCore);
	int i, l=-1, *rows, N;
	table->getSelectedRows(rows, N);
	if (N==0)
	{
		QMessageBox::information(this, "MolSpektAnalysis", "You have to select some rows first!",
								 QMessageBox::Ok);
		return;
	}
	table->blockSignals(true);
	mCore->blockSignals(true);
	for (i=0; i<N; ++i)
	{
		if (l==-1) l = ltc->get_vss(rows[i]) + 1;
		else ltc->set_vss(rows[i], l++);
	}
	mCore->blockSignals(false);
	table->blockSignals(false);
	Changed();
	delete[] rows;
}

void LineTable::Delete()
{
	QModelIndexList list = table->getSelectedIndexes();
	table->blockSignals(true);
	mCore->blockSignals(true);
	for (auto it = list.begin(); it != list.end(); ++it) mCore->setRow(nullptr, it->row());
	mCore->RemoveEmptyRows();
	mCore->blockSignals(false);
	table->blockSignals(false);
	Changed();
}

void LineTable::cutRows(int& numRows, int& numColumns, QString**& Data)
{
    if (NpL > 0)
	{
		int NR = mCore->rowCount();
		int SR = NR - NpL, n, r, *rows, N;
		table->getSelectedRows(rows, N);
		bool Rem[NpL];
		for (n=0; n < NpL; ++n) Rem[n] = false;
		for (n=0; n < N; ++n) if (rows[n] >= SR) Rem[rows[n] - SR] = true;
		for (n=r=0; n < NpL; ++n) if (Rem[n]) ++r;
		NpL -= r;
		delete[] rows;
	}
	TableViewWindow::cutRows(numRows, numColumns, Data);
}

void LineTable::DeleteRows()
{
	int *rows, N;
	table->getSelectedRows(rows, N);
	deleteRows(rows, N);
	delete[] rows;
}

void LineTable::deleteRows(int* rows, int N)
{
	int j, k, n = mCore->rowCount(), c, nc = mCore->columnCount();
	int pr = n - NpL;
	QString B, Buffer;
	table->blockSignals(true);
	mCore->blockSignals(true);
	for (j=0; j<N; ++j)
	{
		if (!(B = mCore->getSourceFile(rows[j])).isEmpty()) Buffer = B;
		mCore->setRow(nullptr, rows[j]);
		if (rows[j] >= pr) NpL--;
		k = rows[j] + 1;
		if (k<n && mCore->getSourceFile(k).isEmpty()) mCore->setSourceFile(k, Buffer);
	}
	for (j=k=0; j <= lRow; ++j) if (nullptr == mCore->getRow(j)) ++k;
	lRow -= k;
	mCore->RemoveEmptyRows();
	mCore->blockSignals(false);
	table->blockSignals(false);
	Changed();
}

void LineTable::DeleteRows(int NL, int* PN, int* vss, int* Jss)
{
	LineTableCore* ltc = reinterpret_cast<LineTableCore*>(mCore);
	int R = mCore->rowCount(), r, C = mCore->columnCount(), c, n, sr = R - NpL;
	table->blockSignals(true);
	mCore->blockSignals(true);
	for (r=0; r<R; ++r)
	{
		for (n=0; n < NL && (mCore->getProgression(r) != PN[n] || (vss != 0 && Jss != 0 && (ltc->get_vss(r) != vss[n] || ltc->getJss(r) != Jss[n]))); ++n) ;
		if (n < NL)
		{
			ltc->setRow(nullptr, r);
			if (r >= sr) --NpL;
		}
	}
	ltc->RemoveEmptyRows();
	mCore->blockSignals(false);
	table->blockSignals(false);
	Changed();
}

void LineTable::SetError(int NL, int* PN, int* vss, int* Jss, double* Err)
{
	LineTableCore* ltc = reinterpret_cast<LineTableCore*>(mCore);
	int R = ltc->rowCount(), r, n;
	table->blockSignals(true);
	mCore->blockSignals(true);
	for (r=0; r<R; ++r)
	{
		for (n=0; n < NL && (ltc->getProgression(r) != PN[n] || ltc->get_vss(r) != vss[n] || ltc->getJss(r) != Jss[n]); ++n) ;
		if (n < NL) ltc->setUncertainty(r, Err[n]);
	}
	mCore->blockSignals(false);
	table->blockSignals(false);
	Changed();
}

Progression LineTable::getSelectedProgression()
{
	LineTableCore* ltc = reinterpret_cast<LineTableCore*>(mCore);
	Progression P;
	int L[1000], r = table->currentIndex().row(), R = ltc->rowCount(), n;
	if (r == -1)
	{
		P.N = 0;
		return P;
	}
	int I = P.Iso = ltc->getIso(r);
	int J = P.Js = ltc->getJs(r);
	int v = P.vs = ltc->get_vs(r), PN = ltc->getProgression(r);
	QString F = ltc->getSourceFile(r);
	for (r=n=0; r<R; ++r) if (ltc->getProgression(r) == PN && ltc->getSourceFile(r) == F && ltc->getIso(r) == I && ltc->get_vs(r) == v && ltc->getJs(r) == J) L[n++] = r;
	P.L = new Line[P.N = n];
	for (n=0; n < P.N; ++n)
	{
		P.L[n].E = ltc->getWaveNumber(L[n]);
		P.L[n].Jss = ltc->getJss(L[n]);
		P.L[n].err = ltc->getUncertainty(L[n]);
		P.L[n].vss = ltc->get_vss(L[n]);
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
	LineTableCore* ltc = reinterpret_cast<LineTableCore*>(mCore);
	int n, mv, Mv = P.L[0].vss + 1, N, LM, M = (P.Js == P.L[0].Jss ? P.N - 1 : P.N / 2 - 1);
	int v[1000], vs, Js, I, PN, RC = mCore->rowCount(), cr, m, nwL;
	double S[1000], B;
	QString F;
	for (n=0, mv = Mv - 2; n < P.N; n++)
	{
		if (P.L[n].vss <= mv) mv = P.L[n].vss - 1;
		else if (P.L[n].vss >= Mv) Mv = P.L[n].vss + 1;
	}
	n = cr = table->currentIndex().row();
	F = ltc->getSourceFile(n);
	vs = ltc->get_vs(n);
	Js = ltc->getJs(n);
	I = ltc->getIso(n);
	PN = ltc->getProgression(n);
	for (++n; n < RC && ltc->getSourceFile(n) == F && ltc->get_vs(n) == vs && ltc->getJs(n) == Js && ltc->getIso(n) == I && ltc->getProgression(n) == PN; ++n) ;
	for (F = ltc->getSourceFile(n), vs = ltc->get_vs(n), Js = ltc->getJs(n), I = ltc->getIso(n), PN = ltc->getProgression(n), N=0; n != cr; ++n)
	{
		if (n == RC) 
		{
			n=0;
			if (cr == 0) break;
		}
		if (ltc->getSourceFile(n) != F || ltc->get_vs(n) != vs || ltc->getJs(n) != Js || ltc->getIso(n) != I || ltc->getProgression(n) != PN)
		{
			LM = ((ltc->getJss(n!=0 ? n : RC) - 1) == Js ? M : 2 * M);
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
				vs = ltc->get_vss(v[m]);
				F = ltc->getComment(v[m]);
				if ((vs < mv || vs > Mv) && F.indexOf("overlap") == -1 && F.indexOf("laser") == -1) ++nwL;
			}
			if (nwL <= 1)
			{
				vs = (n!=0 ? n : RC);
				for (m = vs - N; m < vs; ++m) MarkLines(&m, 1);
				return;
			}
			vs = ltc->get_vs(n);
			Js = ltc->getJs(n);
			F = ltc->getSourceFile(n);
			I = ltc->getIso(n);
			PN = ltc->getProgression(n);
			N=0;
		}
		S[N] = ltc->getSNR(n);
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
	int vs = NE.text().toInt();
	int r, *rows, N;
	LineTableCore* ltc = reinterpret_cast<LineTableCore*>(mCore);
	table->getSelectedRows(rows, N);
	table->blockSignals(true);
	mCore->blockSignals(true);
	for (r=0; r < NR; ++r) ltc->set_vs(rows[r], vs);
	mCore->blockSignals(false);
	table->blockSignals(false);
	Changed();
	delete[] rows;
}

void LineTable::setFC(const int FC)
{
	int r,*rows, N;
	table->getSelectedRows(rows, N);
	LineTableCore* ltc = reinterpret_cast<LineTableCore*>(mCore);
	table->blockSignals(true);
	mCore->blockSignals(true);
	for (r=0; r<N; ++r) ltc->setFineStructureQN(r, FC);
	mCore->blockSignals(false);
	table->blockSignals(false);
	Changed();
	delete[] rows;
}

void LineTable::sortbyvs()
{
	int *SArray = heapSort(sortByvs);
	sortTab(SArray);
}

void LineTable::SetPN()
{
	LineTableCore* ltc = reinterpret_cast<LineTableCore*>(mCore);
	int n, aJ = -1, nJ, ndJ, adJ = -1, aI=-1, nI, av=-1, nv, aF = -1, nF, N = ltc->rowCount(), i, j;
	int *SA = new int[N];
	double nWn, aWn = 0.0;
	QString aFi, nFi;
	bool UTA = ShowUpTerm();
	int *tSO = heapSort(sortForSPN);
	for (n=0; n<N; n++) SA[tSO[n]] = n;
	delete[] tSO;
	MaxPN = 0;
	table->blockSignals(true);
	mCore->blockSignals(true);
	for (n=0; n < N; ++n)
	{
		nJ = ltc->getJs(SA[n]);
		ndJ = nJ - ltc->getJss(SA[n]);
		ndJ *= ndJ;
		nI = ltc->getIso(SA[n]);
		nv = ltc->get_vs(SA[n]);
		nWn = (UTA ? ltc->getUpperEnergy(SA[n]) : 0.0);
		i = (nFi = ltc->getSourceFile(SA[n])).lastIndexOf(QRegExp("[\\/]")) + 1;
		j = nFi.indexOf('.', i);
		nFi = nFi.mid(i, (j>=0 ? j : nFi.length()) - i);
		nF = ltc->getFineStructureQN(SA[n]);
		if (nFi.isEmpty()) ltc->setSourceFile(SA[n], nFi = aFi);
		if (nJ != aJ || ndJ != adJ || nv != av || fabs(nI - aI) > 9.0 || nF != aF || nFi != aFi || fabs(nWn - aWn) > 1e2 * ltc->getUncertainty(SA[n]))
		{
			aJ = nJ;
			av = nv;
			aI = nI;
			aF = nF;
			aFi = nFi;
			aWn = nWn;
			adJ = ndJ;
			++MaxPN;
		}
		ltc->setProgression(SA[n], MaxPN);
	}
	mCore->blockSignals(false);
	table->blockSignals(false);
	delete[] SA;
	Changed();
}

void LineTable::getKnownLevels(int NI, int &mvs, int &mvss, int &mJs, int &mJss, bool ***&uL, bool ***&lL)
{
	//printf("LineTab::getKnownLevels\n");
	LineTableCore* ltc = reinterpret_cast<LineTableCore*>(mCore);
	int n, m, i, o, I, NR = ltc->rowCount();
	for (n=0, mvs = 0, mvss = 0, mJs = 0, mJss = 0; n < NR; ++n)
	{
		if ((m = ltc->get_vs(n)) > mvs) mvs = m;
		if ((m = ltc->get_vss(n)) > mvss) mvss = m;
		if ((m = ltc->getJs(n)) > mJs) mJs = m;
		if ((m = ltc->getJss(n)) > mJss) mJss = m;
	}
	uL = CreateBool(NI, mvs + 1, mJs + 1);
	lL = CreateBool(NI, mvss + 1, mJss + 1);
	for (i=0; i < NI; ++i)
	{
		for (n=0; n <= mvs; ++n) for (m=0; m <= mJs; ++m) uL[i][n][m] = false;
		for (n=0; n <= mvss; ++n) for (m=0; m <= mJss; ++m) lL[i][n][m] = false;
	}
	for (i=0; i < NR; ++i)
	{
		I = (o = ltc->getIso(i)) / 10;
		if (o == 10 * I && I > 0 && I <= NI)
		{
			--I;
			if ((m = ltc->get_vs(i)) >= 0 && (n = ltc->getJs(i)) >= 0) uL[I][m][n] = true;
			if ((m = ltc->get_vss(i)) >= 0 && (n = ltc->getJss(i)) >= 0) lL[I][m][n] = true;
		}
	}
	//printf("NR=%d, mJs=%d, mJss=%d\n", NR, mJs, mJss);
}

void LineTable::findErrors()
{
	int n, N = mCore->columnCount();
	for (n=0; n<N && mCore->getUncertainty(n) != 0.0; ++n) ;
	if (n<N) MarkLines(&n, 1);
}

void LineTable::HeaderItemDoubleClicked(const int index)
{
    if (lastClickedHeaderIndex != index) return;
	const std::vector<LineTableCore::TableCols>& columns = reinterpret_cast<LineTableCore*>(mCore)->getColumnVector();
    switch(columns[index])
    {
		case LineTableCore::CPN:
            sortTab(heapSort(sortByProgression));
            break;
		case LineTableCore::Cvs:
            sortTab(heapSort(sortByvs));
            break;
		case LineTableCore:: CJs:
            sortTab(heapSort(sortByJs));
            break;
		case LineTableCore::Cvss:
            sortTab(heapSort(sortBy_vss));
            break;
		case LineTableCore::CJss:
            sortTab(heapSort(sortByJss));
            break;
		case LineTableCore:: CF:
            sortTab(heapSort(sortByF));
            break;
		case LineTableCore::CWN:
            sortTab(heapSort(sortByFrequency));
            break;
        case LineTableCore::Cerr:
            sortTab(heapSort(sortBy_err));
            break;
        case LineTableCore::CIso:
            sortTab(heapSort(sortByIso));
            break;
        case LineTableCore::CFile:
            sortTab(heapSort(sortByFile));
            break;
        case LineTableCore::CSNR:
            sortTab(heapSort(sortBySNR));
            break;
        case LineTableCore::CDev:
            sortTab(heapSort(sortByDev));
            break;
        case LineTableCore::CC:
            sortTab(heapSort(sortByComment));
            break;
		case LineTableCore::CFCF:
			sortTab(heapSort(sortByFCF));
            break;
		case LineTableCore::CEUp:
			sortTab(heapSort(sortByEUp));
            break;
		case LineTableCore::CEav:
			sortTab(heapSort(sortByEav));
            break;
		case LineTableCore::CEUma:
			sortTab(heapSort(sortByEUma));
            break;
		case LineTableCore::CEdJ:
			sortTab(heapSort(sortByEdJ));
            break;
		case LineTableCore::CCalc:
			sortTab(heapSort(sortByCalc));
            break;
		case LineTableCore::COmC:
			sortTab(heapSort(sortByOmC));
            break;
        default:
            // should not happen
            break;
    }
}

void LineTable::shrinkAllSpectRefs()
{
	table->blockSignals(true);
	mCore->blockSignals(true);
	mCore->shrinkAllSpectRefs();
	table->blockSignals(false);
	mCore->blockSignals(false);
}

void LineTable::sortByProgNumber()
{
	sortTab(heapSort(sortByProgression));
}
