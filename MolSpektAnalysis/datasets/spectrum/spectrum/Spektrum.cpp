//
// C++ Implementation: spektrum
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2006 - 2019
//
// Copyright: See README file that comes with this source code
//
//

#include "Spektrum.h"

#include <qvariant.h>
#include <qstring.h>
#include <qfile.h>
#include <qmessagebox.h>
#include <qstringlist.h>
#include <QFileDialog>
#include <QKeyEvent>
#include <QPainter>
#include <QTextStream>
#include <QRadioButton>
#include <QCheckBox>
#include <QToolButton>

#include <math.h>

#include "inputdialog.h"
#include "tools.h"
#include "utils.h"
#include "TermPlot.h"
#include "elstate.h"
#include "termtable.h"
#include "molecule.h"
#include "linetable.h"
#include "constants.h"
#include "potential.h"
#include "fit.h"
#include "fcftab.h"
#include "fitobject.h"
#include "intprog.h"
#include "gaussian.h"
#include "datensatz.h"
#include "isotab.h"
#include "band.h"
#include "testprog.h"
#include "measuredtermenergies.h"
#include "tableline.h"


//_oD_Sr2"
#define OnlyDoublet false
//#define Strontium

Spektrum::Spektrum(MainWindow *MW) : DiagWindow(MDIChild::Spect, MW), Rauschen(-1.0), m_fitedLineVector()
{
	LaserFTol = 0.3;
    LaserFrequency = 0.0;
	Type = LaserInducedFluorescence;
	xStartLabel->setText("Frequency [cm-1]:    from ");
    xStart->setText("10000");
    xStopLabel->setText(" to ");
    xStop->setText("20000");
    yStartLabel->setText("Intensity:    from ");
    yStart->setText("0");
    yStopLabel->setText(" to ");
    yStop->setText("1000");
	
	setFilter("Spectrum (*.spect)");
	setFileExt(".spect");
	
    LineTablesAT = new LineTable*[MaxLineTables];
	numLineTablesAT = 0;
	DMarkers = false;
    AnzahlMarker = 0;
	LaserLinie = NULL;
    marker = NULL;
    ST = CMaxSearchDev;
    GMH = 5;
	GMarker = NULL;
    MinPeakHeight = 0.5;
    LPI = 1;
    DC = 0;
	AssignmentStatus = 0;
	ShowAssignmentsOnTop = false;
	double *Buffer;
	Buffer = new double[2];
	Buffer[0] = -9999999999.9;
	MBeginn.Line = Buffer;
	Buffer = new double[2];
	Buffer[0] = 999999999.9;
	MStop.Line = Buffer;
	useIntensities = true;
    connect(Bild, SIGNAL(SelectionChanged(QRect*)), this, SLOT(markRegion(QRect*)));
}

Spektrum::~Spektrum()
{
    //printf("Spektrum::destroy\n");
    if (marker != NULL) delete[] marker;
	delete[] MBeginn.Line;
	delete[] MStop.Line;
	if (GMarker != NULL) delete[] GMarker;
	delete[] LineTablesAT;
    for (size_t n =0; n < m_fitedLineVector.size(); ++n) delete m_fitedLineVector[n];
}

void Spektrum::addMarker(bool autoaccept)
{
	//printf("Spektrum::addMarker\n");
	int n, i;
	QString FileName = getFileName();
	LineTable *L, *NM = 0;
	for (n=0; n < numLineTablesAT; n++) LineTablesAT[n]->AcceptAssignments(FileName, false);
	numLineTablesAT = 0;
	for (n=0; n < AnzahlMarker; n++) if (marker[n].Marked) marker[n].added = false;
	for (n=0; n < AnzahlMarker; n++)
	{
		if (marker[n].Marked && !marker[n].added)
		{
			if (marker[n].Mol == 0) 
			{
				//printf("Spektrum::addMarker error: marker[%d].Mol=0!\n", n);
				if (MW != 0) NM = MW->getUnassignedLineTable();
				for (i=0; (i < numLineTablesAT ? LineTablesAT[i] != NM : false); i++) ;
				if (i == numLineTablesAT) LineTablesAT[numLineTablesAT++] = NM;
				if (NM != 0)
				{
					NM->AddMarked(AnzahlMarker, marker, LaserLinie, FileName);
					NM->show();
					MW->setActive(NM);
					MW->setActive(this);
				}
			}
			else if ((L = marker[n].Mol->getLineTable(marker[n].LState, marker[n].UState)) != 0)
			{
				L->AddMarked(AnzahlMarker, marker, LaserLinie, FileName);
				for (i=0; (i < numLineTablesAT ? LineTablesAT[i] != L : false); i++) ;
				if (i == numLineTablesAT) LineTablesAT[numLineTablesAT++] = L;
				if (autoaccept) L->AcceptAssignments(FileName, true);
				if (MW != 0) 
				{
					MW->setActive(L);
					MW->setActive(this);
				}
			}
		}
	}
}

void Spektrum::LineTableSaved(ElState *UState, ElState *LState)
{
	int n;
	for (n=0; n < AnzahlMarker; n++) 
		if (marker[n].Marked && marker[n].UState == UState && marker[n].LState == LState)
			marker[n].Marked = false;
	Paint();
}

void Spektrum::AcceptAssignments()
{
	int n;
	LineTable *L;
	for (n=0; n < AnzahlMarker; n++) if (marker[n].Marked)
	{
		if (marker[n].Mol == 0 && MW != 0) L = MW->getUnassignedLineTable();
		else L = marker[n].Mol->getLineTable(marker[n].LState, marker[n].UState);
		if (L != 0) L->AcceptAssignments(getFileName(), true);
	}
	ClearMarked();
	emit SpectrumChanged(this);
}	

bool Spektrum::GetShowAssignmentsOnTop()
{
	return ShowAssignmentsOnTop;
}

void Spektrum::SetShowAssignmentsOnTop(bool V)
{
	ShowAssignmentsOnTop = V;
	Paint();
}

void Spektrum::editFind()
{
	//printf("editFind\n");
	/*if (AnzahlMarker > 0) 
	{
		printf("Es sind bereits Marker vorhanden.\n");
       	printf("Das Einfügen weiterer würde das Spektrum verfälschen!\n");
   	}*/
   	double x1, x2, x3, y1, y2, y3, a = 0.0, c;
   	double m1, m2=0.0, d, aR;
	const double *Min1, *Min2 = NULL;
	int i, j=0, N, NAM = AnzahlMarker;
	Marker *AMarker = marker;
	N = Daten->GetDSL();
	if (N < 3) return;
	x2 = Daten->GetValue(0, 0);
   	y2 = Daten->GetValue(0, 1);
   	x3 = Daten->GetValue(1, 0);
   	y3 = Daten->GetValue(1, 1);
   	//printf("Find: Nach GetValue(), N=%d\n", N);
   	Datensatz Peaks;
   	for (i=2; i<N; i++)
   	{
    	x1 = x2;
    	y1 = y2;
    	x2 = x3;
    	y2 = y3;
    	x3 = Daten->GetValue(i, 0);
    	y3 = Daten->GetValue(i, 1);
	  	if (y1 < y2 && y2 >= y3 
				&& !Daten->GetMarked(i-1) && !Daten->GetMarked(i) && !Daten->GetMarked(i-2)) 
		{
			//if (y2 < 0.0) printf("y2<0!\n");
			//printf("j=%d, x2=%f\n", mc, x2);
			/*if (AnzahlMarker != 0)
			{
				printf("i=%d, x2=%f\n", i, x2);
				printf("x1=%f, y1=%f, x2=%f, y2=%f, x3=%f, y3=%f\n", x1, y1, x2, y2, x3, y3);
			}*/
			ParabInterpol(x1, y1, x2, y2, x3, y3, c, a);
			Peaks.AddValue(c, a, true);
			j++;
			//printf("i=%d\n", i);
			//if (y2 < 0) printf("x1=%f, x2=%f, x3=%f\ny1=%f, y2=%f, y3=%f\n a=%f, b=%f, c=%f\n", 
		  //x1, x2, x3, y1, y2, y3, a, b, c);
		}
    }
   	//printf("AnzahlMarker=%d, Peaks.GetDSL=%d\n", AnzahlMarker, Peaks.GetDSL());
	AnzahlMarker += Peaks.GetDSL();
   	//printf("j=%d, AnzahlMarker=%d, NAM=%dm DSL=%d\n", j, AnzahlMarker, NAM, Daten.GetDSL());
   	marker = new Marker[AnzahlMarker];
	for (i=NAM; i<AnzahlMarker; i++) Peaks.SetMarker(i - NAM, marker + i); //printf("i=%d\n", i);}
	Daten->InsertML(Peaks);
	//printf("Nach InsertML, i=%d, DSL=%d\n", i, Daten.GetDSL());
	N = Daten->GetDSL();
	if (NAM > 0) for (i=0, j=0; i<N; i++) if (Daten->GetMarked(i)) Daten->SetMarker(i, marker + j++);
	//printf("Nach setMarker, j=%d\n", j);
	if (LaserFrequency != 0)
	{
		for (i=0; (i < AnzahlMarker ? fabs(LaserFrequency - marker[i].Line[0]) >= 0.1 : false); i++) ;
		if (i < AnzahlMarker)
		{
			for (LaserLinie = marker + (i++); (i < AnzahlMarker ? fabs(LaserFrequency - marker[i].Line[0]) < 0.1 : false); i++)
				if (marker[i].Line[1] > LaserLinie->Line[1]) LaserLinie = marker + i;
		}
		else LaserLinie = 0;
	}
	else for (i=1, LaserLinie = marker; i<AnzahlMarker; i++) if (marker[i].Line[1] > LaserLinie->Line[1]) LaserLinie = marker + i;
	//printf("Vor LaserLine\n");
	//printf("LaserLinie->Line[1]=%f\n", LaserLinie->Line[1]);
	//printf("marker[7230].Line[0]=%f\n", marker[7230].Line[0]);
	if (!DMarkers) DisplayMarkers();
	x2 = Daten->GetValue(0, 0);
   	y2 = Daten->GetValue(0, 1);
   	x3 = Daten->GetValue(1, 0);
   	y3 = Daten->GetValue(1, 1);
	//printf("Vor zweiter Haupschleife\n");
	for (i=2, j=0; i < N && j < AnzahlMarker; i++)
	{ 
		//printf("i=%d, N=%d\n", i, N);
		x1 = x2;
       	y1 = y2;
       	x2 = x3;
       	y2 = y3;
       	x3 = Daten->GetValue(i, 0);
       	y3 = Daten->GetValue(i, 1);
		if (y1 > y2 && y2 <= y3)
		{
			//printf("mc=%d, x2=%f\n", mc, x2);
			m1 = m2;
			m2 = y2;
			Min1 = Min2;
			Min2 = Daten->GetPoint(i-1);
			if (x2 > marker[j].Line[0])
			{
				/*if (Min1 != 0 ? Min1[0] > marker[j].Line[0] : false) 
					printf("Fehler: i=%d, j=%d, x2=%f\n", i, j, x2); 
				if (j > 0 && j < 100 || j > 17900) 
					printf("j=%d, M1=%f, L=%f, M2=%f\n", j, Min1[0], marker[j].Line[0], Min2[0]);*/
				a = marker[j].Line[1];
				d = (a - m1) / (a - m2);
				if (d < 1.33 && d > 0.75) marker[j].sOverlap = false;
				else marker[j].sOverlap = true;
				marker[j].LMin = Min1;
				marker[j].RMin = Min2;
				if (marker[j].LMin == NULL) marker[j].HFLM = marker[j].Line[1] - marker[j].RMin[1];
				else marker[j].HFLM = marker[j].Line[1] 
					- (marker[j].LMin[1] < marker[j].RMin[1] ? marker[j].LMin[1] : marker[j].RMin[1]);
				j++;
			}
		}
		/*if (i==15295) 
		{
			printf("Line=%f, RM=%f, LM=%f, i=%d\n", marker[i].Line[0], marker[i].RMin[0],
				   marker[i].LMin[0], i);
			printf("j=%d, Min[j-1]=%f, Min[j]=%f\n", j, Min[j-1][0], Min[j][0]);
		}*/
	}
	//printf("Nach Schleife\n");
	//printf("j=%d, AnzahlMarker=%d\n", j, AnzahlMarker);
	if (j == AnzahlMarker - 1) 
	{
		marker[j].LMin = Min2;
		marker[j].RMin = NULL;
		marker[j].HFLM = marker[j].Line[1] - Min2[1];
	}
	BRauschen();
	//vR = 5.0 * Rauschen;
	aR = 2.5 * Rauschen;
	for (i=0; i < AnzahlMarker; i++) marker[i].SNR = marker[i].HFLM / aR;
	//for (i=0; i<AnzahlMarker; i++) if (marker[i].sOverlap) 
			//if (marker[i].HFLM < vR && (i > 0 ? marker[i-1].HFLM < aR : true) 
						 //&& (i + 1 < AnzahlMarker ? marker[i+1].HFLM < aR : true))
				//marker[i].sOverlap = false;
	if (LaserFrequency == 0.0) 
	{
		LaserFrequency = LaserLinie->Line[0];
		if ((i = windowTitle().indexOf("=")) != -1) 
			setWindowTitle(windowTitle().left(i+1) + QString::number(LaserFrequency, 'f', 4));
	}
	//printf("vor delete AMarker\n");
	if (AMarker != 0) delete[] AMarker;
	//printf("Nach delete AMarker\n");
	if (MW != 0) MW->loadLineTables();
	emit SpectrumChanged(this);
}

void Spektrum::setFileName(QString Name)
{
	setWindowTitle(Name);
	MDIChild::setFileName(Name);
}

bool Spektrum::readData(QString DateiName)
{
    double x=0, y=0, Max=0, Min=0;
    int pos=0, i, l, j=-1, k=-1;
    QFile Datei( DateiName );
	if (!read(&Datei)) return false;
	DateiName = windowTitle();
	i = DateiName.indexOf(".");
	if (i != -1) setName(DateiName.left(i));
	else setName(DateiName);
    DMarkers = false;
    QTextStream stream(&Datei);
	QString TextDaten = stream.readLine();
    if (TextDaten.left(5) == "type=")
	{
		if ((TextDaten = TextDaten.right(TextDaten.length() - 5)) == "LaserInducedFluorescence") 
			Type = LaserInducedFluorescence;
		else if (TextDaten == "Absorption") Type = Absorption;
		else if (TextDaten == "NormalizedAbsorption") Type = NormalizedAbsorption;
		else if (TextDaten == "ThermalEmission") Type = ThermalEmission;
		else Type = TextDaten.toInt();
		TextDaten = stream.readLine();
	}
	if (TextDaten.left(6) == "laser=") 
	{
		LaserFrequency = TextDaten.right(TextDaten.length() - 6).toDouble();
		TextDaten = "";
	}
	else
	{
		LaserFrequency = 0.0;
		if (DateiName.mid(i - 1, 1) == "m" || TextDaten.indexOf("srt") != -1 ) 
			TextDaten = "";
		else TextDaten += "\n";
	}
    TextDaten +=  stream.readAll();
    Datei.close();
    l = TextDaten.length();
    printf("Daten eingelesen. Dateigröße=%d Bytes\n", l);
    QByteArray String = TextDaten.toLatin1();
    //printf("Daten in Ascii-String umgewandelt.\n");
    if (Daten->GetDSL() > 0) Daten->reinit();
    if (AnzahlMarker > 0)
    {
		delete[] marker;
		marker = 0;
		AnzahlMarker = 0;
    }
	//printf("Nach Daten.reinit\n");
	char c;// test[2]={0,0};
    bool n=false, marked = false;
    for (i=0; i < l; i++)
    {
		c = String[i];
		if (c==0) 
		{
	    	//printf("String zu Ende bei i=%d, l=%d\n", i, l);
	    	break;
		}
		//printf("i=%d\n", i);
		switch (pos)
		{
			case 0:
	    		if ((c >= '0' &&  c <= '9') || c=='-') pos=1;
	    		else break;
			case 1:
	    		if (c >= '0' && c <= '9') 
	    		{
					//test[0]=c;
					//printf("x=%e, c='%s'\n", x, test);
					x=10*x+c-'0';
	   			}
	    		else if (c=='.') j=i;
	    		else if (c=='-') n=true;
	    		else 
	    		{
					pos=2;
					if (n==true)
					{
		    			//printf("1.:x=%e\n",x);
		    			x*=(-1);
		    			n=false;
		    			//printf("2.:x=%e\n",x);
					}
					if (j!=-1)
					{
		    			for (k=j; k<i-1; k++) x*=0.1;
		    			j=-1;
		    			//printf("4.:x=%e\n",x);
					}
	    		}
	    		break;
			case 2:
	    		if ((c >= '0' && c <= '9') || c=='-' || c=='.') pos=3;
	    		else break;
			case 3:
	    		if (c >= '0' && c <= '9') y=10*y+c-'0';
	    		else if (c=='.') j=i;
	    		else if (c=='-') n=true;
	    		else 
	    		{
					pos = 4;
					if (n==true)
					{
		    			y*=(-1);
		    			n=false;
					}
					if (j!=-1)
					{
		    			for (k=j; k<i-1; k++) y*=0.1;
		    			j=-1;
					}
	    		}
	    		if (c != 10) break;
			case 4:
	    		if (c == 'm') 
	    		{
					if (!marked) AnzahlMarker++;
					marked = true;
	    		}
	    		else if (c == 10 || c == 13) 
	    		{
					pos=0;
					Daten->AddValue(x, y, marked);
					//printf("x=%f, y=%f\n", x, y);
					//if (marked) printf("Marked: x=%f\n", x);
					if (Max==0 && y==0) Max=Min=y;
					if (y > Max) Max = y;
					if (y < Min) Min = y;
					x=y=0;
					marked = false;
					//if (Daten.GetDSL()==2) return;
	    		}
	    		break;
		}
    }
    //printf("Daten umgewandelt.\n");
    nDatenS = 1;
	x = (Max - Min) / 20;
    Min -= x;
    Max += x;
    YMin = Min;
    YMax = Max;
    j = Daten->GetDSL();
    //printf("GetDSL()=%d\n", j);
    if (j == 0)
    {
		printf("Fehler beim Einlesen der Datei oder leere Datei.\n");
		return false;
    }
    XMin = Daten->GetValue(0, 0);
    //printf("NachGetValue\n");
    if (AnzahlMarker > 0)
    {
      marker = new Marker[AnzahlMarker];
      //printf("AnzahlMarker=%d\n", AnzahlMarker);
      LaserLinie = marker;
      k = 0;
      for (i=0; i<j; i++) if (Daten->GetMarked(i)) 
      {
		  Daten->SetMarker(i, marker + k);
		  marker[k].HFLM = -1.0;
		  if (marker[k].Line[1] > LaserLinie->Line[1]) LaserLinie = marker + k;
		  //printf("marker[%d].Line[1]=%f\n", k, marker[k].Line[1]);
		  k++;
      }
	  if (LaserFrequency != 0.0 && fabs(LaserLinie->Line[0] - LaserFrequency) > 0.1)
	  {
		  LaserLinie = 0;
		  for (k=0; (k < AnzahlMarker ? LaserFrequency - marker[k].Line[0] > 0.1 : false); k++) ;
		  if (k < AnzahlMarker) 
			  for (LaserLinie = marker + k;
						    (k < AnzahlMarker ? marker[k].Line[0] - LaserFrequency <= 0.1 : false); k++)
				  if (LaserLinie->Line[1] < marker[k].Line[1]) LaserLinie = marker + k;
	  }
      //printf("LaserLinie->Line[1]=%f\n", LaserLinie->Line[0]);
    }
    XMax = Daten->GetValue(j - 1, 0);
    if (XMin > XMax)
    {
		Daten->ReverseOrder();
		x = XMin;
		XMin = XMax;
		XMax = x;
    }
    xStart->setText(QString::number(XMin, 'g', 11));
    xStop->setText(QString::number(XMax, 'g', 11));
    yStart->setText(QString::number(Min));
    yStop->setText(QString::number(Max));
    if (LaserFrequency == 0.0) LaserFrequency = (AnzahlMarker > 0 ? LaserLinie->Line[0] : 0.0);
	emit SpectrumChanged(this);
	Paint();
	Saved();
	return true;
}

void Spektrum::setData(double** Data, int numRows)
{
    if (AnzahlMarker > 0)
	{
		delete[] marker;
		marker = 0;
		AnzahlMarker = 0;
	}
	DiagWindow::setData(Data, numRows);
}

bool Spektrum::writeData(QString FileName)
{
	int r, N = Daten->GetDSL();
	QFile Datei(FileName);
	if (!write(&Datei)) return false;
	QTextStream S(&Datei);
	switch (Type)
	{
		case 0:
			S << "type=LaserInducedFluorescence\n";
			S << "laser=" << QString::number(LaserFrequency, 'f', 6) << "\n";
			break;
		case 1:
			S << "type=Absorption\n";
			break;
		case 2:
			S << "type=NormalizedAbsorption\n";
			break;
		case 3:
			S << "type=ThermalEmission\n";
			break;
		default:
			S << "type=" << QString::number(Type) << "\n";
			break;
	} 
	for (r=0; r<N; r++) 
		S << QString::number(Daten->GetValue(r, 0), 'f', 6) << "\t" 
				<< QString::number(Daten->GetValue(r, 1), 'f', 6) 
				<< (Daten->GetMarked(r) ? "\tm\n" : "\n");
	Saved();
	return true;
}

void Spektrum::setType(int nType)
{
	Type = nType;
}

int Spektrum::getType()
{
	return Type;
}

void Spektrum::IsLI()
{
	int i, v, I, J, Js = 0, n = 0, C=0, veu = 0, Jeu = 0;
	double T = 0.0, d, d1, ****ELU = 0;
	QString Text;
	TermTable *TT;
	for (i=0; i < AnzahlMarker; i++) if (marker[i].Marked && marker[i].DisplayData)
	{
		I = marker[i].Iso;
		J = marker[i].Jss;
		v = marker[i].vss;
		Js = marker[i].Js;
		if (ELU == 0) 
		{
			if (marker[i].LState == 0)
			{
				printf("Spektrum::IsLI() error: marker[%d] is not assigned to a loower state!\n", i);
				continue;
			}
			if ((TT = marker[i].LState->getTermTable()) == 0)
			{
				printf("Spektrum::IsLI() error: for the state %s is no term energy data available\n",
					   marker[i].lState.toLatin1().data());
				continue;
			}
			ELU = TT->getData();
			veu = TT->getMaxv();
			Jeu = TT->getMaxJ();
		}
		T += ELU[C][I][v][J] + marker[i].Line[0];
		n++;
	}
	if (n==0)
	{
		QMessageBox::information( this, "MolSpektAnalysis: Show laser excitation", 
								  "There is no properly assigned progression marked!", 
				  QMessageBox::Ok);
		return;
	}
	T /= n;
	isLI(I, v, J, T, d, ELU[C], veu);
	if (J != Js)
	{
		if (J < Js) 
		{
			if ((i=J+2) <= Jeu) isLI(I, n, i, T, d1, ELU[C], veu);
		}
		else if ((i=J-2) >= 0) isLI(I, n, i, T, d1, ELU[C], veu);
		if (fabs(d1) < fabs(d))
		{
			d = d1;
			v = n;
			J = i;
		}
	}
	Text = "Laser excitation: v''=" + QString::number(v) + ", J''=" + QString::number(J) + ", diff="
			+ QString::number(d, 'g', 7);
	QMessageBox::information(this, "MolSpektAnalysis: Show laser excitation", Text, QMessageBox::Ok);
}

bool Spektrum::isLI(const int &I, int &v, const int &J, double T, double &d, double ***ELU, int veu)
{
	T -= LaserFrequency;
	//printf("isLi: T=%f, LaserFrequency=%f\n", T, LaserFrequency);
	for (v=0; (v <= veu ? ELU[I][v][J] < T : false); v++) ;
	if ((v > 0 ? (v <= veu ? ELU[I][v][J] - T > T - ELU[I][v-1][J] : true) : false)) v--;
	if ((d = fabs(T - ELU[I][v][J])) < LaserFTol) return true;
	return false;
}

void Spektrum::DisplayMarkers()
{
    if (DMarkers) 
    {
		DMarkers = false;
		YMax -= 0.33333333 * (YMax - YMin);
    }
    else 
    {
		DMarkers = true;
		YMax += 0.5 * (YMax - YMin);
		if (marker[0].HFLM == -1.0) editFind();
    }
    Paint();
}

void Spektrum::GetMarker(int &AM, Marker *&M)
{
	AM = AnzahlMarker;
	M = marker;
}

void Spektrum::FitLine(double)
{
    /*int i, o, N, j=0, NP;
    N = Daten.GetDSL();
    double x[100], y[100], sig[100], chisq=-1, alamda = -1;
    double **covar = Create1(100, 100), **alpha = Create1(100, 100);
    double achisq = 0, aalamda = -2, MH = MinPeakHight, LH;
    for (i=0; i<N-1 && LP > Daten.GetValue(i, 0); i++);
    printf("LP=%f, Daten.GetValue(i, 0)=%f\n", LP, Daten.GetValue(i, 0));
    if ((i > 0 ? LP - Daten.GetValue(i - 1, 0) < Daten.GetValue(i, 0) - LP : false)) i -= 1;
    printf("LP=%f, Daten.GetValue(i, 0)=%f\n", LP, Daten.GetValue(i, 0));
    if (fabs(LP - Daten.GetValue(i, 0)) > 0.001) 
    {
	printf("Falsches Spektrum geladen!\n");
	return;
    }
    LP = Daten.GetValue(i, 0);
    LH = Daten.GetValue(i, 1);
    if (MH > 0.7 * LH) MH = 0.7 * LH;
    o = i - 55;
    if (o < 0) o = 0;
    if (o + 110 >= N) o = N - 111;
    if (o < 0) 
    {
	printf("Zu kleiner Datensatz!\n");
	return;
    }
    for (i=o, j=0; j < 100 && i < N; i++)
    {
	if (!Daten.GetMarked(i))
	{
	    x[j] = Daten.GetValue(i, 0);
	    y[j] = Daten.GetValue(i, 1);
	    sig[j++] = Rauschen;
	}
    }
    N = 100;
    for (i=1; i<N-1; i++) if (y[i-1] < y[i] && y[i] > y[i+1]) j+=3;
    NP = j/3;
    double a[j];
    int ia[j];
    j=0;
    for (i=1; i<N-1; i++) if (y[i-1] < y[i] && y[i] > y[i+1] && y[i] > MH) 
    {
	ia[j] = 1;
	a[j++] = y[i];
	ia[j] = 1;
	a[j++] = x[i];
	ia[j] = 1;
	a[j++] = 0.03;
    }
    int z=5;
    while (alamda > aalamda || achisq - chisq > 0.01 || (chisq / achisq > 0.999 ? z-- > 0 : true) 
	|| achisq - chisq < 0)
    {
	aalamda = alamda;
	achisq = chisq;
	mrqmin(x - 1, y - 1, sig - 1, N, a - 1, ia - 1, j, covar, alpha, &chisq, &(fgauss), &alamda);
	printf("alamda=%f; chisq=%f\n", alamda, chisq);
    }
    for (i=1; a[i] < LP && i < j - 3; i+=3);
    if (i>2) if (LP - a[i-3] < a[i] - LP) i-=3;
    //Line->Line[0] = a[i];
    //Line->Line[1] = a[i-1];
//	FL[i].Width = a[i+1];
    printf("position=%f, height=%f\n", a[i], a[i-1]);
    Destroy1(covar, N);
    Destroy1(alpha, N);*/
}

void Spektrum::markRegion(QRect *i_regionToMark)
{
    if (!ZoomB->isChecked())
    {
        double dXSF = 1.0 / XSF;
        m_minSelectedFrequency = dXSF * (static_cast<double>(i_regionToMark->left()) - XO);
        m_maxSelectedFrequency = dXSF * (static_cast<double>(i_regionToMark->right()) - XO);
        Paint();
    }
}

void Spektrum::FitGaussianLineProfile(double &o_energy, double &o_intensity, double &o_width, double &o_intensityOffset, double &o_sigma)
{
    if (Rauschen == -1.0) editFind();
    double sig = 1.0 / (Rauschen * Rauschen), cx = m_minSelectedFrequency;
    int nStart, N = 0;
    for (int n=0; n < Daten->GetDSL() && cx < m_maxSelectedFrequency; ++n) if ((cx = Daten->GetValue(n, 0)) >= m_minSelectedFrequency)
    {
        if (N==0) nStart = n;
        if (cx <= m_maxSelectedFrequency) ++N;
    }
    if (N==0)
    {
        o_sigma = -1.0;
        return;
    }
    double *X = new double[N], *Y = new double[N], *Sig = new double[N];
    int n, m;
    for (n = nStart, m=0; m < N; ++n, ++m)
    {
        X[m] = Daten->GetValue(n, 0);
        Y[m] = Daten->GetValue(n, 1);
        Sig[m] = sig;
    }
    Gaussian* line = new Gaussian(X, Y, Sig, N);
    double chiSq = line->LevenbergMarquardt(100, 0.01);
    o_sigma = sqrt(chiSq / (N-1));
    line->GetValues(o_intensity, o_energy, o_width, o_intensityOffset);
    m_fitedLineVector.push_back(line);
    Paint();
}

void Spektrum::PictureClicked(QPoint *P)
{
    printf("Spektrum::PictureClicked\n");
	DiagWindow::PictureClicked(P);
    if (m_minSelectedFrequency > 0.0 || m_maxSelectedFrequency > 0.0)
    {
        m_minSelectedFrequency = m_maxSelectedFrequency = -1.0;
        Paint();
    }
    if (!DMarkers) return;
    int i, j, d, J, v, vs, Jss, Js, vss, Iso, x=P->x(), y=P->y(), C=-1, xd, yd, veu, Jeu, nC, nI, veo;
    int JStep, JStart, Jeo, c, NL, Jd, J1, J2 = 0, J3, JS, FC;
	double CPos, nd, ***ELU, ****UT, T=0.0, sdiff, B, D=0.0, O=0.0, T1, T2 = 0.0, dist = 0.0, X, Y;
	bool used, satellites = false, DispData;
	TermTable *TT = 0;
	ElState *UState, *LState;
	Molecule *Mol;
	IsoTab *IsoT = 0;
	for (i=0; i<AnzahlMarker; i++) if (marker[i].Visible && x >= marker[i].x1 && x <= marker[i].x2
				       && y >= marker[i].y2 && y <= marker[i].y1)
    {
		//printf("marker[%d] clicked\n", i);
		//if (marker[i].sOverlap) printf("marker[%d].sOverlap=true\n", i);
		if (marker[i].Marked) 
		{
			marker[i].Marked = false;
			marker[i].DisplayData = false;
			C=-2;
			LineTable *L = (marker[i].Mol != 0 ? 
					marker[i].Mol->getLineTable(marker[i].LState, marker[i].UState) :
					(MW != 0 ? MW->getUnassignedLineTable() : 0));
			if (L != 0) L->AddMarked(AnzahlMarker, marker, LaserLinie, getFileName());
		}
		else C=i;
		//printf("Ende Block\n");
	}
	if (C==-1) for (i=0; i < AnzahlMarker; i++) 
			if (x >= marker[i].x1 && x <= marker[i].x2 && y >= marker[i].y2 && y <= marker[i].y1)
	{
		xd = x - (marker[i].x2 - marker[i].x1) / 2;
		yd = y - (marker[i].y1 - marker[i].y2) / 2;
		if ((nd = sqrt(xd * xd + yd * yd)) < dist)
		{
			dist = nd;
			C=i;
		}
	}
	if (C >= 0 ? marker[C].Marked : false)
	{
		marker[C].Marked = false;
		marker[C].DisplayData = false;
		C=-2;
	}
	if (C >= 0)
	{
		marker[C].Marked = true;
		if (AssignmentStatus != 1) marker[C].DisplayData = false;
		marker[C].added = false;
		if (AssignmentStatus == 0)
		{
			for (j=0; j < AnzahlMarker - 1 && !(marker[j].Marked && marker[j].DisplayData 
					  && (marker[j].LState != 0 ? ((TT = marker[j].LState->getTermTable()) != 0 ? 
						  marker[j].vss <= (veu = TT->getMaxv()) 
					      && marker[j].Jss <= (Jeu = TT->getMaxJ()) : false) : false)); j++) ;
			//printf("Nach 1M\n");
			if (marker[j].Marked && marker[j].DisplayData)
			{
				if ((ELU = TT->getData()[0]) != 0) 
				{
					Iso = marker[j].Iso;
					vs = marker[j].vs;
					Js = marker[j].Js;
					FC = marker[j].FC;
					LState = marker[j].LState;
					UState = marker[j].UState;
					Mol = marker[j].Mol;
					CPos = ELU[Iso][marker[j].vss][marker[j].Jss] + marker[j].DD + marker[j].Line[0] 
							- marker[C].Line[0];
					if (marker[j].Jss == Js) Jss = Js;
					else Jss = -1;
					for (v=0; (v <= veu ? ELU[Iso][v][Js] < CPos : false); v++) ;
					if (v > veu) v--;
					if (v > 0) if (CPos - ELU[Iso][v-1][Js] < ELU[Iso][v][Js] - CPos) v--; 
					if (Jss == -1)
					{
						if ((Js > 0 ? fabs(CPos - ELU[Iso][v][Js - 1]) : 1000 ) 
										< (Js < Jeu ? fabs(CPos - ELU[Iso][v][Js + 1]) : 1000))
							Jss = Js - 1;
						else Jss = Js + 1;
					}
					for (j=0; j < AnzahlMarker - 1 && !(marker[j].Marked && marker[j].DisplayData 
								 && marker[j].Jss == Jss && marker[j].vss == v); j++) ;
					if (marker[j].Marked && marker[j].DisplayData && marker[j].Jss == Jss 
								&& marker[j].vss == v)
					{
						//printf("Best level v=%d, J=%d already used!\n", v, Jss);
						if (marker[j].Line[0] > marker[C].Line[0]) d = 1;
						else d = -1;
						J = 2 * Js - Jss;
						for (j=0; j < AnzahlMarker - 1 && !(marker[j].Marked && marker[j].DisplayData 
								 && marker[j].Jss == J && marker[j].vss == v); j++) ;
						if ((marker[j].Marked && marker[j].DisplayData && marker[j].Jss == J 
										&& marker[j].vss == v) || J < 0)
						{	
							used = true;
							while (used)
							{
								v += d;
								for (j=0; j < AnzahlMarker - 1 && !(marker[j].Marked 
									 && marker[j].DisplayData 
							 		&& marker[j].Jss == Jss && marker[j].vss == v); j++) ;
								if (!marker[j].Marked || !marker[j].DisplayData 
									|| marker[j].Jss != Jss || marker[j].vss != v)
									used = false;
								if (J < 0) continue;
								for (j=0; j < AnzahlMarker - 1 && !(marker[j].Marked 
									 && marker[j].DisplayData 
							 		&& marker[j].Jss == J && marker[j].vss == v); j++) ;
								if (!marker[j].Marked || !marker[j].DisplayData || marker[j].Jss != J 
												  || marker[j].vss != v)
								{
									if (used) Jss = J;
									else
									{
										used = false;
										if ((d > 0 && J < Jss) || (d < 0 && J > Jss)) Jss = J;
									}
								}
							}
						}
						else if (J >= 0) Jss = J;
					}
					if (v >= 0)
					{
						//printf("Vor 2M\n");
						marker[C].DisplayData = true;
						marker[C].vs = vs;
						marker[C].Js = Js;
						marker[C].vss = v;
						marker[C].FC = FC;
						marker[C].Jss = Jss;
						marker[C].UState = UState;
						marker[C].LState = LState;
						marker[C].lState = LState->getName();
						if (UState != 0) marker[C].uState = UState->getName();
						marker[C].Mol = Mol;
						if (v <= veu && Jss <= Jeu) marker[C].DD = CPos - ELU[Iso][v][Jss];
						else marker[C].DD = 0.0;
						marker[C].Iso = Iso;
						IsoT = Mol->getIso();
						marker[C].IsoName = IsoT->texName[Iso];
						//printf("marker[%d] marked\n", C);
					}
				}
			}
		}
		else if (AssignmentStatus == 1)
		{
			//printf("Assign bands\n");
			DispData = marker[C].DisplayData;
			marker[C].DisplayData = true;
			for (i=C+1; (i < AnzahlMarker ? (!marker[i].Marked || !marker[i].DisplayData) 
							 && marker[i].Line[0] - marker[C].Line[0] < 1e2 : false); i++) ;
			for (j=C-1; (j>=0 ? (!marker[j].Marked || !marker[j].DisplayData) 
							 && marker[C].Line[0] - marker[j].Line[0] < 1e2 : false); j--) ;
			if (2 * marker[C].Line[0] < marker[j].Line[0] + marker[i].Line[0]) i=j;
			if (marker[i].Marked && marker[i].DisplayData)
			{
				//printf("Marked\n");
				Jd = (Jss = marker[i].Jss) - (Js = marker[i].Js);
				marker[C].overlap = (marker[C].sOverlap || (DispData && (marker[C].Mol != marker[i].Mol 
									|| marker[C].Iso != marker[i].Iso || marker[C].FC != marker[i].FC 
									|| marker[C].vs != marker[i].vs || marker[C].vss != marker[i].vss)) ? true : false);
				marker[C].uncertainty = (marker[C].overlap ? 0.02 : 0.005);
				Mol = marker[i].Mol;
				Iso = marker[i].Iso;
				FC = marker[i].FC;
				vss = marker[i].vss;
				vs = marker[i].vs;
				if (!DispData || marker[C].Iso != Iso || marker[C].FC != FC || marker[C].vss != vss || marker[C].vs != vs
					|| marker[C].Jss - marker[C].Js != Jd || marker[C].lState != marker[i].lState
					|| marker[C].uState != marker[i].uState)
				{
					marker[C].Js = -1;
					marker[C].vs = vs;
					marker[C].FC = FC;
					marker[C].vss = vss;
					marker[C].UState = marker[i].UState;
					marker[C].LState = marker[i].LState;
					marker[C].Mol = Mol;
					marker[C].Iso = Iso;
					marker[C].IsoName = marker[i].IsoName;
					marker[C].lState = marker[i].lState;
					marker[C].uState = marker[i].uState;
				}
				if (marker[C].HFLM * 5.0 < marker[i].HFLM || marker[i].satellite) satellites = true;
				JStep = (Mol != 0 ? Mol->getJStep(Iso) : 1);
				if (JStep == 2 && 2 * (Jss / JStep) != Jss) JStart = 1;
				else JStart = 0;
				if (JStart - Jd == -1) JStart += JStep;
				//printf("Vor band\n");
				Band band;
				for (J = 0; J < 1000; J++) band.lines[J] = 0;
				for (j = NL = 0; j < AnzahlMarker; j++) 
					if (marker[j].Marked && marker[j].DisplayData && marker[j].Iso == Iso
						 && marker[j].vss == vss && marker[j].vs == vs && marker[j].FC == FC
						 && marker[j].Jss - marker[j].Js == Jd && marker[j].Js >= 0 && j != C
						 && marker[j].Js < 1000 && marker[j].lState == marker[i].lState && marker[j].uState == marker[i].uState)
				{
					band.lines[marker[j].Js] = marker + j;
					NL++;
					//printf("AdToBand: j=%d, J=%d\n", j, marker[j].Js);
				}
				if (NL == 2)
				{
					for (J1 = 0; band.lines[J1] == 0; J1++) ;
					for (J2 = J1 + 1; band.lines[J2] == 0; J2++) ;
					D = double(band.lines[J2]->Js - band.lines[J1]->Js)
					  /	(band.lines[J2]->Line[0] - band.lines[J1]->Line[0]); 
					marker[C].Js = band.lines[J2]->Js 
								 + rund((marker[C].Line[0] - band.lines[J2]->Line[0]) * D);
					marker[C].Jss = marker[C].Js - Jd;
					marker[C].oc = marker[C].Line[0] - band.lines[J2]->Line[0] 
							     + double(band.lines[J2]->Js - marker[C].Js) / D;
				}
				else if (NL >= 3)
				{
					//printf("Parab\n");
					marker[C].oc = 1e3;
					for (J1 = J2 = -1, J3 = JStart, j = 0; j < NL; J3 += JStep) if (band.lines[J3] != 0)
					{
						if (++j >= 3 && (J3 - J2 > JStep || j==3 || j == NL))
						{
							if (j==3) JS = JStart;
							else JS = J2 + JStep;
							if (j == NL) Jss = 999;
							else Jss = J3;
							ParabInterpol(J1, band.lines[J1]->Line[0], J2, band.lines[J2]->Line[0],
										  J3, band.lines[J3]->Line[0], X, Y);
							B = (J1 != X ? (band.lines[J1]->Line[0] - Y) / ((J1 - X) * (J1 - X))
										 : (band.lines[J2]->Line[0] - Y) / ((J2 - X) * (J2 - X)));
							D = Y + (JS - X) * (JS - X) * B;
							O = Y + (Jss - X) * (Jss - X) * B;
							T = fabs(marker[C].oc);
							//printf("X=%f, Y=%f, B=%f, D=%f, O=%f, T=%f\n", X, Y, B, D, O, T);
							//printf("J1=%d, J2=%d, J3=%d, JS=%d, Jss=%d\n", J1, J2, J3, JS, Jss);
							if (JS <= X && X <= Jss ? 
								((J1 != X ? band.lines[J1]->Line[0] > Y : band.lines[J2]->Line[0] > Y) ?
								 marker[C].Line[0] + T > D || marker[C].Line[0] + T > O 
								 : marker[C].Line[0] - T < D || marker[C].Line[0] - T < O)
								: (D < marker[C].Line[0] + T && marker[C].Line[0] - T < O) 
									|| (O < marker[C].Line[0] + T && marker[C].Line[0] - T < D)) 
								for (J = JS; J <= Jss; J += JStep) 
							{
								if (band.lines[J] != 0) continue;
								//printf("sdiff=%f\n", sdiff);
								if (fabs(sdiff = marker[C].Line[0] - Y+(X-J)*(J-X)*B) < T)
								{
									//printf("sdiff=%f, J=%d\n", sdiff, J);
									T = fabs(marker[C].oc = sdiff);
									marker[C].Js = J;
								}
								else if (j == NL && J > J3) break;
							}
						}
						J1 = J2;
						J2 = J3;
					}
				}
				if (marker[i].Mol != 0)
				{
					//printf("Begin Assign J, marker[C].Js=%d, marker[C].oc=%f\n", 
						//   marker[C].Js, marker[C].oc);
					IsoT = Mol->getIso();
					Mol->getTermData(0, nC, nI, veu, Jeu, UT);
					if (UT != 0)
					{
						//printf("Js=%d, Jd=%d, Jeu=%d, marker[C].oc=%f\n", 
							//marker[C].Js, Jd, Jeu, marker[C].oc);
						ELU = UT[0];
						UT = 0;
						if ((vss >= 0 && vss < veu && marker[C].Js + Jd < Jeu && marker[C].Js + Jd >= 0 ?
								ELU[Iso][vss][marker[C].Js + Jd] != 0.0 : false) || marker[C].Js == -1) 
							marker[C].oc = 10.0;
						if (marker[C].UState != 0) 
							marker[C].Mol->getTermData(marker[C].UState->getStateNum(), 
								nC, nI, veo, Jeo, UT);
						if (vss >= 0 && vss < veu)
						{
							//printf("TermData Loaded, Jeu=%d, JStep=%d, JStart=%d\n", Jeu, JStep, JStart);
							c = (nC > 1 && Jd == 0 ? 1 : 0);
							for (J = JStart, J1 = -1; J < Jeu; J += JStep) 
								if ((band.lines[J - Jd] != 0 && ELU[Iso][vss][J] != 0) || J >= Jeu - JStep)
							{
								//printf("J=%d\n", J);
								if (J2 == -1) Jss = JStart;
								else Jss = J1 + JStep;
								if (J < Jeu - JStep)
								{
									J2 = J1;
									J1 = J;
									if (J2 == -1) continue;
								}
								//printf("J1=%d, J2=%d, Iso=%d, vss=%d, Jss=%d\n", J1, J2, Iso, vss, Jss);
								T1 = ELU[Iso][vss][J1] + band.lines[J1 - Jd]->Line[0];
								if (J2 != -1)
								{
									T2 = ELU[Iso][vss][J2] + band.lines[J2 - Jd]->Line[0];
									B = (T2 - T1) / double((J2 - Jd) * (J2 - Jd + 1) 
											- (J1 - Jd) * (J1 - Jd + 1));
									T = T1 - B * double((J1 - Jd) * (J1 - Jd + 1));
								}
								else 
								{
									B = 0.0;
									if (UT == 0) break;
								}
								//printf("T1=%f, T2=%f, T=%f, B=%f\n", T1, T2, T, B);
								if (UT != 0 && vs < veo && vs >= 0 && J1 < Jeo? 
														UT[c][Iso][vs][J1 - Jd] != 0 : false) 
									T1 -= UT[c][Iso][vs][J1 - Jd];
								else T1 = 0.0;
								if (UT != 0 && J2 != -1 && vs < veu && vs >= 0 && J2 < Jeo) 
									T2 -= UT[c][Iso][vs][J2 - Jd];
								else T2 = 0.0;
								if (T1 != 0.0)
								{
									if (T2 != 0.0)
									{
										D = (T2 - T1) / double((J2 - Jd) * (J2 - Jd + 1) 
												- (J1 - Jd) * (J1 - Jd + 1));
										O = T1 - D * double((J1 - Jd) * (J1 - Jd + 1));
									}
									else
									{
										D = 0.0;
										O = T1;
									}
								}
								//printf("T1=%f, T2=%f, O=%f, D=%f\n", T1, T2, O, D);
								//printf("2.\n");
								for (j = Jss; j <= J; j += JStep)
								{
									if (ELU[Iso][vss][j] == 0.0) 
									{
										if (j == 0) continue;
										else break;
									}
									if (UT != 0 && j - Jd < Jeo && vs < veo && vs >= 0? 
															UT[c][Iso][vs][j - Jd] != 0.0 : false)
										sdiff = ELU[Iso][vss][j] + marker[C].Line[0] 
											- (UT[c][Iso][vs][j - Jd] + O + (j - Jd) * (j - Jd + 1) * D);
									else sdiff = ELU[Iso][vss][j] + marker[C].Line[0] 
												- (T + (j - Jd) * (j - Jd + 1) * B);
									if (fabs(marker[C].oc) > fabs(sdiff) && band.lines[j - Jd] == 0)
									{
										marker[C].oc = sdiff;
										marker[C].Js = j - Jd;
									}
									//printf("J=%d, sdiff=%f\n", j, sdiff);
								}
							}
							//printf("Iso=%d, vss=%d, vs=%d, c=%d\n", Iso, vss, vs, c);
							if (UT != 0) //{
								marker[C].oc = ELU[Iso][vss][marker[C].Js + Jd] + marker[C].Line[0]
											- UT[c][Iso][vs][marker[C].Js];
								//printf("oc=%f, ELU=%f, Line=%f, UT=%f\n", marker[C].oc, 
									//   ELU[Iso][vss][marker[C].Jss], marker[C].Line[0],
										//UT[c][Iso][vs][marker[C].Js]);}
							else marker[C].oc = 0.0;
						}
					}
				}
				if (marker[C].Js == -1)
				{
					if (marker[C].Line[0] < marker[i].Line[0]) marker[C].Js = marker[i].Js - 1;
					else marker[C].Js = marker[i].Js + 1;
				}
				marker[C].Jss = marker[C].Js + Jd;
				marker[C].satellite = satellites;
				if (marker[C].Js >= 0)
				{
					for (J = marker[C].Js; (J > 0 ? band.lines[J] == 0 : false); J--) ;
					if (band.lines[J] != 0) {printf("p: J=%d\n", J); 
						marker[C].DD = (marker[C].Line[0] - band.lines[J]->Line[0]) 
									 / (marker[C].Js - J);}
					for (J = marker[C].Js; (J < 999 ? band.lines[J] == 0 : false); J++) ;
					if (band.lines[J] != 0) {printf("n: J=%d\n", J);
						band.lines[J]->DD = (band.lines[J]->Line[0] - marker[C].Line[0])
								          / (J - marker[C].Js);}
				}
			}
			else if (!DispData)
			{
				marker[C].Js = 50;
				marker[C].Mol = 0;
			}
		}
		if (!marker[C].DisplayData)
		{
			//printf("Vor DCalc\n");
			for (j=C-1; (j>=0 ? !marker[j].Marked : false); j--) ;
			if (j>=0) marker[j].DD = marker[C].Line[0] - marker[j].Line[0];
			for (j=C+1; (j < AnzahlMarker ? !marker[j].Marked : false); j++) ;
			if (j < AnzahlMarker) marker[C].DD = marker[j].Line[0] - marker[C].Line[0];
			else marker[C].DD = -1.0;
			marker[C].LState = 0;
			marker[C].UState = 0;
			marker[C].Mol = 0;
			marker[C].added = false;
			//printf("Nach DCalc\n");
		}
		addMarker(false);
	}
	if (IsoT != 0) delete IsoT;
    Paint();
}

void Spektrum::Paint()
{
	if (AnzahlMarker > 0) emit SpectrumChanged(this);
	DiagWindow::Paint();
}

void Spektrum::PSpektrum(QPainter &P, const QRect &A, bool PrintFN)
{
    DiagWindow::PSpektrum(P, A, PrintFN);
    P.setPen(QColor(0, 0, 255));
    int xm = ScaleYWidth + 2, xM = A.width() - 2, x, yc, yl;
    int ym = ScaleTopOffset + 2, yM = A.height() - ScaleXHeight - 2;
    double Estart, Eend, dXSF = 1.0 / XSF;
    for (size_t n = 0; n < m_fitedLineVector.size(); ++n)
    {
        Gaussian* line = m_fitedLineVector[n];
        line->GetDataRange(Estart, Eend);
        int xStart = ((x = int(XO + XSF * Estart)) > xm ? x : xm);
        int xEnd = ((x = int (XO + XSF * Eend)) < xM ? x : xM);
        if (xStart > xM || xEnd < xm) continue;
        for (x = xStart; x < xEnd; ++x)
        {
           yc = static_cast<int>(YO - YSF * line->GetPoint(dXSF * static_cast<double>(x - XO)));
           if (x > xStart && ((yc >= ym && yc <= yM) || (yl >= ym && yl <= yM)))
           {
               if (yc < ym) yc = ym;
               else if (yc > yM) yc = yM;;
               P.drawLine(x-1, yl, x, yc);
           }
           yl = yc;
        }
    }
}

int Spektrum::PQC()
{
	/*int n, i=0, j, v, J, N=0, I;
	for (n=0; n < AnzahlMarker; n++) if (marker[n].Marked && marker[n].DisplayData) N++;
	if (N < 5) return -1;
	Marker *mB[N], *mP;
	double E = 0.0, r = 0.5 * Rauschen, h[6];
	for (n=0; n < AnzahlMarker; n++) 
		if (marker[n].Marked && marker[n].DisplayData) mB[i++] = marker + n;
	J = mB[0]->Js;
	I = mB[0]->Iso;
	if (J==0 || J==Jeu || AnzahlMarker < 5) return 0;
	for (n=0; n<N; n++) 
	{
		if (mB[n]->DD > ST) return -2;
		E += ELU[I][mB[n]->vss][mB[n]->Jss] + mB[n]->Line[0];
		for (i=0; i<n; i++) if (mB[i]->vss == mB[n]->vss && mB[i]->Jss == mB[n]->Jss) return -3;
	}
	for (n=0, i=N; n<N; n++)
	{
		mP = (mB[n] >= marker + 2 ? (mB[n] < marker + AnzahlMarker - 2 ? mB[n] - 2 
			: marker + AnzahlMarker - 5) : marker);
		for (j=0; j<6; j++) h[j] = mP[j].Line[1];
		if ((h[0] < h[1] && h[1] < h[2] && h[2] < h[3] && h[3] < h[4] && h[4] < h[5]) ||
			(h[0] > h[1] && h[1] > h[2] && h[2] > h[3] && h[3] > h[4] && h[4] > h[5]) ||
			(fabs(h[0] - h[1]) < r && fabs(h[1] - h[2]) < r && fabs(h[2] - h[3]) < r 
				   && fabs(h[3] - h[4]) < r && fabs(h[4] - h[5]) < r)) 
		{
			mB[n]->Marked = false;
			i--;
		}
		else for (j=0; j<n; j++) 
				if (mB[n]->vss > mB[j]->vss + 1 || (mB[n]->Jss == mB[j]->Jss && mB[n]->vss > mB[j]->vss))
		{
			mB[n]->Marked = false;
			i--;
		}
	}
	if (i < 5) return -7;
	if (i < N) 
	{
		for (j=0; mB[j]->Marked; j++);
		for (n=j+1; n<N; n++) if (mB[n]->Marked) mB[j++] = mB[n];
		N = i;
	}
	E /= N;
	E -= LaserFrequency;
	if (J == mB[0]->Jss) 
	{
		if (OnlyDoublet) return -8;
		for (v=1; (v<=veu ? ELU[I][v][J] < E && ELU[I][v][J] > ELU[I][v-1][J] : false); v++);
		if (E - ELU[I][v-1][J] > 0.1 && (v <= veu ? ELU[I][v][J] - E > 0.1 : false)) return -4;
		printf ("v=%d\n", v);
		return 1;
	}
	for (n=0, i=1; i<N; i++) for (j=0; j<i; j++) if (mB[i]->vss == mB[j]->vss) 
	{
		r = mB[i]->HFLM / mB[j]->HFLM;
		if (r > 0.8 && r < 1.25) n++;
	}
	if (n < 3) return -5;
	for (j=J-1; j <= J+1; j+=2) 
	{
		for (v=1; (v<=veu ? ELU[I][v][j] < E && ELU[I][v][j] > ELU[I][v-1][j] : false); v++);
		if (E - ELU[I][v-1][j] < 0.1 || (v <= veu ? ELU[I][v][j] - E < 0.1 : true)) return 2;
	}*/
	return -6;
}

/*IntProg* Spektrum::SearchSLUProg(Marker *M1, Marker *M2)
{
	if (veu==0) return 0;
	if (LaserFrequency == 0.0) LaserFrequency = LaserLinie->Line[0];
	int n, I, J, Js, v, vR, JR, vLQ, vLR, vS, vSR;
	double D, DL, SPos, STE, LS;
	IntProg *AProg = new IntProg, *BProg = 0;
	Marker *MB;
	if (M1->Line[0] > M2->Line[0])
	{
		MB = M1;
		M2 = M2;
		M2 = MB;
	}
	D = M2->Line[0] - M1->Line[0];
	DL = M2->Line[0] - LaserFrequency;
	LS = M1->Line[0] + M2->Line[0];
	for (I=0; I<NI; I++) for (J=0; J<=Jeu; J++) 
	{
		JR = J + 2;
		if (JR <= Jeu) for (v=0, vR=0, vLQ=0, vLR=0; v<=veu; v++)
		{
			SPos = ELU[I][v][J] + D;
			while ((vR <= veu ? ELU[I][vR][JR] < SPos : false)) vR++;
			if ((vR > 0 ? (vR <= veu ? ELU[I][vR][JR] - SPos > SPos - ELU[I][vR-1][JR] : true) : false)) vR--;
			if (fabs(SPos - ELU[I][vR][JR]) > ST) continue;
			SPos = ELU[I][v][J] + DL;
			while ((vLQ <= veu ? ELU[I][vLQ][J] < SPos : false)) vLQ++;
			if ((vLQ > 0 ? (vLQ <= veu ? ELU[I][vLQ][J] - SPos > SPos - ELU[I][vLQ-1][J] : true) : false)) vLQ--;
			while ((vLR <= veu ? ELU[I][vLR][JR] < SPos : false)) vLR++;
			if ((vLR > 0 ? (vLR <= veu ? ELU[I][vLR][JR] - SPos > SPos - ELU[I][vLR-1][JR] : true) : false)) vLR--;
			if (fabs(SPos - ELU[I][vLQ][J]) > 0.1 && fabs(SPos - ELU[I][vLR][JR]) > 0.1) continue;
			AProg->N = AProg->G = 0;
			AProg->Iso = I;
			AProg->FQS = 0.0;
			Js = J + 1;
			STE = 0.5 * (LS + ELU[I][v][J] + ELU[I][vR][JR]);
			for (n=0, vS=0, vSR=0; vS < veu && n < AnzahlMarker; n++)
			{
				SPos = STE - marker[n].Line[0];
				while ((vS <= veu ? ELU[I][vS][J] < SPos : false)) vS++;
				if ((vS > 0 ? (vS <= veu ? ELU[I][vS][J] - SPos > ELU[I][vS-1][J] : true) : false)) vS--;
								 				
			}
		}
	}
}*/

void Spektrum::SetLaserFrequency(double Frequency)
{
	int n;
	LaserFrequency = Frequency;
	LaserLinie = NULL;
	for (n=0; n<AnzahlMarker; n++) 
		if (fabs(marker[n].Line[0] - LaserFrequency) < 0.1 
				  && (LaserLinie != NULL ? marker[n].Line[1] > LaserLinie->Line[1] : true)) 
			LaserLinie = marker + n;
	if ((n = windowTitle().indexOf("=")) != -1) 
		setWindowTitle(windowTitle().left(n+1) + QString::number(LaserFrequency, 'f', 4));
}

double Spektrum::GetLaserFrequency()
{
	return LaserFrequency;
}

Marker *Spektrum::GetLaserLine()
{
	return LaserLinie;
}

void Spektrum::FoundB()
{
    ClearMarked();
    Prog.RotateB();
    ShowFound();
}

void Spektrum::FoundF()
{
    ClearMarked();
    Prog.RotateF();
    ShowFound();
}

void Spektrum::GetMA(Marker **&A, int &N)
{
	int n;
	N = 0;
	for (n=0; n < AnzahlMarker; n++) if (marker[n].DisplayData) N++;
	if (A!=NULL) delete[] A;
	A = new Marker*[N];
	for (n=0, N=0; n < AnzahlMarker; n++) if (marker[n].DisplayData) A[N++] = marker + n;
}

void Spektrum::GetMGMH()
{
	printf("GetMGMH\n");
	int n, N;
	if (GMarker != NULL) delete[] GMarker;
	for (N=0, n=0; n < AnzahlMarker; n++) if (marker[n].HFLM >= MinPeakHeight) N++;
	GMarker = new Marker*[N+2];
	GMarker[0] = &MBeginn;
	//printf("&MStop=%d\n", &MStop);
	GMarker[N+1] = &MStop;
	//printf("GMarker[%d]=%d\n", N+1, GMarker[N+1]);
	for (n=0, N=1; n < AnzahlMarker; n++) 
		if (marker[n].HFLM >= MinPeakHeight) GMarker[N++] = marker + n;
	//printf("GMarker[%d]=%d\n", N, GMarker[N]);
}

void Spektrum::ShowFound()
{
    int N;
    Marker **M;
    Prog.GetMarker(M, N);
    FokusProg(M, N);
}
	
void Spektrum::FokusProg(Marker **M, int &N)
{
	double Min=100000, Max=0, R, yMax = 0;
	int i;
	if (N==0) return;
    for (i=0; i<N; i++)
    {
		if (M[i]->Line[0] < Min) Min = M[i]->Line[0];
		if (M[i]->Line[0] > Max) Max = M[i]->Line[0];
		if (M[i]->Line[1] > yMax) yMax = M[i]->Line[1];
    }
    R = (Max - Min) * 0.1;
    Min -= R;
    Max += R;
	R = yMax - YMin;
	yMax += R;
    if (YMax < 2 * yMax) YMax = 2 * yMax;
    xStart->setText(QString::number(Min, 'g', 11));
    xStop->setText(QString::number(Max, 'g', 11));
    yStart->setText(QString::number(YMin, 'g', 5));
    yStop->setText(QString::number(yMax, 'g', 5));
	if (isVisible()) 
	{
		addMarker(false);
		if (MW != 0) MW->setActive(this);
    	Paint();
	}
}

void Spektrum::FocusSelected()
{
	double Min=100000, Max=0, R, yMax = 0;
	int i;
	if (AnzahlMarker == 0) return;
    for (i=0; i < AnzahlMarker; i++) if (marker[i].Marked)
    {
		if (marker[i].Line[0] < Min) Min = marker[i].Line[0];
		if (marker[i].Line[0] > Max) Max = marker[i].Line[0];
		if (marker[i].Line[1] > yMax) yMax = marker[i].Line[1];
    }
    R = (Max - Min) * 0.1;
    Min -= R;
    Max += R;
	R = yMax - YMin;
	yMax += R;
    if (YMax < 2 * yMax) YMax = 2 * yMax;
    xStart->setText(QString::number(Min, 'g', 11));
    xStop->setText(QString::number(Max, 'g', 11));
    yStart->setText(QString::number(YMin, 'g', 5));
    yStop->setText(QString::number(yMax, 'g', 5));
	addMarker(false);
    Paint();
}

void Spektrum::CSLSettings()
{
    QDialog *D = new QDialog;
	QGridLayout *L = new QGridLayout(D);
	D->setWindowTitle("Settings for assignments");
	QCheckBox *B = new QCheckBox("Use intensities", D);
	B->setChecked(useIntensities);
	L->addWidget(B, 0, 0, 1, 2);
	L->addWidget(new QLabel("Min peak hight:", D), 1, 0);
	QLineEdit *MP = new QLineEdit(QString::number(MinPeakHeight), D), *STE = new QLineEdit(QString::number(ST), D);
	L->addWidget(MP, 1, 1);
	L->addWidget(new QLabel("Freq. tolerancy:", D), 2, 0);
	L->addWidget(STE, 2, 1);
	L->setRowMinimumHeight(3, 20);
	QPushButton *OK = new QPushButton("OK", D), *Cancel = new QPushButton("Cancel", D);
	L->addWidget(OK, 4, 0);
	L->addWidget(Cancel, 4, 1);
	connect(OK, SIGNAL(clicked()), D, SLOT(accept()));
	connect(Cancel, SIGNAL(clicked()), D, SLOT(reject()));
	//InputDialog D(MinPeakHight, ST);
    if (D->exec() == QDialog::Rejected) return;
    MinPeakHeight = MP->text().toDouble();
	ST = STE->text().toDouble();
	useIntensities = B->isChecked();
	//MinPeakHight = D.GetMinPeak();
    //ST = D.GetToleranz();
	GetMGMH();
}

void Spektrum::ClearMarked()
{
    printf("Spektrum::ClearMarked(): numLineTablesAT=%d\n", numLineTablesAT);
	int i;
    for (i=0; i<AnzahlMarker; i++) 
    {
		marker[i].Marked = false;
		marker[i].DisplayData = false;
    }
	QString File = getFileName();
	for (i=0; i < numLineTablesAT; i++) 
	{
		LineTablesAT[i]->AcceptAssignments(File, false);
		if (LineTablesAT[i]->getAnzahlLinien() == 0) LineTablesAT[i]->close();
	}
	numLineTablesAT = 0;
    //if (tabelle != NULL) tabelle->EndMarkSelected();
	printf("Vor Paint()\n");
    Paint();
	printf("Ende ClearMarked\n");
}


void Spektrum::TTransition(int NumWFPoints)
{
    if (AnzahlMarker == 0) 
    {
		QMessageBox::information( this, "MolSpektAnalysis", "Es sind keine Marker geladen!", 
				  QMessageBox::Ok);
		return;
    }
    if (MW == 0)
    {
		QMessageBox::information(this, "MolSpektanalysis", 
				 "Es sind keine Moleküldaten geladen!", QMessageBox::Ok);
		return;
    }
    int I, v, J, i, j = 0, N = 0, g, k, gJ, kJ, l, vs, NI, mveu = 0, veu, Jeu, veo, Jeo, M, NC, NS;
	int s, c;
	int numMolecules = MW->getNumMolecules(), lvs;
    double SLines = 0.0, SLinesg, SLinesk, SMarker = 0.0, SBDk, SBDg, ****ELl, ****ELu, ocs, locs;
    Molecule *Mol;
	for (i=0; i<AnzahlMarker; i++) if (marker[i].Marked) N++;
	printf("TestProgression: N=%d\n", N);
	Marker *ML[N];
	int vssg[N], vssk[N];
    TestProg *TProg[10], *P = new TestProg, *PBuff;
    for (i=0; i<AnzahlMarker; i++) if (marker[i].Marked) 
    {
		ML[j++] = marker + i;
		SMarker += marker[i].Line[0];
    }
    SMarker /= N;
    for (i=0; i<10; i++) 
    {
		TProg[i] = new TestProg;
		TProg[i]->vss = new int[N];
		TProg[i]->Jss = new int[N];
		TProg[i]->vs = new int[N];
		TProg[i]->Js = new int[N];
		TProg[i]->Iso = new int[N];
		TProg[i]->DD = new double[N];
		TProg[i]->SBD = 10 * N;
		TProg[i]->lState = new ElState*[N];
		TProg[i]->uState = new ElState*[N];
    }
    P->vss = new int[N];
    P->Jss = new int[N];
    P->vs = new int[N];
    P->Js = new int[N];
    P->Iso = new int[N];
    P->DD = new double[N];
	P->lState = new ElState*[N];
	P->uState = new ElState*[N];
    double Diff[N-1], NE, SDiff, SDiffg, SDiffk, DDg[N], DDk[N];
	ElState *lState, *uState;
    for (i=1; i<N; i++) Diff[i-1] = ML[i]->Line[0] - ML[0]->Line[0];
    //printf("Vor Hauptschleife\n");
    for (M=0; M < numMolecules; M++)
	{
		(Mol = MW->getMolecule(M))->getTermData(0, NC, NI, veu, Jeu, ELl);
		//printf("Nach getTermData, NC=%d, NI=%d, veu=%d, Jeu=%d\n", NC, NI, veu, Jeu);
		if (veu > mveu) mveu = veu;
		veu--;
		Jeu--;
		for (I=0; I<NI; I++) for (J=0; J<Jeu; J++) for (v=0; v<veu; v++)
    	{
			//printf("J=%d, v=%d\n", J, v);
			P->vss[0] = j = g = v;
			if (veu <= (k = v + 4)) k = veu - 1;
			SLinesg = SLinesk = SLines = ELl[0][I][v][J];
			kJ = J - 2;
			gJ = J + 2;
			vssg[0] = vssk[0] = -1;
			for (i=1; i<N; i++) 
			{
		    	//printf("i=%d\n", i);
		    	NE = ELl[0][I][v][J] - Diff[i - 1];
		    	//if (I==0 && J==36 && v==0) printf("Diff[%d]=%f\n", i, Diff[i]);
		    	while (ELl[0][I][j][J] > NE && j > 0 ) j--;
		    	if (ELl[0][I][j+1][J] - NE < NE - ELl[0][I][j][J]) 
		    	{
					P->vss[i] = j + 1;
					SLines += ELl[0][I][j+1][J];
					//if (I==0 && J==36 && v==12) printf("j+1=%d\n", j+1);
	    		}
	    		else 
	    		{
					P->vss[i] = j;
					SLines += ELl[0][I][j][J];
	 			}
				if (kJ >= 0)
				{
					while (ELl[0][I][k][kJ] > NE && k > 0) k--;
					if (ELl[0][I][k+1][kJ] - NE < NE - ELl[0][I][k][kJ]) vssk[i] = k + 1;
					else vssk[i] = k;
					if (fabs(NE - ELl[0][I][vssk[i]][kJ]) > fabs(NE - ELl[0][I][P->vss[i]][J]))
					{ 
						vssk[i] = -1;
						SLinesk += ELl[0][I][P->vss[i]][J];
					}
					else SLinesk += ELl[0][I][vssk[i]][kJ];
				}
				if (gJ <= Jeu)
				{
					while (ELl[0][I][g][gJ] > NE && g > 0) g--;
					if (ELl[0][I][g+1][gJ] - NE < NE - ELl[0][I][g][gJ]) vssg[i] = g + 1;
					else vssg[i] = g;
					if (fabs(NE - ELl[0][I][vssg[i]][gJ]) > fabs(NE - ELl[0][I][P->vss[i]][J]))
					{
						vssg[i] = -1;
						SLinesg += ELl[0][I][P->vss[i]][J];
					}
					else SLinesg += ELl[0][I][vssg[i]][gJ];
				} 
			}
			SLines /= N;
			SLinesg /= N;
			SLinesk /= N;
			SDiff = SMarker + SLines;
			SDiffg = SMarker + SLinesg;
			SDiffk = SMarker + SLinesk;
			P->SDiff = SDiff;
			P->SBD = SBDg = SBDk = 0.0;
			//printf("Mitte der Hauptschleife\n");
			for (i=0; i < N; i++) 
			{
		    	//printf("P->vss[%d]=%d\n", i, P->vss[i]);
				P->DD[i] = SDiff - ELl[0][I][P->vss[i]][J] - ML[i]->Line[0];
				//printf("vssg[%d]=%d\n", i, vssg[i]);
				if (vssg[i] > -1 && gJ <= Jeu) 
					DDg[i] = SDiffg - ELl[0][I][vssg[i]][gJ] - ML[i]->Line[0];
				else DDg[i] = SDiffg - ELl[0][I][P->vss[i]][J] - ML[i]->Line[0];
				//printf("vssk[%d]=%d\n", i, vssk[i]);
				if (vssk[i] > -1 && kJ >= 0) DDk[i] = SDiffk - ELl[0][I][vssk[i]][kJ] - ML[i]->Line[0];
				else DDk[i] = SDiffk - ELl[0][I][P->vss[i]][J] - ML[i]->Line[0];
	    		//if (I==0 && J==36 && v==12) printf("P->DD[%d]=%f\n", i, P->DD[i]);
	    		P->SBD += fabs(P->DD[i]);
				SBDg += fabs(DDg[i]);
				SBDk += fabs(DDk[i]);
			}
			//printf("Nach MWB\n");
			if (SBDg < P->SBD && SBDg < TProg[9]->SBD && gJ <= Jeu)
			{
				for (i=0, lState = Mol->getStateP(0); i<N; i++)
				{
					TProg[9]->Iso[i] = I;
					TProg[9]->Js[i] = J + 1;
					TProg[9]->DD[i] = DDg[i];
					TProg[9]->lState[i] = lState;
					if (vssg[i] != -1) 
					{
						TProg[9]->vss[i] = vssg[i];
						TProg[9]->Jss[i] = gJ;
					}
					else
					{
						TProg[9]->vss[i] = P->vss[i];
						TProg[9]->Jss[i] = J;
					}
				}
				TProg[9]->SBD = SBDg;
				TProg[9]->SDiff = SDiffg;
				for (i=9; i>0; i--)
	    		{
					if (TProg[i]->SBD < TProg[i-1]->SBD)
					{
		    			PBuff = TProg[i];
		    			TProg[i] = TProg[i-1];
		    			TProg[i-1] = PBuff;
					}
					else break;
	    		}
			}
			//printf("Nach g\n");
			if (SBDk < P->SBD && SBDk < TProg[9]->SBD && kJ >= 0)
			{
				for (i=0, lState = Mol->getStateP(0); i<N; i++)
				{
					TProg[9]->Iso[i] = I;
					TProg[9]->Js[i] = J - 1;
					TProg[9]->DD[i] = DDk[i];
					TProg[9]->lState[i] = lState;
					if (vssk[i] != -1) 
					{
						TProg[9]->vss[i] = vssk[i];
						TProg[9]->Jss[i] = kJ;
					}
					else
					{
						TProg[9]->vss[i] = P->vss[i];
						TProg[9]->Jss[i] = J;
					}
				}
				TProg[9]->SBD = SBDk;
				TProg[9]->SDiff = SDiffk;
				for (i=9; i>0; i--)
	    		{
					if (TProg[i]->SBD < TProg[i-1]->SBD)
					{
		    			PBuff = TProg[i];
		    			TProg[i] = TProg[i-1];
		    			TProg[i-1] = PBuff;
					}
					else break;
	    		}
			}
			//printf("Nach k\n");	
			if (P->SBD < TProg[9]->SBD)
			{
		    	for (i=0, lState = Mol->getStateP(0); i<N; i++)
	    		{
					P->Jss[i] = J;
					P->Iso[i] = I;
					P->Js[i] = -1;
					P->lState[i] = lState;
	    		}
	    		PBuff = P;
				P = TProg[9];
	    		TProg[9] = PBuff;
	    		for (i=9; i>0; i--)
	    		{
					if (TProg[i]->SBD < TProg[i-1]->SBD)
					{
		    			PBuff = TProg[i];
		    			TProg[i] = TProg[i-1];
		    			TProg[i-1] = PBuff;
					}
					else break;
	    		}
			}
		}
    }
    //printf("Nach Hauptschleife\n");
    Marker **MBuff, *FoundLines[mveu];
	j = 2 * (veu + 1);
    int FoundCount, vm, vn, jd=2, n, aJss = -1, lJ[j], nLev, li[j];
    double ZS, lE[j], FCF[j], IS, FS = 0.0, QF, bQF, BI, MI;
	bool more = true, Rem[j];
	Potential *lPot, *uPot;
    for (i=0; i<10; i++) 
    {
		if (TProg[i]->SBD < 10 * N)
		{
			FoundCount = 0;
			(Mol = TProg[i]->lState[0]->getMolecule())->getTermData(0, NC, NI, veu, Jeu, ELl);
			lPot = Mol->getPot(0);
			veu--;
			Jeu--;
			while (more)
			{
				j = AnzahlMarker - 2;
		   		if (TProg[i]->Js[0] == -1) more = false;
				if (aJss == TProg[i]->Jss[0])
				{
					more = false;
					for (n=1; TProg[i]->Jss[n] == aJss; n++) ;
					aJss = TProg[i]->Jss[n];
				}
				else aJss = TProg[i]->Jss[0];
				for (v=0; v<veu; v++)
		    	{
					if(ELl[0][TProg[i]->Iso[0]][v][aJss] == 0.0) break;
					ZS = TProg[i]->SDiff - ELl[0][TProg[i]->Iso[0]][v][aJss];
					while (ZS < marker[j].Line[0] && j>0) j--;
					if (marker[j+1].Line[0] - ZS < ZS - marker[j].Line[0]) j++;
					if ((fabs(marker[j].Line[0] - ZS) < ST || 
		    	    	 (LaserLinie != NULL ? fabs(marker[j].Line[0] - LaserLinie->Line[0]) < 0.5 
						: false)) && !marker[j].Marked && marker[j].HFLM > MinPeakHeight)
					{
		    			//printf("Gefunden: v=%d, J=%d, WN=%f\n", i, Jiu, marker[j].Line[0]);
		    			FoundLines[FoundCount] = marker + j;
		    			FoundLines[FoundCount]->DD = ZS - FoundLines[FoundCount]->Line[0];
						//FoundLines[FoundCount]->SNR = FoundLines[FoundCount]->HFLM / Rauschen;
		    			FoundLines[FoundCount]->Iso = TProg[i]->Iso[0];
		    			FoundLines[FoundCount]->vss = v;
						FoundLines[FoundCount]->Js = TProg[i]->Js[0];
						FoundLines[FoundCount]->LState = TProg[i]->lState[0];
						FoundLines[FoundCount]->UState = TProg[i]->uState[0];
		    			FoundLines[FoundCount++]->Jss = aJss;
					}
					if (j==AnzahlMarker-1) j--;
    			}
			}
		    MBuff = new Marker*[2*(N+FoundCount)];
		    l = FoundCount - 1;
	    	k = 0;
	    	for (j=0; j<N+FoundCount; j++)
	    	{
				if (k < N) vm = TProg[i]->vss[k];
				else vm = -1;
				if (l >= 0) vn = FoundLines[l]->vss;
				else vn = -1;
				if (vm > vn)
				{
		    		ML[k]->vss = TProg[i]->vss[k];
		    		ML[k]->Jss = TProg[i]->Jss[k];
		    		ML[k]->Iso = TProg[i]->Iso[k];
					ML[k]->Js = TProg[i]->Js[k];
		    		ML[k]->DD = TProg[i]->DD[k];
					//ML[k]->SNR = ML[k]->HFLM / Rauschen;
					ML[k]->LState = TProg[i]->lState[k];
					ML[k]->UState = TProg[i]->uState[k];
		    		MBuff[j] = ML[k++];
				}
				else MBuff[j] = FoundLines[l--];
	    	}
			vm = FoundCount + N;
			if (TProg[i]->Js[0] == -1)
			{
	    		if (MBuff[0]->Jss + 2 > Jeu) jd = -2;
	    		else jd = 2;
	    		for (j=0, k=vm, l=0; jd == 2 || j < vm; j++)
	    		{
					if (j > vm)
					{
			    		j = l = 0;
		    			k = vm;
		    			jd = -2;
		    			if (MBuff[0]->Jss - 2 < 0) break;
					}
					if (j == vm) break;
					ZS = MBuff[j]->Line[0] + ELl[0][MBuff[j]->Iso][MBuff[j]->vss][MBuff[j]->Jss]
			     		- ELl[0][MBuff[j]->Iso][MBuff[j]->vss][MBuff[j]->Jss + jd];
					//printf("MBuff[%d]->Line[0]=%f, ZS=%f\n", j, MBuff[j]->Line[0], ZS);
					while (marker[l].Line[0] < ZS && l < AnzahlMarker - 1) l++;
					if (l > 0) if (ZS - marker[l-1].Line[0] < marker[l].Line[0] - ZS) l--;
					//printf("Schleifenmitte, l=%d\n", l);
					if ((marker[l].HFLM > 0.5 * MBuff[j]->HFLM || MBuff[j] == LaserLinie) && 
							(fabs(ZS - marker[l].Line[0]) < ST || marker + l == LaserLinie))
					{
			    		MBuff[k++] = marker + l;
		    			marker[l].vss = MBuff[j]->vss;
		    			marker[l].Jss = MBuff[j]->Jss + jd;
		    			marker[l].Iso = MBuff[j]->Iso;
		    			marker[l].DD = ZS - marker[l].Line[0];
						//marker[l].SNR = marker[l].HFLM / Rauschen;
		    			marker[l].vs = -1;
		    			marker[l].Js = MBuff[j]->Jss + jd / 2;
						marker[l].LState = MBuff[j]->LState;
						marker[l].UState = MBuff[j]->UState;
		    			//printf("Gefunden: Linie bei %f\n", marker[l].Line[0]);
					}
					else for (n=0; n<N; n++) if (ML[n] == MBuff[j]) j = vm + 1;
					//printf("Schleifenende\n");
	    		}
	    		if (j==vm) 
	    		{
					for (j=0; j<vm; j++) MBuff[j]->Js = MBuff[vm]->Js;
					vm = k;
	    		}
	    		else for (j=0; j<vm; j++) MBuff[j]->Js = MBuff[j]->Jss;
			}
			//printf("Vor MBuff\n");
			J = MBuff[0]->Js;
			I = MBuff[0]->Iso;
			NE = ELl[0][I][MBuff[0]->vss][MBuff[0]->Jss] + MBuff[0]->Line[0];
			NS = Mol->getNumStates();
			if (J == MBuff[0]->Jss) 
			{
				for (v=0; v < veu; v++)  
				{
					lE[v] = ELl[0][I][v][J];
					if (lE[v] == 0.0) break;
					lJ[v] = J;
					li[v] = -1;
				}
				nLev = v;
				for (v=0; v < vm; v++) li[MBuff[v]->vss] = v;
			}
			else for (nLev = v = 0; v < veu; v++)
			{
				for (j=0, li[nLev] = li[nLev + 1] = -1; j < vm; j++) if (MBuff[j]->vss == v)
				{
					if (MBuff[j]->Jss < J) li[nLev] = j;
					else li[nLev + 1] = j;
				}
				if (J > 0 ? (lE[nLev] = ELl[0][I][v][lJ[nLev] = J-1]) > 0.0 : false) nLev++;
				if (J < Jeu ? (lE[nLev] = ELl[0][I][v][lJ[nLev] = J+1]) > 0.0 : false) nLev++;
			}
			for (s=1, bQF = 10.0, lvs = vs = -1, locs = 1e12; s < NS; s++) 
				if (MBuff[0]->Jss != J || Mol->getStateP(s)->getLambda() > 0) 
			{
				Mol->getTermData(s, NC, NI, veo, Jeo, ELu);
				if (I >= NI) continue;
				uPot = Mol->getPot(s);
				veo--;
				Jeo--;
				if (NC == 0) continue;
				if (MBuff[0]->Jss == J)
				{
					c = (NC > 1 ? 1 : 0);
					for (vs=0; (vs <= veo ? ELu[c][I][vs][J] < NE : false); vs++) ;
					if ((vs > 0 ? (vs <= veo ? ELu[c][I][vs][J] - NE > NE - ELu[c][I][vs-1][J] 
						: true) : false)) vs--;
					if (fabs(ELu[c][I][vs][J] - NE) > 50.0) vs = -1;
				}
				else
				{
					
					for (vs = c = 0; (vs <= veo ? ELu[0][I][vs][J] < NE : false); vs++) ;
					if ((vs > 0 ? (vs <= veo ? ELu[0][I][vs][J] - NE > NE - ELu[0][I][vs-1][J]
						: true) : false)) vs--;
					if (fabs(ELu[0][I][vs][J] - NE) > 50.0) vs = -1;
				}
				if (vs != -1)
				{
					//printf("Vor FCFTest\n");
					if (lPot != 0 && uPot != 0 && useIntensities)
					{
                        lPot->getFastFCF(I, uPot, J, c, ELu[c][I][vs][J], nLev, lE, lJ, FCF, NumWFPoints);
						//printf("Nach getFastFCF\n");
						for (j=0, IS = FS = 0.0; j < nLev; j++) 
							if (li[j] != -1 ? MBuff[li[j]] != LaserLinie 
											  && !MBuff[li[j]]->overlap : false)
						{
							IS += MBuff[li[j]]->HFLM;
							FS += FCF[j];
						}
						FS /= IS;
						for (j=0, QF = 0.0; j < nLev; j++) 
						{
							if (li[j] == -1) QF += FCF[j];
							else QF += fabs(FCF[j] - FS * MBuff[li[j]]->HFLM);
						}
					}
					else QF = 0.5;
					if (/*lvs == -1 ||*/ QF < bQF)
					{
						//printf("Fall QF<bQF\n");
						lvs = vs;
						BI = 2.0 * MinPeakHeight;
						MI = 0.5 * FS * MinPeakHeight;
						for (n=0, ocs = 0.0; n < nLev; n++) if ((j = li[n]) != -1)
						{
							if (MBuff[j]->HFLM < BI && FCF[n] < MI) Rem[j] = false;
							else
							{
								MBuff[j]->vs = vs;
								ocs += (MBuff[j]->oc = ELl[0][I][MBuff[j]->vss][MBuff[j]->Jss] 
										+ MBuff[j]->Line[0] 
										- (MBuff[j]->Jss == J ? ELu[c][I][vs][J] : ELu[0][I][vs][J]));
								Rem[j] = true;
							}
						}
						locs = ocs;
						uState = Mol->getStateP(s);
						bQF = QF;
						/*for (j = 0; j < nLev; j++) if (li[j] != -1 ? !Rem[li[j]] : false)
								for (n=j, lc = 0, BI = MI = 0.0; n < nLev; n++)
						{
							if (li[n] == -1)
							{
								if (FCF[n] > BI) 
								{
									BI = FCF[j];
									if (BI > MI) break;
									if (lc == 2 ? FCF[l2] < BI : false) lc = 1;
								}
							}
							else if (FCF[n] > BI)
							{
								BI = 0.0;
								if (FCF[n] > MI) MI = FCF[n];
								switch (++lc)
								{
									case 1:
										l1 = n;
										break;
									case 2:
										l2 = n;
										break;
									case 3:
										Rem[li[l1]] = Rem[li[l2]] = true;
									default:
										Rem[li[n]] = true;
										break;
								}
							}
						}*/
					}
					else if (QF == bQF)
					{
						//printf("Fall QF==bQF\n");
						for (j=0, ocs = 0.0; j<vm; j++)
							ocs += ELl[0][I][MBuff[j]->vss][MBuff[j]->Jss] + MBuff[j]->Line[0] 
								 - (MBuff[j]->Jss == J ? ELu[c][I][vs][J] : ELu[0][I][vs][J]);
						if (ocs < locs)
						{
							lvs = vs;
							for (j=0, ocs = 0.0; j<vm; j++)
							{
								MBuff[j]->vs = vs;
								ocs += (MBuff[j]->oc = ELl[0][I][MBuff[j]->vss][MBuff[j]->Jss] 
									 + MBuff[j]->Line[0] - (MBuff[j]->Jss == J ? ELu[c][I][vs][J] 
									                      					: ELu[0][I][vs][J]));
								Rem[j] = true;
							}
							locs = ocs;
							uState = Mol->getStateP(s);
						}
					}
					//printf("Nach FCFTest\n");
				}
			}
			if (lvs == -1) for (j=0; j<vm; j++) 
			{
				MBuff[j]->vs = -1;
				MBuff[j]->oc = 0.0;
				MBuff[j]->UState = 0;
			}
			else 
			{
				for (j=0; (j < vm ? Rem[j] : false); j++) ;
				for (n=j+1; n < vm; n++) if (Rem[n]) MBuff[j++] = MBuff[n];
				vm = j;
				for (j=0; j<vm; j++) MBuff[j]->UState = uState;
			}
			//printf("Nach MBuff\n");
	    	Prog.Insert(MBuff, vm, 10 - i);
		}
		delete[] TProg[i]->vss;
		delete[] TProg[i]->Jss;
		delete[] TProg[i]->vs;
		delete[] TProg[i]->Js;
		delete[] TProg[i]->Iso;
		delete[] TProg[i]->DD;
		delete[] TProg[i]->lState;
		delete[] TProg[i]->uState;
		delete TProg[i];
    }
    delete[] P->vss;
    delete[] P->Jss;
    delete[] P->vs;
    delete[] P->Js;
    delete[] P->Iso;
    delete[] P->DD;
	delete[] P->lState;
	delete[] P->uState;
    delete P;
	//printf("Vor ShowFound\n");
	ShowFound();
}

IntProg *Spektrum::SearchLongestProg(Marker *M)
{
	//printf("SearchLongestProgression, Spektrum=%s\n", getName().toAscii().data());
	IntProg PB[10], *PR = new IntProg, PT, PD;
	int I, v, J, N, NAL, nJ, n, m, Mo, NM = MW->getNumMolecules(), NC, NI, veu, Jeu, NP=0, uc;
	bool cbR, cbQ, cbP, uL[1000], FCFA;
	double T, aT, uT[1000], d, oT, ****Data;
	Marker **AL = NULL, *ML = 0;
	Molecule *Mol;
	GetMA(AL, NAL);
	//printf("NI=%d, veu=%d, Jeu=%d\n", NI, veu, Jeu);
	//printf("M->Line[0]=%f\n", M->Line[0]);
	/*bool P;
	if (fabs(M->Line[0] - 16856.127137) < 0.000001) P = true;
	else P = false;*/
	for (Mo=0; Mo < NM; Mo++)
	{
		Mol = MW->getMolecule(Mo);
		Mol->getTermData(0, NC, NI, veu, Jeu, Data);
		FCFA = Mol->FCFavailable();
		veu--;
		Jeu--;
		if (Data == 0) continue;
		for (uc=0; uc < NC; uc++) for (I=0; I<NI; I++) for (v=0; v<=veu; v++) for (J=0, aT = 0.0; J<=Jeu; J++)
		{
			//printf("I=%d, NI=%d, v=%d, veu=%d, J=%d, Jeu=%d\n", I, NI, v, veu, J, Jeu);
			if (J > 0 && v > 0 ? Data[uc][I][v][J] <= Data[uc][I][v][J-1] || Data[uc][I][v][J] <= Data[uc][I][v-1][J] 
				: false) break;
			PT.clear();
			PT.Mol = Mol;
			PT.lState = Mol->getStateP(0);
			PT.FC = (NC > 1 ? uc : -1);
			T = Data[uc][I][v][J] + M->Line[0];
			cbR = (J >= 2 ? isLI(I, n, nJ = J - 2, T, d, Data[uc], veu) : false);
			cbP = (J <= Jeu - 2 ? isLI(I, n, nJ = J + 2, T, d, Data[uc], veu) : false);
			if ((cbQ = isLI(I, n, J, T, d, Data[uc], veu)) 
						  || (PT.satellite = cbSat(I, v, J, *M, AL, NAL, ML)))
			{
				//if (cbP || cbQ || cbR) PT.satellite = false;
				cbQ = true;
				if (J >= 2) cbR = true;
				if (J <= Jeu - 2) cbP = true;
				//if (PT.satellite) printf("I=%d, v=%d, J=%d: satellite\n", I, v, J);
			}
			else PT.satellite = false;
			if (PT.satellite && ML != 0 && ML->UState != 0 && ML->UState->getLambda() == 0) cbQ = false;
			if (cbR || cbQ || cbP || Type == ThermalEmission)
			{
				//if (I==0 && J==197 && v==11) printf("I=%d, v=%d, J=%d\n", I, v, J);
				PT.Iso = I;
				N = SML(PT, I, J, T, aT = 2 * ST, veu, Data[uc]);
				//if (I==0 && J==197 && v==11) printf("N=%d\n", N);
				if ((N==0 || (N < 3 && !cbR && !cbP)) && !FCFA)  continue;
				//printf("T=%f, ELU[%d][%d][%d]=%f, M->Line[0]=%f\n", T, I, v, J, 
				//ELU[I][v][J], M->Line[0]);
				for (n=0; n<N; n++) 
				{
					uL[n] = true;
					uT[n] = PT.marker[n]->Line[0] + Data[uc][I][PT.vss[n]][J];
				}
				//printf("T=%f\n", T);
				oT = T;
				while (true)
				{
					aT = T;
					for (m=n=0, T=0.0; n<N; n++) if (uL[n]) 
					{
						m++;
						T += uT[n];
						//printf("uT[%d]=%f ", n, uT[n]);
					}
					T /= m;
					//printf("T=%f, m=%d\n", T, m);
					if (fabs(T - oT) > ST)
					{
						m = 0;
						break;
					}
					if (m==0) break;
					if (fabs(aT - T) < 0.0000001) break;
					for (n=0; n<N; n++) uL[n] = (fabs(uT[n] - T) < ST ? true : false);
				}
				if (m == 0) continue; 
				for (n=0; (n<N ? uL[n] : false); n++) ;
				for (m=n+1; m<N; m++) if (uL[m])
				{
					PT.marker[n] = PT.marker[m];
					PT.vss[n] = PT.vss[m];
					PT.Jss[n] = PT.Jss[m];
					PT.overlap[n++] = PT.overlap[m];
				}
				PT.N = N = n;
				//if (v==0 && J == 68) printf("J'=68: N=%d\n", N);
				//if (v==25 && I==5 && J==91) for (n=0; n<N; n++)
					//printf("marker[%d]-marker=%d\n", n, PT.marker[n] - marker);
				//printf("I=%d, NI=%d, v=%d, veu=%d, J=%d, Jeu=%d\n", I, NI, v, veu, J, Jeu);
				PT.Np = N;
				if (2 * N > PB[9].N || (PT.satellite && !PB[9].satellite && N>1)) 
				{
					if (cbR) 
					{
						//if (I==0 && J==197 && v==11) printf("cbR, v=%d, N=%d\n", v, N);
						SML(PT, I, nJ = J - 2, T, ST, veu, Data[uc]);
						if (PT.N > PT.Np)
						{
							PT.Js = J-1;
							//if (PT.Js == 23) for (n=0; n < PT.N; n++) 
								//	printf("Js=23: vss=%d, Jss=%d, f=%f\n", 
									//	   PT.vss[n], PT.Jss[n], PT.marker[n]->Line[0]);
							PD.copy(PT);
							eNGL(PD, AL, NAL, Data[uc]);
							//if (I==0 && J==197 && v==11) printf("PT.G=%d\n", PT.G);
							if ((PD.G >= PB[9].G /*&& PT.G >= 3*/) || (PD.satellite && !PB[9].satellite))
							{
								chkDoubletts(PD);
								//if (I==0 && J==197 && v==11) printf("NGD=%d\n", PD.NGD);
								if (PD.isIncluded(M))
								{	
									cFQS(PD);
									for (n = NP - 1; (n>=0 ? PD.isbetter(PB + n) : false); n--) ;
									for (m = (NP < 10 ? NP : 9), n++; m > n; m--) PB[m].copy(PB[m-1]);
									if (n < 10) PB[n].copy(PD);
									if (NP < 10) NP++;
								}
							}
						}
						PT.N = PT.Np;
					}
					if (cbP)
					{
						//if (I==5 && J==91 && v==25) printf("cbP, v=%d, NAL=%d\n", v, NAL);
						//if (v==25 && I==5 && J==91) for (n=0; n<PT.Np; n++)
							//printf("marker[%d]-marker=%d\n", n, PT.marker[n] - marker);
						SML(PT, I, nJ = J + 2, T, ST, veu, Data[uc]);
						//if (I==5 && v==25 && J == 91) printf("J''=91: N=%d\n", PT.N);
						if (PT.N > PT.Np)
						{
							PT.Js = J+1;
					//if (v==23 && P && J==93) printf("marker[5]-marker=%d\n", PT.marker[5] - marker);
							PD.copy(PT);
							eNGL(PD, AL, NAL, Data[uc]);
							//if (v==25 && J==91 && I==5) 
								//printf("J=91, nach eNGL: PT.N=%d, PT.G=%d\n", PT.N, PT.G); 
							if (PD.G >= PB[9].G || (PD.satellite && !PB[9].satellite))
							{
								chkDoubletts(PD);
								//if (v==25 && J==91 && I==5) 
							//printf("J=91, nach chkDoubletts: PD.N=%d, PD.NGD=%d\n", PD.N, PD.NGD); 
								if (PD.isIncluded(M))
								{
									cFQS(PD);
									for (n = NP - 1; (n>=0 ? PD.isbetter(PB + n) : false); n--) ;
									for (m = (NP < 10 ? NP : 9), n++; m > n; m--) PB[m].copy(PB[m-1]);
									if (n < 10) PB[n].copy(PD);
									if (NP < 10) NP++;
								}
							}
						}
						PT.N = PT.Np;
					}
					if (cbQ && ((N >= PB[9].N) || (PT.satellite && !PB[9].satellite))
									   && !OnlyDoublet)
					{
						//if (P && v==23 && J==93) printf("cbQ\n");
						PT.NGD = 0;
						PT.Js = J;
						eNGL(PT, AL, NAL, Data[uc]);
						//printf("Nach eNGL\n");
						if (PT.G >= PB[9].G)
						{
							cFQS(PT);
							for (n = NP - 1; (n>=0 ? PD.isbetter(PB + n) : false); n--) ;
							for (m = (NP < 10 ? NP : 9), n++; m > n; m--) PB[m].copy(PB[m-1]);
							if (n < 10) PB[n].copy(PT);
							if (NP < 10) NP++;
							//printf("cbQ, I=%d, v=%d, J=%d, N=%d\n", I, v, J, N);
						}
					}
					//printf("Nach cbQ\n");
				}
			}
		}
	}
	//printf("SLP: vor Ende\n");
	//printf("PB->N=%d, PB->Js=%d\n", PB->N, PB->Js);
	delete[] AL;
	if (PB[0].N == 0)
	{
		delete PR;
		return NULL;
	}
	//if (PB->FQS == 0.0) cFQS(*PB);
	double ED;
	int oNC, oNI, veo, Jeo, c, bs, bc;
	double ****oData;
	for (n=0; n < NP; n++)
	{
		/*if ((m = Assignvs(PB[n])) < 2)
		{
			if (m==0) 
			{
				PR->copy(PB[n]);
				break;
			}
			if (PB[n].QF > PR->QF) PR->copy(PB[n]); 
		}*/
		N = PB[n].Mol->getNumStates();
		PB[n].Mol->getTermData(0, NC, NI, veu, Jeu, Data);
		uc = (PB[n].FC != -1 ? PB[n].FC : 0);
		for (m=c=0, PB[n].uT = 0.0; m < PB[n].N; m++) if (!PB[n].overlap[m])
		{
			PB[n].uT += PB[n].marker[m]->Line[0] 
					  + Data[uc][PB[n].Iso][PB[n].vss[m]][PB[n].Jss[m]];
			c++;
		}
		PB[n].uT /= c;
		PB[n].QF = 10.0;
		PB[n].uState = 0;
		PB[n].vs = bs = -1;
		for (m=1; m<N; m++)
		{
			PB[n].Mol->getTermData(m, oNC, oNI, veo, Jeo, oData);
			if (PB[n].Js >= Jeo || PB[n].Iso >= oNI) continue;
			c = (PB[n].FC != -1 ? PB[n].FC : (PB[n].Js == PB[n].Jss[0] && oNC > 1 ? 1 : 0));
			for (v=0; (v < veo ? oData[c][PB[n].Iso][v][PB[n].Js] < PB[n].uT : false);
				 v++) ;
			if (v < veo ? (ED = fabs(oData[c][PB[n].Iso][v][PB[n].Js] - PB[n].uT)) 
				< PB[n].QF : false)
			{
				PB[n].QF = ED;
				PB[n].vs = v;
				PB[n].uState = PB[n].Mol->getStateP(m);
				bs = m;
				bc = c;
			}
			if (v>0 ? ((ED = fabs(oData[c][PB[n].Iso][v-1][PB[n].Js] - PB[n].uT)) 
				< PB[n].QF) : false)
			{
				PB[n].QF = ED;
				PB[n].vs = v-1;
				PB[n].uState = PB[n].Mol->getStateP(m);
				bs = m;
				bc = c;
			}
		}
		if (bs > 0)
		{
			PB[n].Mol->getTermData(bs, oNC, oNI, veo, Jeo, oData);
			for (m=0; m < PB[n].N; m++) 
				PB[n].oc[m] = Data[uc][PB[n].Iso][PB[n].vss[m]][PB[n].Jss[m]] 
							+ PB[n].marker[m]->Line[0]
							- oData[bc][PB[n].Iso][PB[n].vs][PB[n].Js];	
		}
		PB[n].QF += 1.0;
		PB[n].QF /= (PB[n].NGD > 0 ? PB[n].G * PB[n].NGD : PB[n].G);
		for (m=0; m < PR->N; m++) PR->SNR[m] = PR->marker[m]->HFLM / Rauschen;
	}
	for (n=1, m=0, ED = PB[0].QF; n < NP; n++) if (PB[n].QF < ED)
	{
		m=n;
		ED = PB[n].QF;
	}
	PR->copy(PB[m]);
	for (n=0; n < NP; n++) Prog.Insert(PB + n, (n<m ? (n==m ? 0 : n+1) : n));
	//if (PB->satellite) printf("PB: satellite\n");
	//printf("Ende SLP\n");
	return PR;
}

bool Spektrum::chkPrep()
{
	printf("chkPrep\n");
	if (MW != 0 ? MW->getNumMolecules() == 0 : true)
	{
		
		return false;
	}
	if (Daten->GetDSL() == 0)
	{
		QMessageBox::information(this, "MolSpektanalysis", 
				 "Es muss zuerst ein Spektrum geladen werden!", QMessageBox::Ok);
		return false;
	}
	if (AnzahlMarker == 0) editFind();
	else if (!DMarkers) DisplayMarkers();
	printf("Ende chkPrep\n");
	return true;
}

void Spektrum::rSortP(IntProg &P)
{
	int n, b=1;
	Marker *B;
	bool bb;
	while (b!=-1)
	{
		b=-1;
		for (n=1; n<P.N; n++)
		{
			if (P.marker[n-1]->Line[0] > P.marker[n]->Line[0])
			{
				B = P.marker[n-1];
				P.marker[n-1] = P.marker[n];
				P.marker[n] = B;
				b = P.vss[n-1];
				P.vss[n-1] = P.vss[n];
				P.vss[n] = b;
				b = P.Jss[n-1];
				P.Jss[n-1] = P.Jss[n];
				P.Jss[n] = b;
				bb = P.overlap[n-1];
				P.overlap[n-1] = P.overlap[n];
				P.overlap[n] = bb;
			}
		}
	}
}

void Spektrum::getUnassignedLines(Marker **&M, int &N)
{
	//printf("getUnassignedLines\n");
	int n;
	//printf("GMarker[317]=%d\n", GMarker[317]);
	for (N=0, n=1; GMarker[n] != &MStop; n++) if (!GMarker[n]->DisplayData && GMarker[n] != LaserLinie)
				N++;
	//printf("N=%d\n", N);
	M = new Marker*[N];
	for (n=1, N=0; GMarker[n] != &MStop; n++) 
		if (!GMarker[n]->DisplayData && GMarker[n] != LaserLinie) M[N++] = GMarker[n];
	//printf("Ende getUnassignedLines\n");
}

void Spektrum::assignBand()
{
	int n, NL = 0;
	Marker **lines = new Marker*[1000];
	double E1, E2, ED1, ED2, EDD, SE;
	for (n=0; n < AnzahlMarker && NL < 1000; n++) if (marker[n].Marked) lines[NL++] = marker + n;
	if (NL < 2)
	{
		QMessageBox::information(this, "MolSpektAnalysis", 
				"Error: For spectra of this type at least two lines have to be selected to use this assingment function!");
		delete[] lines;
		return;
	}
	E1 = lines[NL - 1]->Line[0];
	ED1 = E1 - lines[NL - 2]->Line[0];
	EDD = (NL > 2 ? lines[NL - 2]->Line[0] - lines[NL - 3]->Line[0] - ED1 : 0.0);
	for (SE = E1 + ED1 + EDD; (NL < 1000 ? (lines[NL] = getMarker(SE)) != 0 : false); NL++)
	{
		E2 = E1;
		E1 = lines[NL]->Line[0];
		ED2 = ED1;
		ED1 = E2 - E1;
		EDD = ED2 - ED1;
		SE = E1 + ED1 + EDD;
	}
	E1 = lines[0]->Line[0];
	ED1 = E1 - lines[1]->Line[0];
	EDD = (NL > 2 ? lines[1]->Line[0] - lines[2]->Line[0] - ED1 : 0.0);
	for (SE = E1 + ED1 + EDD; (NL < 1000 ? (lines[NL] = getMarker(SE)) != 0 : false); NL++)
	{
		E2 = E1;
		E1 = lines[NL]->Line[0];
		ED2 = ED1;
		ED1 = E2 - E1;
		EDD = ED2 - ED1;
		SE = E1 + ED1 + EDD;
	}
	for (n=0; n < NL; n++) lines[n]->Marked = true;
	FokusProg(lines, NL);
	delete[] lines;
}

Marker* Spektrum::getMarker(double SearchEnergy, double Tol)
{
	int n;
	if (Tol <= 0.0) Tol = ST;
	for (n=0; (n < AnzahlMarker ? marker[n].Line[0] < SearchEnergy : false); n++) ;
	if (n>0 ? SearchEnergy - marker[n-1].Line[0] < marker[n].Line[0] - SearchEnergy : false) n--;
	return (fabs(SearchEnergy - marker[n].Line[0]) < Tol ? marker + n : 0);
}

void Spektrum::SingleSLP()
{
	//printf("SingleSLP\n");
    if (!chkPrep()) return;
	/*if (Type == ThermalEmission || Type == NormalizedAbsorption)
	{
		assignBand();
		return;
	}*/
	Marker **ML = NULL, *LL = NULL, *AL;
	int n, N;
	IntProg *P = NULL;
	getMarkedLines(ML, N);
	if (N == 0) getUnassignedLines(ML, N);
	for (AL = NULL; P == NULL; AL = NULL)
	{
		for (n=0; n<N; n++) 
			if ((LL != NULL ? ML[n]->HFLM <= LL->HFLM && ML[n] != LL : true)
						  && (AL != NULL ? ML[n]->HFLM > AL->HFLM : true)) AL = ML[n];
		if (AL != NULL) P = SearchLongestProg(LL = AL);
		else break;
	}
	if (P == NULL)
	{
		QMessageBox::information(this, "MolSpektanalysis", 
				 "Es konnte keine neue Progression gefunden werden!", QMessageBox::Ok);
		return;
	}
	MarkProgression(*P);
	FokusProg(P->marker, P->N);
	delete P;
}

void Spektrum::MultiSLP()
{
	printf("MultiSLP\n");
	if (!chkPrep()) return;
	Marker **ML = NULL;
	int N, n;
	getMarkedLines(ML, N);
	printf("Vor Hauptschleife\n");
	if (N>0) 
	{
		IntProg *PB;
		for (n=0; n<N; n++)
		{
			PB = SearchLongestProg(ML[n]);
			if (PB == 0) continue;
			MarkProgression(*PB);
			Prog.Insert(PB->marker, PB->N, PB->G);
			PB->marker = NULL;
			delete PB;
		}
	}	
	else
	{
		IntProg **P = new IntProg*[1000];
		MSLP(P, N);
		printf("Anzahl der gefundenen Progressionen: %d\n", N);
		for (n=0; n<N; n++)
		{
			ClearMarked();
			MarkProgression(*P[n]);
			Prog.Insert(P[n]->marker, P[n]->N, P[n]->G);
			P[n]->marker = NULL;
			delete P[n];
		}
		printf("Nach Prog.Insert\n");
		delete[] P;
	}
	ShowFound();
}

void Spektrum::AutoSLP()
{
	if (MW != 0 ? MW->getNumMolecules() == 0 : true)
	{
		QMessageBox::information(this, "MolSpektanalysis", 
				 "Es müssen zuerst Moleküldaten mit Termenergien oder Dunham-Koeffizienten für den Grundzustand geladen werden!",
	  				QMessageBox::Ok);
		return;
	}
	QString Dir = QFileDialog::getExistingDirectory(this, 
			"Wähle das Verzeichnis mit den zu Durchsuchenden Spektren.");
	if (Dir.isEmpty()) return;
	ASLP(Dir);
}

void Spektrum::ASLP(QString Dir)
{
	int i, n, N;
	printf("ASLP: Dir=%s\n", Dir.toLatin1().data());
	QString FP, F, Ext;
	IntProg **P = new IntProg*[1000];
	QDir Path(Dir);
	QFileInfoList SList = Path.entryInfoList();
	for (i=0; i<SList.size(); i++)
	{
		printf("Neu Datei\n");
		FP = Dir + "/" + (F = SList[i].fileName());
		if (SList[i].isDir()) {if (F.left(1) != ".") ASLP(FP); printf("Dir: %s\n", F.toLatin1().data());}
		else if (SList[i].isFile()) 
		{
			printf("Datei gefunden: %s\n", F.toLatin1().data());
			Ext = F.right(F.length() - F.indexOf('.') - 1);
			if (Ext == "dat" || Ext == ".spect")
			{
				readData(FP);
				printf("Nach Open\n");
				if (Daten->GetDSL() == 0) continue;
				if (AnzahlMarker == 0) editFind();
				else if (!DMarkers) DisplayMarkers();
				MSLP(P, N);
				for (n=0; n<N; n++)
				{
					ClearMarked();
					MarkProgression(*P[n]);
					addMarker(true);
					delete P[n];
				}
			}
		}
		else printf("Keine Datei: %s\n", F.toLatin1().data());
	}
	delete[] P;
}

void Spektrum::MSLP(IntProg **&P, int &NP)
{
	//printf("Beginn MSLP\n");
	Marker **ML = NULL, *MB, **AL = NULL;
	IntProg *B;
	bool S = true;
	int n, m, N, NAL;
	double ***ELU;
	//printf("Vor Marker\n");
	for (n=0; n < AnzahlMarker; n++) marker[n].b = true;
	//printf("Vor getUnnassignedLines\n");
	getUnassignedLines(ML, N);
	//printf("Nach getUnassignedLines\n");
	MB = ML[0];
	while(MB != NULL)
	{
		MB = NULL;
		for (n=1; n<N; n++) if (ML[n]->HFLM > ML[n-1]->HFLM)
		{
			MB = ML[n];
			ML[n] = ML[n-1];
			ML[n-1] = MB;
		}
	}
	NP = 0; 
	while (S)
	{
		//printf("Beginn\n");
		S = false;
		ClearMarked();
		for (n=0; n<NP; n++) MarkProgression(*P[n]);
		for (n=0; n<N; n++)
		{
			//printf("n=%d, N=%d\n", n, N);
			if (NP == 1000) break;
			if (!ML[n]->DisplayData && ML[n]->b) 
			{
				if ((P[NP] = SearchLongestProg(ML[n])) != NULL) 
				{
					//printf("Gefunden: vs=%d, Js=%d, N=%d, oc=%f\n", 
						 //  P[NP]->vs, P[NP]->Js, P[NP]->N, P[NP]->oc[0]);
					if (P[NP]->G >= 3 && P[NP]->NGD >= 2)
					{
						if (fabs(P[NP]->oc[0]) > 5.0) P[NP]->uState = 0;
						MarkProgression(*P[NP++]);
					}
				}
				//else printf("Nichts gefunden\n");
				S = true;
				ML[n]->b = false;
			}
		}
		//printf("Mitte\n");
		for (B = P[0]; B != NULL;) for (B = NULL, n=1; n < NP; n++) if (P[n]->isbetter(P[n-1]))
		{
			B = P[n-1];
			P[n-1] = P[n];
			P[n] = B;
		}
		//printf("Vor Schleife\n");
		for (n = NP - 1; n >= 0; n--)
		{
			ClearMarked();
			for (m=0; m < NP; m++) if (m != n && P[m] != NULL) MarkProgression(*P[m]);
			//printf("Nach MP, n=%d, NP=%d\n", n, NP);
			GetMA(AL, NAL);
			//printf("Nach GetMa\n");
			eNGL(*P[n], AL, NAL, ELU = P[n]->lState->getTermTable()->getData()[0]);
			//printf("Nach eNGL: vs=%d, Js=%d, N=%d\n", P[n]->vs, P[n]->Js, P[n]->N);
			if (P[n]->G < 3)
			{
				//printf("P[%d]->G=%d, NP=%d\n", n, P[n]->G, NP);
				delete P[n];
				P[n] = NULL;
			}
			else if (P[n]->N != P[n]->Np) 
			{
				chkDoubletts(*P[n]);
				//printf("Nach chkDoubletts: vs=%d, Js=%d, N=%d\n", P[n]->vs, P[n]->Js, P[n]->N);
				if (P[n]->NGD < 2)
				{
					//printf("P[%d]->G=%d, NP=%d\n", n, P[n]->G, NP);
					delete P[n];
					P[n] = NULL;
				}
			}
			//printf("Ende Schleife\n");
		}
		//printf("Nach Schleife\n");
		for (n=0; (n < NP ? P[n] != NULL : false); n++) ;
		for (m=n+1; m < NP; m++) if (P[m] != NULL) P[n++] = P[m];
		NP = n;
		//printf("Ende\n");
	}
	//printf("Ende MSLP\n");
	delete[] ML;
	delete[] AL;
}

void Spektrum::AverageWeighted(int NSpectra, QString *SpectFiles, int NLines, int **Lines, ElState *LS, 
						 ElState *US, int Iso, int Comp, int vs, int Js)
{
	if (LS == 0 || US == 0 || NLines <= 0 || Lines == 0 || Iso < 0 || Comp < 0 || vs < 0 || Js < 0)
		 return;
	TermTable *UT = US->getTermTable(), *LT = LS->getTermTable();
	if (UT == 0 || LT == 0 ? true 
		   : Comp >= UT->getNumComp() || Iso >= UT->getNumIso() || vs > UT->getMaxv() 
		   || Js > UT->getMaxJ())
		return;
	double *I=0, *R=0, F1, F2, ****UD = UT->getData(), ****LD = LT->getData(), **lines, Ml, *B;
	double Lf[NLines], V;
	int J, v, n, i, j, k, l, Mv = LT->getMaxv(), MJ = LT->getMaxJ(), DSL = 0, **il, N;
	int Lp[NLines];
	if (Lines[0][1] != Js)
	{
		lines = Create(Ml = 2*(Mv+1), 2);
		il = CreateInt(Ml, 2);
		for (v = Mv, n=0; v >= 0; v--) for (J = Js + 1; J >= Js - 1; J-=2) 
				if (J <= MJ ? LD[0][Iso][v][J] != 0.0 : false)
		{
			lines[n][0] = UD[Comp][Iso][vs][Js] - LD[0][Iso][v][J] - 1.0;
			lines[n][1] = lines[n][0] + 2.0;
			n++;
		}
	}
	else
	{
		lines = Create(Ml = Mv + 1, 2);
		il = CreateInt(Ml, 2);
		for (v = Mv, n=0; v >= 0; v--) if (LD[0][Iso][v][Js] != 0.0)
		{
			lines[n][0] = UD[Comp][Iso][vs][Js] - LD[0][Iso][v][Js] - 1.0;
			lines[n][1] = lines[n][0] + 2.0;
			n++;
		}
	}
	B = lines[0];
	while (B!=0) for (i=1, B=0; i<n; i++) if (lines[i-1][0] > lines[i][0])
	{
		B = lines[i-1];
		lines[i-1] = lines[i];
		lines[i] = B;
	}
	for (i=j=0; i < NLines; i++) 
		if (Lines[i][0] >= 0 && Lines[i][0] <= Mv && Lines[i][1] >= 0 && Lines[i][1] <= MJ)
	{
		Lf[j] = - LD[0][Iso][Lines[i][0]][Lines[i][1]];
		if (Lf[j] != 0.0) Lf[j++] += UD[Comp][Iso][vs][Js];
	}
	NLines = j;
	for (i=0; i < NSpectra; i++)
	{
		Daten->reinit();
		readData(SpectFiles[i]);
		if (Daten->GetDSL() < DSL) continue;
		if (R == 0)
		{
			DSL = Daten->GetDSL() - AnzahlMarker;
			R = new double[DSL];
			I = new double[DSL];
			for (j=k=l=0; j < Daten->GetDSL(); j++) if (!Daten->GetMarked(j))
			{
				I[k] = 0.0;
				R[k] = Daten->GetValue(j, 0);
				if (R[k] > lines[l][0]) il[l][0] = j;
				else 
				{
					if (R[k] > lines[l][1]) il[l][1] = j;
					l++;
				}
				k++;
			}
			for (j=k=0; j < NLines; j++)
			{
				while (R[k] < Lf[j] && k < DSL - 1) k++;
				while (R[k] > Lf[j]) k--;
				if (Lf[j] - R[k] > R[k+1] - Lf[j]) Lp[j] = k+1;
				else Lp[j] = k;
			}
		}
		for (j=l=N=0, F1 = F2 = 0.0; j < DSL; j++) 
		{
			if (j < il[l][0])
			{
				V = Daten->GetValue(j, 1);
				F1 += V*V;
				N++;
			}
			else while (j < il[l][0])
			{
				while (j < il[l][1]) j++;
				l++;
			}
		}
		for (j=0; j < NLines; j++) F2 += Daten->GetValue(Lp[j], 1);
		F1 = F2 * double(N*N) / (F1 * F1);
		for (j=k=0; k < DSL; j++) if (!Daten->GetMarked(j)) I[k++] += Daten->GetValue(j, 1);
	}
	Daten->reinit();
	for (i=0; i < DSL; i++) Daten->AddValue(R[i], I[i], false);
	delete[] R;
	delete[] I;
	Destroy(lines, Ml);
	Destroy(il, Ml);
}

void Spektrum::MarkProgression(IntProg &P)
{
	int n;
	IsoTab *Iso = (P.Mol != 0 ? P.Mol->getIso() : 0);
	QString IsoName = (Iso != 0 ? Iso->texName[P.Iso] : "");
	for (n=0; n<P.N; n++) 
	{
		P.marker[n]->DD = P.DD[n];
		P.marker[n]->Iso = P.Iso;
		P.marker[n]->Js = P.Js;
		P.marker[n]->Jss = P.Jss[n];
		P.marker[n]->oc = P.oc[n];
		P.marker[n]->overlap = P.overlap[n];
		P.marker[n]->satellite = P.satellite;
		//P.marker[n]->SNR = P.SNR[n];
		P.marker[n]->vs = P.vs;
		P.marker[n]->vss = P.vss[n];
		P.marker[n]->DisplayData = true;
		P.marker[n]->Marked = true;
		P.marker[n]->LState = P.lState;
		P.marker[n]->UState = P.uState;
		P.marker[n]->Mol = P.Mol;
		P.marker[n]->lState = (P.lState != 0 ? P.lState->getName() : "");
		P.marker[n]->uState = (P.uState != 0 ? P.uState->getName() : "");
		P.marker[n]->IsoName = IsoName;
		P.marker[n]->uncertainty = 0.0;
		P.marker[n]->FC = P.FC;
	}
	if (Iso != 0) delete Iso;
	//if (P.satellite) printf("MarkProgression: satellite\n");
}
	
void Spektrum::FindSat()
{
    double R = 10.0;
    //int OMEGA = 1;
    //printf("Beginn FindSatelites\n");
	if (AnzahlMarker == 0)
	{
		printf("FindSatelites kann nicht ausgeführt werden, da keine Marker vorhanden sind!\n");
		return;
	}
	int r = 0, MJ, j, fjd, k, i, di=1, n=0, N=0, Nv, NL, NH, FPInd = 0, Jeu, vs, st = 20, en = 20;
    Molecule *Mol;
	ElState *LState, *UState;
	for (i=0; i<AnzahlMarker; i++) if (marker[i].Marked && marker[i].DisplayData) 
	{
		N+=2;
		j=i;
	}
	double ****TD = 0, ***ELU;
	vs = marker[j].vs;
	LState = marker[j].LState;
	UState = marker[j].UState;
	if ((Mol = marker[j].Mol) != 0) Mol->getTermData(0, i, j, k, Jeu, TD);
	if (TD == 0 || i==0)
	{
		printf("Spektrum::FindSat() error: The necessary molecular data is missing!\n");
		return;
	}
	else ELU = TD[0];
    int AL = (int)((double)(4 * N * AnzahlMarker) * R / 
		   (marker[AnzahlMarker - 1].Line[0] - marker[0].Line[0]));
    Marker *FoundProg[41][N], *TEOMarker[AL], *MBuff;
	printf("N = %d\n", N);
	if (N < 12) 
	{
		QMessageBox::information( this, "MolSpektAnalysis", 
					"Diese Funktion benötigt mindestens vier zugeordnete und markierte Übergänge!", 
					QMessageBox::Ok);
		return;
	}
    double MaxDiff = 0.0, B, TEO[AL], MinSP, MaxSP, ProgPos, DBuff1, DBuff2;
    double Intens, AIntens, IntensB;
    double xa[41], ya[41], x, y, d=0.0, diff[41][3], SPunkt, XSum;
	int va[N], aJss, Jssa[AL], Iso, TEOPos[AL], TEOCount[N], TEOGC;
	int iBuff, BInd, Bk = 0, BDens, Dens, D[N];
    bool F = true;
    for (i=0; i<41; i++) 
    {
		ya[i] = 0.0;
		diff[i][0] = diff[i][1] = diff[i][2] = 0.0;
		for (n=0; n<N; n++) FoundProg[i][n] = NULL;
    }
    n = 0;
    for (i=0; i<AnzahlMarker; i++) if (marker[i].Marked && marker[i].DisplayData)
    {
		FoundProg[20][n++] = marker + i;
		ya[20] += marker[i].Line[0] + ELU[marker[i].Iso][marker[i].vss][marker[i].Jss];
		B = fabs(marker[i].DD);
		if (B > MaxDiff) MaxDiff = B;
	}
	NL = n;
	for (MBuff = marker + i; MBuff != 0; MBuff = 0) for (i=1; i<n; i++) 
			if (FoundProg[20][i-1]->vss < FoundProg[20][i]->vss)
	{
		MBuff = FoundProg[20][i-1];
		FoundProg[20][i-1] = FoundProg[20][i];
		FoundProg[20][i] = MBuff;
	}
	for (i=0; i<n; i++) va[i] = FoundProg[20][i]->vss;
    for (Nv=1, n=1; n<NL; n++) if (va[n] != va[n-1]) Nv++; 
	MaxDiff *= 2;
    if (MaxDiff > ST) MaxDiff = ST;
    ya[20] /= n;
    if (1 < (fjd = (FoundProg[20][0]->Jss < FoundProg[20][0]->Js ? 
		    FoundProg[20][0]->Js - FoundProg[20][0]->Jss
		: FoundProg[20][0]->Jss - FoundProg[20][0]->Js)))
    {
		r = 2;
		MJ = FoundProg[20][0]->Jss;
    }
    else MJ = FoundProg[20][0]->Js;
    N = NL;
    NH = N / 2;
    if (NH < 3) NH = 3;
    Iso = FoundProg[20][0]->Iso;
    j = MJ + 1;
    xa[20] = MJ * (MJ + 1);
    //printf("Vor Hauptschleife, N=%d, MaxDiff=%f\n", N, MaxDiff);
	while (r==0 || di!=-1 || (F  && FPInd != 0 && j != 0 )) 
    {
		if (j > Jeu)
		{
	    	di = -1;
	    	j = MJ - 1;
		}
		if (r==0 && di == -1 && (!F || FPInd == 0 || j==0))
		{
	    	if (fjd == 0) 
	    	{
				N = 2 * Nv;
				NH *= 2;
				for (i=N-1, n=NL-1; n >= 0; n--) for (k=0; k<2; k++) va[i--] = va[n];
	    	}
	    	else 
	    	{
				N = Nv;
				NH /= 2;
				for (n=1, k=0; n<N; n++)
				{
		    		k++;
					if (va[k] == va[k-1]) k++;
					va[n] = va[k];
				}
	    	}
	    	r = 1;
	    	di = 1;
	    	F = true;
	    	j = MJ;
	    	st = en = 20;
	    	for (i=0; i<41; i++) diff[i][0] = diff[i][1] = diff[i][2] = 0.0;
		}
		FPInd = 20 + j - MJ;
		x = j * (j + 1);
		printf("j=%d\n", j);
		TEOGC = 0;
		DBuff1 = DBuff2 = 0;
		if (st < en)
		{		    
	    	if (FPInd > 20) 
			{
	           ProgPos = ya[FPInd-1] + (diff[FPInd-1][0] + diff[FPInd-1][1])  * (x - xa[FPInd-1]);
			   printf("FPInd=%d, ProgPos=%f, ya[FPInd-1]=%f, diff[FPInd-1][0]=%f, diff[FPInd-1][1]=%f, x=%f, xa[FPInd-1]=%f\n", 
					  FPInd, ProgPos, ya[FPInd-1], diff[FPInd-1][0], diff[FPInd-1][1], x, xa[FPInd-1]);
			}
	    	else 
			{
				ProgPos = ya[FPInd+1] - (diff[FPInd+2][0] - diff[FPInd+3][1]) * (xa[FPInd+1] - x);
				printf("FPInd=%d, ProgPos=%f, ya[FPInd+1]=%f, diff[FPInd+2][0]=%f, diff[FPInd+3][1]=%f, x=%f, xa[FPInd+1]=%f\n", 
					  FPInd, ProgPos, ya[FPInd+1], diff[FPInd+2][0], diff[FPInd+3][1], x, xa[FPInd+1]);
			}
		}
		else ProgPos = 0.0;
		for (i=0, n=0; n<N; n++)
		{
		    //printf("Iso=%d, va[%d]=%d, ProgPos=%f\n", Iso, n, va[n], ProgPos);
		    if (r != 1) aJss = j - FoundProg[20][n]->Js + FoundProg[20][n]->Jss;
	    	else if (fjd == 1) aJss = j;
	    	else if (n - 2 * (n / 2) == 0) aJss = j - 1;
	    	else aJss = j + 1;
	    	if (aJss >= 0 && aJss <= Jeu)
	    	{
				if (ProgPos != 0.0) 
				{
		    		MinSP = ProgPos - ELU[Iso][va[n]][aJss] - 1.0;
		    		MaxSP = MinSP + 2.0;
				}
				else if (r == 1)
				{
					MinSP = ya[FPInd] - ELU[Iso][va[n]][aJss] 
							- (aJss <= j && FPInd > 0 ? ((d = ya[FPInd] - ya[FPInd-1]) <= 4.5 ?
							d - 0.5 : 5.0) : 5.0);
		    		MaxSP = ya[FPInd] - ELU[Iso][va[n]][aJss]
							+ (aJss >= j && FPInd < 40 ? ((d = ya[FPInd+1] - ya[FPInd]) <= 4.5 ?
							d - 0.5 : 5.0) : 5.0);
					printf("j=%d, aJss=%d\n", j, aJss);
					printf("ya[%d]=%f, ya[%d]=%f, ya[%d]=%f\n", 
						   FPInd-1, ya[FPInd-1], FPInd, ya[FPInd], FPInd+1, ya[FPInd+1]);
					printf("d=%f, MinSP=%f, MaxSP=%f\n", d, MinSP, MaxSP);
				}
				else 
				{
		    		MinSP = FoundProg[20][n]->Line[0] - 10.0;
		    		MaxSP = FoundProg[20][n]->Line[0] + 10.0;
				}
				//printf("Vor 1. while\n");
				while(marker[i].Line[0] > MinSP && i > 0) i--;
				while(marker[i].Line[0] < MinSP && i < AnzahlMarker - 1) i++;
				TEOCount[n] = 0;
				//printf("Vor 2. while\n");
				while(marker[i].Line[0] < MaxSP && i < AnzahlMarker - 1) 
				{
		    		if (r==0 && j==37) 
						printf("marker[%d].Line[0]=%f, HFLM=%f\n", 
							   i, marker[i].Line[0], marker[i].HFLM);
					if (marker[i].HFLM > MinPeakHeight)
		    		{
						//printf("AL=%d, TEOGC=%d\n", AL, TEOGC);
						TEOCount[n]++;
						TEO[TEOGC] = marker[i].Line[0] + ELU[Iso][va[n]][aJss];
						TEOPos[TEOGC] = n;
						Jssa[TEOGC] = aJss;
						TEOMarker[TEOGC++] = marker + i;
						//printf("aJss=%d\n", aJss);
		    		}
		    		i++;
				}
				D[n] = 0;
	    	}
		}
		//printf("Nach dem Falten\n");
		iBuff = 1;
		while (iBuff != -1) 
		{
	    	iBuff = -1;
	    	for (i=1; i<TEOGC; i++) if (TEO[i-1] > TEO[i])
	    	{
				DBuff1 = TEO[i-1];
				TEO[i-1] = TEO[i];
				TEO[i] = DBuff1;
				iBuff = TEOPos[i-1];
				TEOPos[i-1] = TEOPos[i];
				TEOPos[i] = iBuff;
				iBuff = Jssa[i-1];
				Jssa[i-1] = Jssa[i];
				Jssa[i] = iBuff;
				MBuff = TEOMarker[i-1];
				TEOMarker[i-1] = TEOMarker[i];
				TEOMarker[i] = MBuff;
	    	}
		}
		//printf("Nach dem Sortieren\n");
		BInd = BDens = Dens = k = 0;
		IntensB = Intens = XSum = 0.0;
		for (i=0; i<TEOGC; i++)
		{
	    	D[TEOPos[i]]++;
	    	if (D[TEOPos[i]] == 1) Dens++;
	    	Intens += TEOMarker[i]->HFLM;
	    	SPunkt = (XSum += TEO[i]) / (i - k + 1);
	    	while (SPunkt - TEO[k] > ST || TEO[i] - SPunkt > ST) 
	    	{
				D[TEOPos[k]]--;
				if (D[TEOPos[k]] == 0) Dens--;
				Intens -= TEOMarker[k]->HFLM;
				XSum -= TEO[k++];
				SPunkt = XSum / (i - k + 1);
	    	}
	    	if (j==36 && r==1) printf("i: TEO[%d]=%f, k: TEO[%d]=%f\n", i, TEO[i], k, TEO[k]);
			AIntens = Intens / (i - k + 1);
	    	if (Dens > BDens || (Dens == BDens && AIntens > IntensB))
	    	{
				BDens = Dens;
				BInd = i;
				IntensB = AIntens;
				Bk = k;
	    	}
		}
		//printf("Nach der Progressionssuche\n");
		DBuff1 = 0.0;
		for (i=Bk; i<=BInd; i++) DBuff1 += TEO[i];
		DBuff1 /= BInd - Bk + 1;
		for (n=0; n<N; n++) FoundProg[FPInd][n] = NULL;
		//printf("j=%d, N=%d, BDens=%d, NH=%d, IntensB=%f, MinPeakHight=%f, DBuff1=%f, ProgPos=%f, ST = %f, TEOGC=%d\n",
			//   j, N, BDens, NH, IntensB, MinPeakHight, DBuff1, ProgPos, ST, TEOGC);
		if (BDens >= 4 || (fabs(DBuff1 - ProgPos) < ST && ProgPos != 0.0))
		{	    
		    F = true;
		    xa[FPInd] = x;
	    	if (FPInd < 20) st = FPInd;
	    	else en = FPInd;
	    	y = 0;
	    	for (i=Bk; i<=BInd; i++)
	    	{
				FoundProg[FPInd][TEOPos[i]] = TEOMarker[i];
				if (!(TEOMarker[i]->Marked && TEOMarker[i]->DisplayData))
				{
		    		TEOMarker[i]->DD = DBuff1 - TEO[i];
		    		TEOMarker[i]->Iso = Iso;
		    		TEOMarker[i]->vss = va[TEOPos[i]];
		    		TEOMarker[i]->Jss = Jssa[i];
		    		TEOMarker[i]->Js = j;
		    		TEOMarker[i]->vs = vs;
		    		TEOMarker[i]->Marked = true;
		    		TEOMarker[i]->DisplayData = true;
					TEOMarker[i]->satellite = true;
					//TEOMarker[i]->SNR = TEOMarker[i]->HFLM / Rauschen;
					TEOMarker[i]->LState = LState;
					TEOMarker[i]->UState = UState;
					TEOMarker[i]->Mol = Mol;
					TEOMarker[i]->overlap = TEOMarker[i]->sOverlap;
					TEOMarker[i]->FC = -1;
					TEOMarker[i]->uncertainty = (TEOMarker[i]->overlap ? 0.02 : 0.005);
					TEOMarker[i]->added = false;
				}
				y += TEOMarker[i]->Line[0] + ELU[Iso][va[TEOPos[i]]][Jssa[i]];
				//if (j==14) 
					//printf("Jss[%d]=%d, E=%f\n", i, TEOMarker[i]->Jss, TEOMarker[i]->Line[0]);
	    	}
	    	ya[FPInd] = y / (BInd - Bk + 1);
	    	if (FPInd > 20) 
	    	{
				diff[FPInd][0] = (ya[FPInd] - ya[FPInd-1]) / (xa[FPInd] - xa[FPInd-1]);
				if (diff[FPInd-1][0] != 0.0) diff[FPInd][1] = diff[FPInd][0] - diff[FPInd-1][0];
				//printf("diff[%d][0]=%f, diff[%d][1]=%f\n",
		       		//				FPInd, diff[FPInd][0], FPInd, diff[FPInd][1]);
	    	}
	    	else if (FPInd < 20)
	    	{
				diff[FPInd+1][0] = (ya[FPInd+1] - ya[FPInd]) / (xa[FPInd+1] - xa[FPInd]);
				if (diff[FPInd+2][0] != 0.0) 
					diff[FPInd+2][1] = diff[FPInd+2][0] - diff[FPInd+1][0];
				//printf("diff[%d][0]=%f, diff[%d][1]=%f\n",
		       		//			FPInd+1, diff[FPInd+1][0], FPInd+2, diff[FPInd+2][1]);
	    	}
		}
		else if (r != 1 || FPInd != 20) F = false;
		if ((di == 1 && !F ) || FPInd == 40 || j == Jeu - 1)
		{
		    di = -1;
		    j = MJ - 1;
	    	F = true;
		}
		else j += di;
    }
    //printf("Ende FindSatelites\n");
	addMarker(false);
	Paint();
}

int Spektrum::SML(IntProg &P, int &I, int &J, double &T, double &EL, int veu, double ***ELU)
{
	//printf("SML\n");
	int n, v, N = P.N;
	double E, d;
	//if (I==0 && J >= 31 && J <= 33) printf("Beginn SML\n");
	if (T < ELU[I][0][J] || ELU[I][1][J] < ELU[I][0][J]) return 0;
	for (v = veu; (v > 0 ? T <= ELU[I][v][J] || ELU[I][v][J] < ELU[I][v-1][J] || ELU[I][v][J] == 0.0
			: false); v--) ;
	//printf("ST=%f\n", ST);
	for (n=1; v >= 0; v--)
	{
		E = T - ELU[I][v][J];
		while (GMarker[n]->Line[0] < E) n++;
		if (GMarker[n]->Line[0] - E > E - GMarker[n-1]->Line[0]) n--;
		if ((fabs(d = GMarker[n]->Line[0] - E)) < EL) 
		{
			P.vss[P.N] = v;
			P.Jss[P.N] = J;
			P.overlap[P.N] = false;
			P.DD[P.N] = d;
			P.marker[P.N++] = GMarker[n];
		}
		/*if (I==0 && J == 198) {printf("v=%d, d=%f\n", v, d);
			printf("ELU[%d][%d][%d]=%f, GMarker[%d]=%f, T=%f\n", 
				   I, v, J, ELU[I][v][J], n, GMarker[n]->Line[0], T);
			printf("E=%f, GMarker[%d]=%f, GMarker[%d]=%f\n", 
				   E, n-1, GMarker[n-1]->Line[0], n+1, GMarker[n+1]->Line[0]);}*/
	}
	//printf("Ende SML\n");
	return P.N - N;
}

void Spektrum::FindEmissionLines(Molecule* mol, int Iso, ElState* EState, double MinH, double Tol, int MLPB)
{
	ElState *XState = mol->getStateP(0);
	TermTable *XTerm = XState->getTermTable(), *ETerm = EState->getTermTable();
	double ****XData = XTerm->getData(), ****EData = ETerm->getData(), freq;
	int *CompZ = XTerm->getCompZ(), MXv = XTerm->getMaxv(), MJ = ETerm->getMaxJ(), NEC = ETerm->getNumComp();
	int vs, vss, MEv = ETerm->getMaxv(), J = XTerm->getMaxJ(), F, ec, fc, xc, NXC = XTerm->getNumComp();
	int NFC = int(2.0 * EState->getS()) + 1, i=0, n, N, eJs = EState->getJStart(Iso, 0), fJs = EState->getJStart(Iso, 1);
	int JStep = mol->getJStep(Iso), *Js = new int[3 * MJ], *Jss = new int[3 * MJ], M = MLPB / 2, m;
	double *DD = new double[3 * MJ];
	IsoTab *IsoT = mol->getIso();
	bool QLines = EState->getLambda() > 0, SL, LSL;
	Marker **Band = new Marker*[3 * MJ];
	QString XName = XState->getName(), EName = EState->getName();
	if (NXC < NFC && NEC < NFC) NFC = 1;
	if (J <= MJ) MJ = J - 1;
	for (F = 0; F < NFC; F++)
	{
		xc = (NXC >= NFC ? F : 0);
		ec = (NEC >= 2 * NFC && QLines ? 2 * F : (NEC >= NFC ? F : 0));
		fc = (NEC >= 2 * NFC ? ec + 1 : ec);
		for (vs = 0; vs <= MEv; vs++) for (vss = 0; vss < MXv; vss++)
		{
			for (N=0, J = (eJs > 0 ? eJs : eJs + JStep); J <= MJ; J++)
			{
				for (freq = EData[ec][Iso][vs][J] - XData[xc][Iso][vss][J-1]; i < AnzahlMarker - 1 && marker[i].Line[0] < freq; i++) ;
				while (i > 0 && marker[i].Line[0] > freq) i--;
				if (i+1 < AnzahlMarker ? marker[i+1].Line[0] - freq < freq - marker[i].Line[0] : false) i++;
				if ((DD[N] = fabs(freq - marker[i].Line[0])) <= Tol && marker[i].HFLM > MinH)
				{
					Js[N] = J;
					Jss[N] = J-1;
					//printf("vs=%d, vss=%d, Js=%d, Jss=%d\n", vs, vss, Js[N], Jss[N]);
					Band[N++] = marker + i;
				}
			}
			for (J = eJs; J < MJ; J++)
			{
				for (freq = EData[ec][Iso][vs][J] - XData[xc][Iso][vss][J+1]; i < AnzahlMarker - 1 && marker[i].Line[0] < freq; i++) ;
				while (i > 0 && marker[i].Line[0] > freq) i--;
				if (i+1 < AnzahlMarker ? marker[i+1].Line[0] - freq < freq - marker[i].Line[0] : false) i++;
				if ((DD[N] = fabs(freq - marker[i].Line[0])) <= Tol && marker[i].HFLM > MinH)
				{
					Js[N] = J;
					Jss[N] = J+1;
					//printf("vs=%d, vss=%d, Js=%d, Jss=%d\n", vs, vss, Js[N], Jss[N]);
					Band[N++] = marker + i;
				}
			}
			if (QLines) for (J = fJs; J < MJ; J++)
			{
				for (freq = EData[fc][Iso][vs][J] - XData[xc][Iso][vss][J]; i < AnzahlMarker - 1 && marker[i].Line[0] < freq; i++) ;
				while (i > 0 && marker[i].Line[0] > freq) i--;
				if (i+1 < AnzahlMarker ? marker[i+1].Line[0] - freq < freq - marker[i].Line[0] : false) i++;
				if ((DD[N] = fabs(freq - marker[i].Line[0])) <= Tol && marker[i].HFLM > MinH)
				{
					Js[N] = Jss[N] = J;
					Band[N++] = marker + i;
				}
			}
			if (N < MLPB) continue;
			for (n=1, m=0, LSL = SL = false; n<=N; n++) 
			{
				if (n<N ? Js[n] <= Js[n-1] + 2 * JStep && Js[n] > Js[n-1] : false) m++;
				else
				{
					if (m>=M) LSL = SL = true;
					else if (n<N ? Js[n] - Js[n-1] > 2 * JStep || Js[n] < Js[n-1] : true)
					{
						if (!LSL) 
						{
							for (m=n-1; (m>0 ? Js[m] - Js[m-1] <= 2 * JStep && Js[m] > Js[m-1] : false); m--) 
								Jss[m] = -1;
							Jss[m] = -1;
						}
						else LSL = false;
					}
					m=0;
 				}
			}
			if (SL) for (n=0; n<N; n++)
			{
				if (Jss[n] == -1) continue;
				if (Band[n]->DisplayData)
				{
					if (Band[n]->Mol == mol && Band[n]->FC == F && Band[n]->Iso == Iso && Band[n]->UState == EState
						&& Band[n]->Js == Js[n] && Band[n]->Jss == Jss[n] && Band[n]->vs == vs && Band[n]->vss == vss) continue;
					else Band[n]->overlap = true;
				}
				else Band[n]->overlap = Band[n]->sOverlap;
				Band[n]->DD = DD[n];
				Band[n]->DisplayData = true;
				Band[n]->FC = CompZ[F];
				Band[n]->Iso = Iso;
				Band[n]->IsoName = IsoT->texName[Iso];
				//printf("2. vs=%d, vss=%d, Js=%d, Jss=%d\n", vs, vss, Js[n], Jss[n]);
				Band[n]->Js = Js[n];
				Band[n]->Jss = Jss[n];
				Band[n]->LState = XState;
				Band[n]->lState = XName;
				Band[n]->Marked = true;
				Band[n]->Mol = mol;
				Band[n]->UState = EState;
				Band[n]->uState = EName;
				Band[n]->satellite = false;
				Band[n]->uncertainty = (Band[n]->overlap ? 0.02 : 0.005);
				Band[n]->vs = vs;
				Band[n]->vss = vss;
			}
		}
	}
	delete IsoT;
	delete[] Js;
	delete[] Jss;
	delete[] DD;
	delete[] Band;
	FocusSelected();
}

double Spektrum::getAUT(IntProg &P)
{
	if (P.N == 0) return 0.0;
	int n;
	double T = 0.0, ***ELU = P.lState->getTermTable()->getData()[0];
	for (n=0; n<P.N; n++) T += P.marker[n]->Line[0] + ELU[P.Iso][P.vss[n]][P.Jss[n]];
	return T / P.N;
}

void Spektrum::continueBand()
{
	if (band == 0)
	{
		printf("Spektrum::continueBand error: band=0!\n");
		return;
	}
	/*int n, N = band->MaxI - band->MinI + 1, SDir = -1, i;
	if (N==1) return;
	ElState *LState = band->lines[band->CL1]->LState, *UState = band->lines[band->CL1]->UState; 
	//double ***ELU = (LState != 0 ? LState->
	double Y[6], PPos, PI, NC, FID, SID;
	for (i = band->MinI + SDir; SDir != 1 || band->MaxI < 999; i += SDir)
	{
		if (N==2)
		{
			NC = 2;
			
			
		}
	}*/
}

void Spektrum::ContinueProgressions()
{
	IntProg Prog;
	double ****lTerm, ****uTerm, UT[1000], AT;
	int nC, nIso, nv, nJ, n, m, unv, unIso, unJ, c;
	for (n=0; n < AnzahlMarker; n++) marker[n].b = false;
	for (n=0; n < AnzahlMarker; n++) 
		if (marker[n].Marked && marker[n].DisplayData && marker[n].LState != 0 && !marker[n].b)
	{
		Prog.clear();
		Prog.Iso = marker[n].Iso;
		Prog.Js = marker[n].Js;
		Prog.vs = marker[n].vs;
		Prog.FC = marker[n].FC;
		Prog.lState = marker[n].LState;
		Prog.uState = marker[n].UState;
		Prog.Mol = marker[n].Mol;
		Prog.Mol->getTermData(marker[n].LState->getStateNum(), nC, nIso, nv, nJ, lTerm);
		if (lTerm == 0 || Prog.Iso < 0 || Prog.Iso >= nIso || marker[n].vss < 0 || marker[n].vss >= nv
		   || marker[n].Jss < 0 || marker[n].Jss >= nJ) continue;
		Prog.uT = lTerm[0][Prog.Iso][marker[n].vss][marker[n].Jss] + marker[n].Line[0];
		Prog.Np = SML(Prog, Prog.Iso, marker[n].Jss, Prog.uT, ST, nv - 1, lTerm[0]);
		if (marker[n].Jss != Prog.Js && (m = 2.0 * Prog.Js - marker[n].Jss) >= 0)
		{
			SML(Prog, Prog.Iso, m, Prog.uT, ST, nv - 1, lTerm[0]);
			chkDoubletts(Prog);
		}
		for (m=0, AT = 0.0; m < Prog.N; m++) AT += UT[m] = lTerm[0][Prog.Iso][Prog.vss[m]][Prog.Jss[m]];
		if (Prog.uState != 0)
		{
			Prog.Mol->getTermData(Prog.uState->getStateNum(), nC, unIso, unv, unJ, uTerm);
			c = (Prog.Js == marker[n].Jss && nC > 1 ? 1 : 0);
			if (uTerm != 0 && Prog.Iso < unIso && Prog.vs < unv && Prog.Js < unJ ? 
						 uTerm[c][Prog.Iso][Prog.vs][Prog.Js] != 0.0 : false)
				for (m=0; m < Prog.N; m++) Prog.marker[m]->oc = UT[m] - uTerm[c][Prog.Iso][Prog.vs][Prog.Js];
			else for (m=0; m < Prog.N; m++) Prog.marker[m]->oc = 0.0;
		}
		else for (m=0; m < Prog.N; m++) Prog.marker[m]->oc = 0.0;
		for (m=0, AT /= Prog.N; m < Prog.N; m++)
		{
			Prog.marker[m]->b = true;
			Prog.marker[m]->DD = AT - UT[m];
			Prog.marker[m]->Iso = Prog.Iso;
			Prog.marker[m]->IsoName = marker[n].IsoName;
			Prog.marker[m]->Js = Prog.Js;
			Prog.marker[m]->Jss = Prog.Jss[m];
			Prog.marker[m]->lState = marker[n].lState;
			Prog.marker[m]->LState = Prog.lState;
			Prog.marker[m]->Mol = Prog.Mol;
			if ((Prog.marker[m]->DisplayData && Prog.marker[m] != marker + n) || Prog.marker[m]->sOverlap)
				Prog.marker[m]->overlap = true;
			else Prog.marker[m]->overlap = false;
			Prog.marker[m]->Marked = true;
			Prog.marker[m]->DisplayData = true;
			Prog.marker[m]->satellite = marker[n].satellite;
			//Prog.marker[m]->SNR = Prog.marker[m]->HFLM / Rauschen;
			Prog.marker[m]->uState = marker[n].uState;
			Prog.marker[m]->UState = Prog.uState;
			Prog.marker[m]->vs = Prog.vs;
			Prog.marker[m]->vss = Prog.vss[m];
			Prog.marker[m]->FC = Prog.FC;
			Prog.marker[m]->uncertainty = (Prog.marker[m]->overlap ? 0.02 : 0.005);
		}
	}
	addMarker(false);
	if (MW != 0) MW->setActive(this);
	Paint();
}

void Spektrum::AssignBands()
{
	if (AssignmentStatus == 1) AssignmentStatus = 0;
	else AssignmentStatus = 1;
}

void Spektrum::assignBandByDoubletPartners()
{
	int n, l, N, m, NMol, bMol = 0, I, NI, bI = 0, v, bv = 0, Nv, bJ = 0, J, cJ, JD, LT, bLT = 0, bJD = 0, NJ, c, bc = 0, Nc, bDN = 0, bLN = 0, DN, LN;
	int JStart, JEnd, JStep;
	double D, bD = 0.0, ****XData, E;
	Molecule *Mol;
	for (n=N=0; n < AnzahlMarker; n++) if (marker[n].Marked && marker[n].DisplayData) N++;
	if (N==0)
	{
		QMessageBox::information(this, "MolSpektAnalysis", 
								 "Error: There are no lines marked which belong to neighboring rotational levels of a band.");
		return;
	}
	IntProg *Progs = new IntProg[N];
	for (n=N=0; n < AnzahlMarker; n++) if (marker[n].Marked && marker[n].DisplayData) Progs[N++].marker[0] = marker + n;
	for (m=0, NMol = MW->getNumMolecules(); m < NMol; m++)
	{
		Mol = MW->getMolecule(m);
		Mol->getTermData(0, Nc, NI, Nv, NJ, XData);
		for (c=0; c < Nc; c++) for (I=0; I < NI; I++) 
			for (v=0, JStep = Mol->getJStep(I), JStart = Mol->getStateP(0)->getJStart(I, 0); v < Nv; v++) 
			for (JD = -JStep; JD <= JStep; JD += 2 * JStep) for (LT = -2; LT <= 2; LT += 4)
		{
			JEnd = NJ - (LT > 0 ? LT : 0) - (JD > 0 ? (N-1) * JStep : 0);
			for (J = JStart - (LT < 0 ? LT : 0) - (JD < 0 ? (N-1) * JD : 0); J < JEnd; J += JStep)
			{
				for (n = DN = LN = 0, cJ = J, D = 0.0; n<N; n++, cJ += JD)
				{
					Progs[n].N = 1;
					LN += SML(Progs[n], I, l = cJ + LT, E = XData[c][I][v][cJ] + Progs[n].marker[0]->Line[0], ST, Nv - 1, XData[c]);
					for (l=1; l < Progs[n].N; l++) if (Progs[n].vss[l] == v)
					{
						DN++;
						D += Progs[n].DD[l] * Progs[n].DD[l];
					}
				}
				if (DN > bDN || (DN == bDN && (LN > bLN || (LN == bLN && D < bD))))
				{
					bDN = DN;
					bLN = LN;
					bD = D;
					bMol = m;
					bI = I;
					bc = c;
					bv = v;
					bJ = J;
					bJD = JD;
					bLT = LT;
				}
			}
		}
	}
	Mol = MW->getMolecule(bMol);
	IsoTab *Iso = Mol->getIso();
	for (n=0, J = bJ; n<N; n++, J += bJD)
	{
		Progs[n].marker[0]->FC = bc;
		Progs[n].marker[0]->Iso = bI;
		Progs[n].marker[0]->IsoName = Iso->texName[bI];
		Progs[n].marker[0]->Js = J + bLT / 2;
		Progs[n].marker[0]->Jss = J;
		Progs[n].marker[0]->LState = Mol->getStateP(0);
		Progs[n].marker[0]->lState = Progs[n].marker[0]->LState->getName();
		Progs[n].marker[0]->Mol = Mol;
		Progs[n].marker[0]->satellite = false;
		Progs[n].marker[0]->uState = "";
		Progs[n].marker[0]->UState = 0;
		Progs[n].marker[0]->vs = -1;
		Progs[n].marker[0]->vss = bv;
	}
	delete Iso;
	delete[] Progs;
	ContinueProgressions();
}

bool Spektrum::excitedLevel(IntProg &P)
{
	int s, S = P.Mol->getNumStates(), v, mv, c, NC, NIso, mJ;
	double d, dev = 5.0, ****Data;
	for (s=1; s<S; s++)
	{
		P.Mol->getTermData(S, NC, NIso, mv, mJ, Data);
		if (NC == 0 || NIso <= P.Iso || P.Js > mJ) continue;
		for (v=0, c = (P.Js == P.Jss[0] && NC > 1 ? 1 : 0); 
				(v <= mv ? Data[c][P.Iso][v][P.Js] - dev < P.uT : false); v++)
			if (fabs(d = P.uT - Data[c][P.Iso][v][P.Js]) < dev)
		{
			dev = fabs(P.oc[0] = d);
			P.uState = P.Mol->getStateP(s);
			P.vs = v;
		}
	}
	return (dev < 5.0 ? true : false);
}

int Spektrum::Assignvs(IntProg &P)
{
	int n, i, M, v, NC, NI, mvs, mJs, S, NS, bv = -1, bS = 0, R=2, JB;
	double ***T, ****Data, bmoc = 1e20, IS;
	bool FCFC = false;
	Molecule *Mol = P.Mol;
	FCFTab *FCFT;
	if (Mol == 0)
	{
		for (n=0; n <= P.N; n++) P.oc[n] = 0.0;
		P.vs = -1;
		return 2;
	}
	Mol->getTermData(0, NC, NI, mvs, mJs, Data);
	//Potential *uPot;
	double ***ELU = Data[0], WF;
	NS = Mol->getNumStates();
	double **OC = Create(NS, P.N), moc[NS], FCFf[NS], **FCF = Create(NS, M = 2 * (JB = mvs + 1));
	int J = (P.Js == P.Jss[0] ? P.Js : P.Js - 1), s[M];
	for (n=0, IS = 0.0; n < P.N; n++) IS += P.marker[n]->HFLM;
	for (S=1; S < NS; S++)
	{
		moc[S] = 0.0;
		for (n=0; n<M; n++) FCF[S][n] = 0.0;
		//uPot = Mol->getPot(S);
		Mol->getTermData(S, NC, NI, mvs, mJs, Data);
		if (NC == 0) continue;
		if (Mol->getStateP(S)->getLambda() == 0 && P.Js == P.Jss[0]) continue;
		T = (P.Js != P.Jss[0] || NC == 1 ? Data[0] : Data[1]);
		for (v=0; (v < mvs ? T[P.Iso][v][P.Js] - ELU[P.Iso][P.vss[0]][P.Jss[0]] 
 				- P.marker[0]->Line[0] < 0 : false); v++) ;
		if ((v > 0 ? (v < mvs ? T[P.Iso][v][P.Js] - ELU[P.Iso][P.vss[0]][P.Jss[0]] 
				- P.marker[0]->Line[0] > P.marker[0]->Line[0] - T[P.Iso][v-1][P.Js] 
				+ ELU[P.Iso][P.vss[0]][P.Jss[0]] : true) : false)) v--;
		//if (fabs(T[P.Iso][v][P.Js] - ELU[P.Iso][P.vss[0]][P.Jss[0]] - P.marker[0]->Line[0]) < 5.0)
		for (n=0; n < P.N; n++)
			moc[S] += fabs(OC[S][n] = P.marker[n]->Line[0] - T[P.Iso][v][P.Js] 
									+ ELU[P.Iso][P.vss[n]][P.Jss[n]]);
		if (((P.Js == P.Jss[0] ? P.G < 5 : P.NGD < 2) ? moc[S] < P.N 
			: moc[S] < 5 * P.N) ? (FCFT = Mol->getFCF(0, S)) != 0 : false)
		{
			/*for (m=0; J <= P.Js + 1; J+=2, JB = m)
				for (n=0; (n <= mvs ? ELU[P.Iso][n][J] != 0 : false); n++)
			{
				Ja[m] = J;
				En[m++] = ELU[P.Iso][n][J];
			}
			lPot->getFastFCF(P.Iso, uPot, P.Js, T[P.Iso][v][P.Js], m, En, Ja, FCF[S]);*/
			FCFT->getFCF(P.Js - P.Jss[0], P.Iso, P.Js, v, JB, FCF[S]);
			if (P.Js == P.Jss[0] ? P.G < 5 : P.NGD < 2)
			{
				for (n=0; n<M; n++) s[n] = n;
				for (i=0; i!=-1; ) for (i=-1, n=1; n<M; n++) if (FCF[S][s[n]] > FCF[S][s[n-1]])
				{
					i = s[n];
					s[n] = s[n-1];
					s[n-1] = i;
				}
				for (n=0; n < P.N && moc[S] != 1e50; n++) 
				{
					if (s[n] < JB)
					{
						v = s[n];
						J = (P.Js == P.Jss[0] ? P.Js : P.Js - 1);
					}
					else
					{
						v = s[n] - JB;
						J = P.Js + 1;
					}
					for (i=0; (i < P.N ? v != P.vss[i] || J != P.Jss[i] : false); i++) ;
					if (i == P.N) moc[S] = 1e50;
				}
			}
			if (moc[S] < 1e50)
			{
				if (!FCFC) FCFC = true;
				for (n=i=0, WF = 0.0; n < P.N; n++) 
				{
					s[v = (P.Jss[n] <= P.Js ? P.vss[n] : P.vss[n] + JB)] = n;
					if (!P.overlap[n]) 
					{
						WF += FCF[S][v] / P.marker[n]->HFLM;
						i++;
					}
				}
				WF /= P.G;
				for (n=0, FCFf[S] = 0.0; n < M; n++)
				{
					FCF[S][n] /= WF;
					FCFf[S] += (s[n] >= 0 ? fabs(P.marker[s[n]]->HFLM - FCF[S][n]) 
											: (FCF[S][n] > 2.0 * MinPeakHeight ? FCF[S][n] : 0.0));
				}
				moc[S] *= FCFf[S];
			}
		}
		else moc[S] *= (IS * 0.5);
		if (moc[S] < bmoc)
		{
			bmoc = moc[S];
			bv = v;
			bS = S;
		}
	}
	if (P.Js == P.Jss[0] ? P.G < 5 : P.NGD < 2)
	{
		if (FCFC && FCFf[bS] < 0.25 * IS) R=1;
		else R=2;
	}
	else
	{
		if (FCFC && FCFf[bS] < 0.25 * IS) R=0;
		else R=1;
	}
	if (bmoc < 10 * P.N)
	{
		P.uState = Mol->getStateP(bS);
		P.vs = bv;
	}
	else
	{
		bmoc = double(10 * P.N);
		P.uState = 0;
		P.vs = -1;
	}
	P.QF = P.N * IS * (P.NGD + P.G) / bmoc;
	for (n=0; n < P.N; n++) 
	{
		P.oc[n] = OC[bS][n];
		if (P.overlap[n] && FCF[bS][n] > 2.0 * MinPeakHeight) P.marker[n] = 0;
	}
	UpdateProg(P);
	Destroy(OC, NS);
	Destroy(FCF, NS);
	return R;
}

void Spektrum::cFQS(IntProg &P)
{
	int n;
	P.FQS = 0.0;
	for (n=0; n < P.N; n++) //{
		P.FQS += sqr(P.DD[n]);
		//printf("P.DD[%d]=%f\n", n, P.DD[n]);}
}

void Spektrum::chkDoubletts(IntProg &P)
{
	//printf("chkDoubletts\n");
	int n, m=P.Np;
	double r;
	P.NGD = 0;
	//if (P.Js==197) printf("Js=197, N=%d\n", N);
	for (n=0; n < P.Np; n++)
	{
		if (m < P.N && n > 0) if (P.vss[n-1] == P.vss[m]) m++;
		for (; (m < P.N ? P.vss[n] < P.vss[m] : false); m++) 
			if (P.overlap[m] || 
						 (((m > P.Np ? P.vss[m-1] > P.vss[m] + 1 || P.marker[m-1] == NULL : true) 
										|| (m+1 < P.N ? P.vss[m+1] < P.vss[m] - 1 : true))
										&& P.marker[m] != LaserLinie)) 
		{
			//if (P.Js==113) 
				//printf("Js=113, Lösche %f mit v''=%d bei 0\n", P.marker[m]->Line[0], P.vss[m]);
			P.marker[m] = NULL; 
		}
		if ((m < P.N ? P.vss[n] > P.vss[m] : true))
		{
			if (P.overlap[n] || (((n > 0 ? P.vss[n-1] > P.vss[n] + 1 || P.marker[n-1] == NULL : true) 
						 || (n+1 < P.Np ? P.vss[n+1] < P.vss[n] - 1 : true))
								&& P.marker[n] != LaserLinie)) 
			{
				/*if (P.Js==113)
				{
					printf("Js=113, Lösche %f mit v''=%d bei 1\n", P.marker[n]->Line[0], P.vss[n]);
				printf("m=%d, n=%d, P.N=%d, vss[n]=%d, vss[m]=%d\n", m, n, P.N, P.vss[n], P.vss[m]);
				}*/
				P.marker[n] = NULL;
			}
			continue;
		}
		/*if (P.Js==79) 
		{
			printf("Js=79: NGD=%d vollst. Doublett mit %f und %f\n", P.NGD, 
				   P.marker[n]->Line[0], P.marker[m]->Line[0]);
			if (P.marker[n]->sOverlap) printf("sOverlap bei %f\n", P.marker[n]->Line[0]);
			if (P.marker[m]->sOverlap) printf("sOverlap bei %f\n", P.marker[m]->Line[0]);
		}*/
		if (!P.overlap[n] && !P.overlap[m])
		{
			r = (Rauschen + P.marker[m]->HFLM) / (P.marker[n]->HFLM + Rauschen);
		//if (P.Js==79) printf("Js=79: kein ov. mit Doublett mit %f, r=%f\n", P.marker[n]->Line[0], r);
			if (r > 0.75 && r < 1.33) P.NGD++;
		}
		else if (P.overlap[m] && P.overlap[n] && (((m > P.Np ? P.vss[m-1] > P.vss[m] + 1 : true) 
				  && (m+1 < P.N ? P.vss[m+1] < P.vss[m] - 1 : true)) 
				  || ((n > 0 ? P.vss[n-1] > P.vss[n] + 1 : true) 
				  && (n+1 < P.Np ? P.vss[n+1] < P.vss[n] - 1 : true)))) 
		{
			/*if (P.Js==113) 
				printf("Js=113, Lösche %f und %f mit v''=%d und %d bei 2\n", P.marker[n]->Line[0], 
					   P.marker[m]->Line[0], P.vss[n], P.vss[m]);*/
			P.marker[n] = P.marker[m] = NULL; 
		}
	}
	if (m < P.N && P.Np > 0) if (P.vss[P.Np - 1] == P.vss[m]) m++;
	//if (P.Js==163) 
		//printf("N=%d, m=%d, P.vss[N-1]=%d, P.vss[m]=%d\n", P.Np, m, P.vss[P.Np - 1], P.vss[m]);
	for (; m < P.N; m++) 
		if (P.overlap[m] || (((m > P.Np ? P.vss[m-1] > P.vss[m] + 1 || P.marker[m-1] == NULL : true) 
			|| (m+1 < P.N ? P.vss[m+1] < P.vss[m] - 1 : true)) && P.marker[m] != LaserLinie)) 
	{
	//if (P.Js==113) printf("Js=113, Lösche %f mit v''=%d bei 4\n", P.marker[m]->Line[0], P.vss[m]);
		P.marker[m] = NULL; 
	}		
	UpdateProg(P);
	//printf("Ende chkDoubletts\n");
}

void Spektrum::UpdateProg(IntProg &P)
{
	int n, m;
	for (n=0; (n < P.N ? P.marker[n] != NULL : false); n++) ;
	if (n == P.Np - 1) P.Np--;
	for (m=n+1; m < P.N; m++) 
	{
		if (P.marker[m] != NULL) 
		{
			P.marker[n] = P.marker[m];
			P.vss[n] = P.vss[m];
			P.Jss[n] = P.Jss[m];
			P.overlap[n] = P.overlap[m];
			P.DD[n] = P.DD[m];
			P.SNR[n] = P.SNR[m];
			P.oc[n++] = P.oc[m];
		}
		if (m == P.Np - 1) P.Np = n;
	}
	P.N = n;
}

void Spektrum::eNGL(IntProg &P, Marker **&AL, int NAL, double ***ELU)
{
	//printf("eNGL\n");
	int n, m=0, i, j;
	double M;
	Marker *o;
	NAL--;
	P.G = 0;
	P.uT = 0.0;
	/*for (n=0; n < P.N; n++)
	{
		printf("P.N=%d, P.Iso=%d, P.Js=%d, P.Jss[%d]=%d, P.vss[%d]=%d\n", 
			P.N, P.Iso, P.Js, n, P.Jss[n], n, P.vss[n]);
		printf("ELU=%f\n", ELU[P.Iso][P.vss[n]][P.Jss[n]]);
		printf("AnzahlMarker=%d, Index=%d\n", AnzahlMarker, P.marker[n] - marker);
		printf("P.marker=%f\n", P.marker[n]->Line[0]);
	}*/
	for (n=0; n < P.N; n++) P.uT += ELU[P.Iso][P.vss[n]][P.Jss[n]] + P.marker[n]->Line[0];
	P.uT /= P.N;
	for (n=0; n < P.N; n++)
	{
		if (fabs(P.DD[n] = P.uT - ELU[P.Iso][P.vss[n]][P.Jss[n]] - P.marker[n]->Line[0]) > ST)
		{
			P.marker[n] = NULL;
			continue;
		} 
		//if (P.Iso==5 && P.Js==92) printf("n=%d, P.N=%d\n", n, P.N);
		o = (P.marker[n] >= marker + 2 ? (P.marker[n] < marker + AnzahlMarker - 2 ? P.marker[n] - 2 
			: marker + AnzahlMarker - 5) : marker);
		/*if (P.Iso==5 && P.Js==92) 
			printf("o-marker=%d, P.marker[n]=%d\n", o - marker, P.marker[n] - marker);
		//if (P.N==9 && P.Js==94 && NAL == 1068) printf("o[0].Line[1]=%f\n", o[0].Line[1]);*/
		if ((o[0].Line[1] < o[1].Line[1] && o[1].Line[1] < o[2].Line[1] 
				   && o[2].Line[1] < o[3].Line[1] && o[3].Line[1] < o[4].Line[1]) 
				   || (o[0].Line[1] > o[1].Line[1] && o[1].Line[1] > o[2].Line[1] 
				   && o[2].Line[1] > o[3].Line[1] && o[3].Line[1] > o[4].Line[1]))
		{
			P.marker[n] = NULL;
			continue;
		}
		for (i=j=0, M = 0.7 * P.marker[n]->HFLM; i<5; i++) 
			if (o[i].HFLM > M && !o[i].DisplayData) j++;
		if (j>2) P.overlap[n] = true;
		//if (P.Iso==5 && P.Js==92) printf("n=%d, P.N=%d\n", n, P.N);
		if (P.marker[n] == LaserLinie) P.overlap[n] = false;
		else if (P.marker[n]->sOverlap) P.overlap[n] = true;
		else if ((P.marker[n] > marker ? (P.marker[n] - 1)->Line[1] > P.marker[n]->Line[1] : false) 
						|| (P.marker[n] < marker + AnzahlMarker - 1 ? 
						(P.marker[n] + 1)->Line[1] > P.marker[n]->Line[1] : false)) 
			P.overlap[n] = true;
		else 
		{
			for (i=0, P.overlap[n] = false; i < P.N; i++) 
				if (P.marker[n] == P.marker[i] && n!=i) P.overlap[n] = true;
			if (!P.overlap[n] && NAL >= 0)
			{
				//if (P.N==9 && P.Js==94 && NAL == 1068) printf("Vor While\n");
				while ((m < NAL ? AL[m] < P.marker[n] : false)) m++;
				while ((m > 0 ? AL[m] > P.marker[n] : false)) m--;
				//if (P.Iso==5 && P.Js==92) printf("Nach While\n");
				if (AL[m] == P.marker[n]) P.overlap[n] = true;
			}
			if (!P.overlap[n]) { P.G++;
					//if (P.Iso==5 && P.Js==92) 
					//printf("J'=94, gute Linie: %f\n", P.marker[n]->Line[0]);
					}
			
		}
	}
	UpdateProg(P);
	//printf("eNGL Ende\n");
}

void Spektrum::getMarkedLines(Marker **&M, int &N)
{
	//printf("getMarkedLines\n");
	int n;
	if (M != NULL) delete[] M;
	for (N=0, n=0; n < AnzahlMarker; n++) if (marker[n].Marked) N++;
	if (N == 0)
	{
		M = NULL;
		return;
	}
	M = new Marker*[N];
	for (n=0, N=0; n < AnzahlMarker; n++) if (marker[n].Marked) M[N++] = marker + n;
	//printf("Ende getMarkedLines\n");
}

void Spektrum::BRauschen()
{
    /*int i, j, n, N;
    if (AnzahlMarker > 10000) N = 1000;
    else N = AnzahlMarker / 10;
    double B[N];
	int I[N];
	int f;
    for (i=0, n=N; i<AnzahlMarker; i++)
    {
		for (j=(n > 0 ? n : 1); j < N && marker[i].HFLM > B[j]; j++) {B[j - 1] = B[j]; I[j-1]=I[j];}
		if (j > 1) 
		{
	    	B[j-1] = marker[i].HFLM;
			I[j-1] = i;
	    	if (n > 0) n--;
		}
		else if (marker[i].HFLM > B[0]) 
		{
			B[0] = marker[i].HFLM;
			I[0] = i;
		}
    }
    Rauschen = B[0];
    f = I[0];
	printf("Rauschen=%f, freq=%f, RM=%f, LM=%f, i=%d\n", Rauschen, marker[f].Line[0], marker[f].RMin[0],
		    marker[f].LMin[0], f);
    MinPeakHight = 1.5 * Rauschen;*/
	printf("BRauschen\n");
	double Min, Max, SF;
	int i, j, P, ID[1000], m, M = AnzahlMarker / 50, R=0;
	int ML = int((Daten->GetValue(Daten->GetDSL() - 1, 0) - Daten->GetValue(0, 0)) / (100.0 * ST));
	for (i=1, Min=Max=marker[0].HFLM; i<AnzahlMarker; i++)
	{
		if (marker[i].HFLM < Min) Min = marker[i].HFLM;
		else if (marker[i].HFLM > Max && marker + i != LaserLinie) Max = marker[i].HFLM;
	}
	SF = 1000.0 / (Max - Min);
	for (m = AnzahlMarker; m > M; SF *= 2.0)
	{
		for (i=0; i<1000; i++) ID[i] = 0;
		for (i=0; i<AnzahlMarker; i++) 
		{
			P = int(SF * (marker[i].HFLM - Min));
			if (P < 1000 && P >= 0) ID[P]++;
		}
		for (i=m=0; i<1000; i++) if (ID[i] > m) m = ID[i];
		printf("SF=%f, m=%d, M=%d\n", SF, m, M);
	}
	for (M=0, i=5; i<995; i++) 
	{
		m = 0;
		for (j=i-5; j<i+5; j++) m += ID[j];
		//printf("i=%d, m=%d, R=%d, M=%d\n", i, m, R, M);
		if (m >= M)
		{
			M = m;
			R = i;
		}
	}
	Rauschen = Min + 1.4142 * (double(R) + 0.5) / SF;
	printf("Rauschen=%f, SF=%f, R=%d, M=%d\n", Rauschen, SF, R, M);
	for (i=999, M=0; i >= 0 && M < ML; i--) M += ID[i];
	MinPeakHeight = Min + (double(i) + 0.5) / SF;
	if (M > 0 && AnzahlMarker > ML)
	{
		j = ML / 10;
		while (fabs(ML - M) > j)
		{
			if (M > ML) MinPeakHeight += SF *= 0.67;
			else MinPeakHeight -= SF *= 0.67;
			for (i=M=0; i < AnzahlMarker; i++) if (marker[i].HFLM >= MinPeakHeight) M++;
			printf("ML=%d, M=%d, MinPeakHight=%f\n", ML, M, MinPeakHeight);
		}
		MinPeakHeight *= 1.6;
	}
	
	MinPeakHeight *= 1.5;
	
	GetMGMH();
	printf("MinPeakHeight=%f\n", MinPeakHeight);
	/*double ex, ex1, x, S, S1, N, N1, xs;
	for (S=S1=N=N1=x=0.0; x < 100.0; x+=0.0001)
	{
		N += x*(ex = exp(-1.0 * (xs = x * x)));
		N1 += x*(ex1 = exp(-.11111 * xs));
		S += ex * xs * x;
		S1 += ex1 * xs * x;
	}
	printf("R = %f, R1 = %f\n", sqrt(S / N), sqrt(S1/N1));*/
	printf("aGPpI=%f\n", 
		   double(M) * ST / (Daten->GetValue(Daten->GetDSL() - 1, 0) - Daten->GetValue(0, 0)));
}

bool Spektrum::cbSat(int &I, int &v, int &J, Marker &L, Marker **&A, int &N, Marker *&ML)
{
	int n, D;
	for (n=0; n<N; n++) if (A[n]->vss == v && A[n]->Iso == I)
	{
		D = A[n]->Jss - J;
		if (D >= -2 && D <= 2 && L.HFLM < A[n]->HFLM) 
		{
			ML = A[n];
			return true;
		}
	}
	return false;
}

void Spektrum::WriteIntensDist()
{
	if (Daten->GetDSL() == 0) return;
	QString F = QFileDialog::getSaveFileName(this,  "Choose a file name for the intensity distribution", "", "Data Files (*.dat)"); 
	if (F.isEmpty()) return;
	QFile Datei(F);
	if (!Datei.open(QIODevice::WriteOnly)) 
    {
		QString Fehlermeldung = "The file " + F + " can not be opened.";
		QMessageBox::information( this, "MolSpektAnalysis", Fehlermeldung, QMessageBox::Ok);
		return;
    }
	QTextStream S(&Datei);
	if (AnzahlMarker == 0) editFind();
	double Min, Max, SF, EDSF;
	int i, P, ID[100000];
	for (i=1, Min=Max=marker[0].HFLM; i<AnzahlMarker; i++)
	{
		if (marker[i].HFLM < Min) Min = marker[i].HFLM;
		else if (marker[i].HFLM > Max && marker + i != LaserLinie) Max = marker[i].HFLM;
	}
	SF = 100000.0 / (Max - Min);
	EDSF = 1.0 / SF;
	for (i=0; i<100000; i++) ID[i] = 0;
	for (i=0; i<AnzahlMarker; i++) 
	{
		P = int(SF * (marker[i].HFLM - Min));
		if (P < 100000 && P >= 0) ID[P]++;
	}
	for (i=0; i<100000; i++) 
		S << QString::number((double(i) + 0.5) * EDSF + Min, 'g', 7) << "	" << ID[i] << "\n";
}

int *Spektrum::getIntensDist(int N, double min, double max)
{
	double SF = double(N) / (max - min);
	int *R = new int[N], n, i;
	for (n=0; n<N; n++) R[n] = 0;
	if (AnzahlMarker == 0) editFind();
	for (n=0; n < AnzahlMarker; n++)
	{
		i = int(SF * (marker[n].HFLM - min));
		if (i >= 0 && i < N) R[i]++;
	}
	return R;
}

void Spektrum::getMinMaxIntensities(double &Min, double &Max)
{
	int n;
	if (AnzahlMarker == 0) editFind();
	for (n=1, Min = Max = marker[0].HFLM; n < AnzahlMarker; n++)
	{
		if (marker[n].HFLM < Min) Min = marker[n].HFLM;
		else if (marker[n].HFLM > Max) Max = marker[n].HFLM;
	}
}
	
void Spektrum::getRMPH(double &R, double &M)
{
	if (AnzahlMarker == 0) editFind();
	R = Rauschen;
	M = MinPeakHeight;
}

void Spektrum::getProgProfile(ElState *UState, ElState *LState, int Iso, int vs, int Js, int Comp, 
					    		int *&vss, int *&Jss, double *&I, int &N)
{
	int n, m;
	for (n=N=0; n < AnzahlMarker; n++) 
		if ((marker[n].UState == UState || UState == 0) && marker[n].LState == LState 
				   && marker[n].Iso == Iso && (marker[n].vs == vs || vs == -1) && marker[n].Js == Js 
				   && int(fabs(Js - marker[n].Jss)) == Comp)
			N++;
	if (N==0)
	{
		vss = Jss = 0;
		I = 0;
		return;
	}
	vss = new int[N];
	Jss = new int[N];
	I = new double[N];
	for (n=N=0; n < AnzahlMarker; n++)
		if ((marker[n].UState == UState || UState == 0) && marker[n].LState == LState 
				   && marker[n].Iso == Iso && (marker[n].vs == vs || vs == -1) && marker[n].Js == Js 
				   && int(fabs(Js - marker[n].Jss)) == Comp)
	{
		vss[N] = marker[n].vss;
		Jss[N] = marker[n].Jss;
		I[N++] = (marker + n != LaserLinie ? marker[n].HFLM : -1.0);
	}
	for (n=0; n<N; n++) if (I[n] == -1.0)
	{
		if (Comp == 1) for (m=0; (m<N ? vss[m] != vss[n] || m==n : false); m++) ;
		else m=N;
		I[n] = (n < m ? I[n] : (n>0 ? (n<N-1 ? 0.5 * (I[n-1] + I[n+1]) : I[n-1]) 
			     : (n<N-1 ? I[n+1] : 1.0)));
		break;
	}
}

bool Spektrum::addByProgression(ElState *LState, int Iso, int Js, int Comp, Spektrum *RSpektrum, 
								double Ri, double Ra, double Res, int RIso, int RJs, int RComp, 
								double **Result)
{
	int n, m, c, N, NA, *vss, *Jss, *Rvss, *RJss, MJ = 0, Mv = 0, Jd = Js - RJs;
	double *Intensity, ****TermData = 0, *O, *F, *RI;
	bool S;
	if (LState == 0)
	{
		QMessageBox::information(this, "MolSpektAnalysis", 
			"For this function a valid lower state has to be given!", QMessageBox::Ok);
		return false;
	}
	if (RSpektrum != 0) 
		RSpektrum->getProgProfile(0, LState, RIso, -1, RJs, RComp, Rvss, RJss, RI, NA);
	else NA = 0;
	//printf("Nach Profile: NA=%d\n", NA);
	if (NA == 0)
	{
		QMessageBox::information(this, "MolSpektAnalysis", 
		"For the given parameters no progression can be found inside the given reference spectrum!",
	  				QMessageBox::Ok);
		return false;
	}
	if (Comp == RComp) 
	{
		Intensity = new double[N = NA];
		vss = new int[N];
		Jss = new int[N];
		for (n=0; n<N; n++)
		{
			Intensity[n] = RI[n];
			vss[n] = Rvss[n];
			Jss[n] = RJss[n] + Jd;
		}
	}
	else
	{
		if (Comp == 0) 
		{
			vss = new int[NA];
			Jss = new int[NA];
			Intensity = new double[NA];
			for (n=0; n < NA; n++) Jss[n] = 0;
			for (n=N=0; n < NA; n++)
			{
				for (m=0; (m<N ? vss[m] != Rvss[n] : false); m++) ;
				if (m==N) vss[N++] = Rvss[n];
				Jss[m]++;
				Intensity[m] += RI[n];
			}
			for (n=0; n<N; n++)
			{
				Intensity[n] /= Jss[n];
				Jss[n] = Js;
			}
		}
		else
		{
			vss = new int[N = 2 * NA];
			Jss = new int[N];
			Intensity = new double[N];
			for (m=n=0; n < NA; n++) for (c=0; c<2; c++)
			{
				vss[m] = Rvss[n];
				Jss[m] = RJss[n];
				Intensity[m++] = RI[n];
			}
		}
	}
	TermTable *TT = LState->getTermTable();
	//printf("NA=%d, N=%d\n", NA, N);
	if (TT != 0)
	{
		TermData = TT->getData();
		MJ = TT->getMaxJ();
		Mv = TT->getMaxv();
		for (n = NA = 0; n<N; n++) if (vss[n] <= Mv && Jss[n] <= MJ) NA++;
	}
	//printf("NA=%d, Mv=%d, MJ=%d\n", NA, Mv, MJ);
	if (NA > 0 && Iso < TT->getNumIso())
	{
		O = new double[NA];
		F = new double[NA];
		for (n = NA = 0; n<N; n++) if (vss[n] <= Mv && Jss[n] <= MJ) 
				if ((vss[n] > 0 ? 
					TermData[0][Iso][vss[n]][Jss[n]] > TermData[0][Iso][vss[n] - 1][Jss[n]] 
					: true) && (Jss[n] > 0 ? 
					TermData[0][Iso][vss[n]][Jss[n]] > TermData[0][Iso][vss[n]][Jss[n] - 1] : true))
		{
			if ((O[NA] = TermData[0][Iso][vss[n]][Jss[n]]) == 0.0) continue;
			//printf("O[%d]=%f\n", NA, O[NA]);
			F[NA++] = Intensity[n];
		}
		//printf("Vor add, NA=%d\n", NA);
		S = add(O, F, NA, Ri, Ra, Res, Result);
		delete[] O;
		delete[] F;
	}
	else S = false;
	delete[] vss;
	delete[] Jss;
	delete[] RJss;
	delete[] Rvss;
	delete[] RI;
	delete[] Intensity;
	return S;
}

bool Spektrum::addByProgFCF(ElState *UState, ElState *LState, int Iso, int vs, int Js, int Comp, 
                            double T, double SR, double Res, double **Result, int NumWFPoints, double uE, double **uWF, double *TStr,
                            int NWFC, double *const/* o_intensF*/, const bool i_simulate)
{
	if (UState == 0 || LState == 0)
	{
		QMessageBox::information(this, "MolSpektanalysis", 
			"For this function valid lower and upper states have to be given!", QMessageBox::Ok);
		return false;
	}
	Potential *UPot = UState->getPotential(), *LPot = LState->getPotential();
	if (UPot == 0 || LPot == 0)
	{
		QMessageBox::information(this, "MolSpektanalysis", 
			"To use this function potentials for both electronic states have to be available!",
    		QMessageBox::Ok);
		return false;
	}
	TermTable *TT = LState->getTermTable(), *TTU = UState->getTermTable();
	if (TT == 0 || (TTU == 0 && (uWF == 0 || uE == 0.0)))
	{
		QMessageBox::information(this, "MolSpektanalysis", 
			"To use this function term energies have to be available for both states!",
    		QMessageBox::Ok);
		return false;
	}
	int Mv = TT->getMaxv(), MJ = TT->getMaxJ(), NIso = TT->getNumIso(), uc = 0;
	double ****TData = TT->getData();
	if (Iso >= NIso)
	{
		printf("Iso=%d, NIso=%d", Iso, NIso);
		QMessageBox::information(this, "MolSpektanalysis", 
			"Error: for the given isotopomer no ground state term energies are available!",
    		QMessageBox::Ok);
		return false;
	}
	if (uE == 0.0)
	{
		int UMv = TTU->getMaxv(), UMJ = TTU->getMaxJ(), UNIso = TTU->getNumIso(), NC = TTU->getNumComp();
		double ****UTData = TTU->getData();
		if (Iso >= UNIso)
		{
			printf("Iso=%d, UNIso=%d\n", Iso, UNIso);
			QMessageBox::information(this, "MolSpektanalysis", 
				"Error: for the given isotopomer no term energies are available for the excited state!",
				QMessageBox::Ok);
			return false;
		}
		if (vs > UMv || Js > UMJ)
		{
			QMessageBox::information(this, "MolSpektanalysis", 
				"Error: for the given excited states level no term energy is available!",
				QMessageBox::Ok);
			return false;
		}
		if (Js > (Comp == 1 ? MJ : MJ + 1))
		{
			QMessageBox::information(this, "MolSpektanalysis", 
				"Error: for the given J no ground state data is available!",
				QMessageBox::Ok);
			return false;
		}
		uc = (1 < NC ? 1 - Comp : 0);
		uE = UTData[uc][Iso][vs][Js];
	}
	int NL, Nv, NJ, v, n, Jss[2*(Mv+1)];
	double **TPop, *O, *F, E[2*(Mv+1)], FCF[2*(Mv+1)];
	double Ri = uE - SR, Ra = uE + SR;
	//printf("Mv=%d\n", Mv);
	if (Comp == 0) for (n=0; (n <= Mv ? TData[0][Iso][n][Js] != 0.0 : false); n++)
	{
		Jss[n] = Js;
		E[n] = TData[0][Iso][n][Js];
	}
	else for (n=v=0; (v <= Mv ? TData[0][Iso][v][(Js > 0 ? Js - 1 : Js + 1)] != 0.0 : false); v++)
	{
		Jss[n] = Js - 1;
		if (Jss[n] >= 0) E[n++] = TData[0][Iso][v][Js-1];
		Jss[n] = Js + 1;
		if (Jss[n] <= MJ ? TData[0][Iso][v][Jss[n]] != 0.0 : false) E[n++] = TData[0][Iso][v][Js+1];
	}
	//printf("Vor getFastFCF, Comp=%d\n", Comp);
    LPot->getFastFCF(Iso, UPot, Js, uc, uE, NL = n, E, Jss, FCF, NumWFPoints, uWF, TStr, NWFC);
	TT->getThermPopulation(TPop, Nv, NJ, T, Iso, 0);
	O = new double[NL];
	F = new double[NL];
    //if (o_intensF != 0) *o_intensF = 0;
	if (Comp == 0) for (v=0; v < n; v++) 
	{
		O[v] = TData[0][Iso][v][Js];
		F[v] = TPop[v][Js] * FCF[v];
        //if (o_intensF != 0) *o_intensF += F[v] * F[v];
	}
	else for (n=v=0; n < NL; v++)
	{
		//printf("v=%d, Jss[%d]=%d\n", v, n, Jss[n]);
		if (Jss[n] == Js - 1)
		{
			//printf("TPop[%d]=%g, FCF[%d]=%g\n", n, TPop[v][Jss[n]], n, FCF[n]);
			F[n] = TPop[v][Jss[n]] * FCF[n];
			O[n++] = TData[0][Iso][v][Js - 1];
            //if (o_intensF != 0) *o_intensF += F[v] * F[v];
		}
		if (n == NL) break;
		if (Jss[n] == Js + 1)
		{
			F[n] = TPop[v][Jss[n]] * FCF[n];
			O[n++] = TData[0][Iso][v][Js + 1];
            //if (o_intensF != 0) *o_intensF += F[v] * F[v];
		}
	}
    //if (o_intensF != 0) *o_intensF = sqrt(1.0 / *o_intensF);
	//printf("Vor add, v=%d, Nv=%d, NL=%d\n", v, Nv, NL);
    bool S = (i_simulate ? false : add(O, F, NL, Ri, Ra, Res, Result));
	//printf("Vor delete O\n");
	delete[] O;
	//printf("Vor delete F\n");
	delete[] F;
	//printf("Vor Destroy TPop\n");
	Destroy(TPop, Nv);
	//printf("Vor Ende\n");
	return S;
}

void Spektrum::autoAddByProgFCF(ElState *UState, ElState *LState, int Iso, int Comp, double T, 
                                double SR, double Res, int NumWFPoints)
{
	printf("Spektrum::autoAddByProgFCF, Errornous function which may give wrong results, BR gets used uninitialized\n");
	if (QMessageBox::information(this, "MolSpektanalysis", 
			"Spektrum::autoAddByProgFCF: This function contains errors and may not give usefull results!", 
			QMessageBox::Ok | QMessageBox::Cancel) == QMessageBox::Cancel) return;
	if (UState == 0 || LState == 0)
	{
		QMessageBox::information(this, "MolSpektanalysis", 
			"For this function valid lower and upper states have to be given!", QMessageBox::Ok);
		return;
	}
	Potential *UPot = UState->getPotential(), *LPot = LState->getPotential();
	if (UPot == 0 || LPot == 0)
	{
		QMessageBox::information(this, "MolSpektanalysis", 
			"To use this function potentials for both electronic states have to be available!",
    		QMessageBox::Ok);
		return;
	}
	TermTable *TT = LState->getTermTable(), *TTU = UState->getTermTable();
	if (TT == 0 || TTU == 0)
	{
		QMessageBox::information(this, "MolSpektanalysis", 
			"To use this function term energies have to be available for both states!",
    		QMessageBox::Ok);
		return;
	}
	int Mv = TT->getMaxv(), MJ = TT->getMaxJ(), NIso = TT->getNumIso();
	int UMv = TTU->getMaxv(), UMJ = TTU->getMaxJ(), UNIso = TTU->getNumIso(), NC = TTU->getNumComp();
	double ****TData = TT->getData(), ****UTData = TTU->getData();
	if (Iso >= NIso)
	{
		QMessageBox::information(this, "MolSpektanalysis", 
			"Error: for the given isotopomer no term energies are available!",
    		QMessageBox::Ok);
		return;
	}
	if (Iso >= UNIso)
	{
		QMessageBox::information(this, "MolSpektanalysis", 
			"Error: for the given excited states levels no term energy is available!",
    		QMessageBox::Ok);
		return;
	}
	int NL, NUL, nu, Nv, NJ, v, n, m, uc = (1 < NC ? 1 - Comp : 0), Jss[2*(Mv+1)], **A, vs, Js;
	int NP = int(2 * SR / Res) + 1, JStart = UState->getJStart(Iso, Comp), clc, TJ;
	int JStep = (UState->getMolecule() != 0 ? UState->getMolecule()->getJStep(Iso) : 1), nI, JD;
	double **TPop, E[2*(Mv+1)], FCF[2*(Mv+1)], O[2*(Mv+1)], F[2*(Mv+1)], **UT;
	double Ri, Ra, **Result = Create(NP, 2), ED = 0.0, LED, AM, LM = 0.0, AI, AR, BR=0.0;
	if (UMJ > MJ) UMJ = MJ + Comp;
	double ***LF = Create(UMJ + 1, 10, 2), AvI[UMJ+1];
	double *AP = new double[UMJ + 1], BL;
	NUL = (UMJ + 1) * (UMv + 1);
	A = CreateInt(NUL, 2);
	UT = Create(NUL, 4);
	TT->getThermPopulation(TPop, Nv, NJ, T, Iso, 0);
	for (MJ = JStart; (MJ + JStep < NJ ? TPop[0][MJ] < TPop[0][MJ + JStep] : false); MJ += JStep) ;
	printf("Vor Hauptschleife\n");
	for (nu = vs = 0; vs <= UMv; vs++) 
	{
		for (Js = JStart, clc = 0; (Js <= UMJ ? UTData[uc][Iso][vs][Js] != 0 : false); 
				   Js += JStep)
		{
			//printf("Js=%d\n", Js);
			Ri = UTData[uc][Iso][vs][Js] - SR;
			Ra = UTData[uc][Iso][vs][Js] + SR;
			//printf("Mv=%d, MJ=%d, UMv=%d, UMJ=%d\n", Mv, MJ, UMv, UMJ);
			if (Comp == 0) for (n=0; (n <= Mv ? TData[0][Iso][n][Js] != 0.0 : false); n++)
			{
				Jss[n] = Js;
				E[n] = TData[0][Iso][n][Js];
			}
			else 
				for (n=v=0, TJ = (Js > 0 ? Js - 1 : Js + 1); 
								 (v <= Mv ? TData[0][Iso][v][TJ] != 0.0 : false); v++)
			{
				//printf("n=%d, v=%d\n", n, v);
				Jss[n] = Js - 1;
				if (Jss[n] >= 0) E[n++] = TData[0][Iso][v][Js-1];
				Jss[n] = Js + 1;
				if (Jss[n] <= MJ ? TData[0][Iso][v][Jss[n]] != 0.0 : false) 
					E[n++] = TData[0][Iso][v][Js+1];
			}
			//printf("Vor getFastFCF\n");
            LPot->getFastFCF(Iso, UPot, Js, uc, UTData[uc][Iso][vs][Js], NL = n, E, Jss, FCF, NumWFPoints);
			if (Comp == 0) for (v=0; v < n; v++) 
			{	
				O[v] = TData[0][Iso][v][Js];
				F[v] = TPop[v][Js] * FCF[v];
			}	
			else for (n=v=0; n < NL; v++)
			{
				if (Jss[n] == Js - 1)
				{
					F[n] = TPop[v][Jss[n]] * FCF[n];
					O[n++] = TData[0][Iso][v][Js - 1];
				}
				if (Jss[n] == Js + 1)
				{
					F[n] = TPop[v][Jss[n]] * FCF[n];
					O[n++] = TData[0][Iso][v][Js + 1];
				}
			}
			//printf("Vor add\n");
			add(O, F, NL, Ri, Ra, Res, Result);
			for (n=0; n < 10; n++) LF[Js][n][1] = 0.0;
			if (Type == Absorption) for (AM = Result[0][1], AvI[Js] = 0.0, n=1, nI = 0; n < NP - 1; n++)
			{
				if (Result[n-1][1] <= Result[n][1] && Result[n][1] > Result[n+1][1])
				{
					LM = AM;
					AM = Result[n][1];
				}
				else if (Result[n-1][1] >= Result[n][1] && Result[n][1] < Result[n+1][1])
				{
					ParabInterpol(Result[n-1][0], Result[n-1][1], Result[n][0], Result[n][1],
								  Result[n+1][0], Result[n+1][1], AR, AI);
					AI = (LM > AM ? LM - AI : AM - AI);
					if (AI > LF[Js][9][1])
					{
						for (m=8; (m>=0 ? AI > LF[Js][m][1] : false); m--)
						{
							LF[Js][m+1][0] = LF[Js][m][0];
							LF[Js][m+1][1] = LF[Js][m][1];
						}
						LF[Js][++m][0] = AR;
						LF[Js][m][1] = AI;
					}
					AvI[Js] += AI;
					nI++;
				}
			}
			else for (AM = Result[0][1], AvI[Js] = 0.0, n=1, nI = 0; n < NP - 1; n++)
			{
				if (Result[n-1][1] >= Result[n][1] && Result[n][1] < Result[n+1][1])
				{
					LM = AM;
					AM = Result[n][1];
				}
				else if (Result[n-1][1] <= Result[n][1] && Result[n][1] > Result[n+1][1])
				{
					ParabInterpol(Result[n-1][0], Result[n-1][1], Result[n][0], Result[n][1],
								  Result[n+1][0], Result[n+1][1], AR, AI);
					AI = (LM < AM ? AI - LM : AI - AM);
					if (AI > LF[Js][9][1])
					{
						for (m=8; (m>=0 ? AI > LF[Js][m][1] : false); m--)
						{
							LF[Js][m+1][0] = LF[Js][m][0];
							LF[Js][m+1][1] = LF[Js][m][1];
						}
						LF[Js][++m][0] = AR;
						LF[Js][m][1] = AI;
					}
					AvI[Js] += AI;
					nI++;
				}
			}
			//printf("Vor Ende\n");
			if (nI > 0) AvI[Js] /= nI;
		}
		for (n = BL = 0; n<9; n++) for (m=0, AP[MJ] = n; m<9; m++)
		{
			Js = MJ + 2 * JStep;
			JD = (JStep == 1 ? 2 * Js + 1 : 4 * Js + 6);
			for (AP[MJ + JStep] = m; Js <= UMJ; JD += 2 * JStep, Js += JStep)
			{
				
		
		
				LED = ED;
				ED = BR - UTData[uc][Iso][vs][Js];
				if (fabs(ED - LED) > 0.03)
				{
					if (clc <= 2) nu -= clc;
					clc = 1;
				}
				else clc++;
				A[nu][0] = vs;
				A[nu][1] = Js;
				UT[nu][0] = BR;
				//UT[nu][1] = BI;
				//UT[nu][2] = BI / SI;
				//UT[nu++][3] = BI / AvI;
			}
		}
		if (clc <= 2) nu -= clc;
		printf("vs=%d, UMv=%d, nu=%d\n", vs, UMv, nu);
	}
	MeasuredTermEnergies *MTT = MW->CreateMeasuredTermTable();
	if (MTT != 0)
	{
		//UState->addTermTable(MTT);
		MTT->addData(Iso, Comp, nu, UT, A);
		MTT->setName(getName() + "TermEnergies");
		MTT->setSource("own work, " + getName());
		MTT->show();
	}
	Destroy(TPop, Nv);
	Destroy(Result, NP);
	Destroy(A, NUL);
	Destroy(UT, NUL);
}

bool Spektrum::addByBand(ElState *UState, ElState *LState, int Iso, int vs, int vss, int Comp, 
						 double T, double SR, double Res, double **Result)
{
	printf("Spektrum::addByBand\n");
	if (UState == 0 || LState == 0)
	{
		QMessageBox::information(this, "MolSpektanalysis", 
			"For this function valid lower and upper states have to be given!", QMessageBox::Ok);
		return false;
	}
	TermTable *TU = UState->getTermTable(), *TL = LState->getTermTable();
	if (TU == 0 || TL == 0)
	{
		QMessageBox::information(this, "MolSpektanalysis", 
			"This function needs term energy tables for both electronic states!", QMessageBox::Ok);
		return false;
	}
	double ****TDU = TU->getData(), ****TDL = TL->getData(), **TPop, *O, *F;
	int MvU = TU->getMaxv(), MvL = TL->getMaxv(), MJU = TU->getMaxJ(), MJL = TL->getMaxJ(), Nv, NJ, NL;
	int n, J, Jd, JS = UState->getJStart(Iso, Comp);
	Molecule *Mol = LState->getMolecule();
	if (Mol != 0) Jd = Mol->getJStep(Iso);
	else Jd = 1;
	if (Iso >= TU->getNumIso() || Iso >= TL->getNumIso() || MvU < vs || MvL < vss)
	{
		QMessageBox::information(this, "MolSpektanalysis", 
			"Error: for the given isotopomer not enough term energies are available!",
    		QMessageBox::Ok);
		return false;
	}
	TL->getThermPopulation(TPop, Nv, NJ, T, Iso, 0);
	if (Comp == 0)
	{
		for (NL = 1, J = JS + Jd; (J <= MJU ? TDU[0][Iso][vs][J] > TDU[0][Iso][vs][J - Jd] : false)
				   && (J <= MJL ? TDL[0][Iso][vss][J] > TDL[0][Iso][vss][J - Jd] : false); 
				   J += Jd, NL++) ;
		O = new double[NL];
		F = new double[NL];
		for (n=0, J = JS; n < NL; n++, J += Jd)
		{
			O[n] = TDU[0][Iso][vs][0] - TDU[0][Iso][vs][J] - TDL[0][Iso][vss][0] + TDL[0][Iso][vss][J];
			F[n] = TPop[vss][J];
		}
	}
	else
	{
		for (NL = 1, J = JS + Jd; (J <= MJU ? TDU[0][Iso][vs][J] > TDU[0][Iso][vs][J - Jd] : false)
				   && (J < MJL ? TDL[0][Iso][vss][J+1] > TDL[0][Iso][vss][J - Jd + 1] : false); 
				   J += Jd, NL++) ;
		NL *= 2;
		if (JS == 0) NL--;
		if (J == MJL)
		{
			if (J <= MJU) NL++;
			if (Jd == 1 && J < MJU) NL++;
		}
		O = new double[NL];
		F = new double[NL];
		for (n=0, J = JS; n < NL; J += Jd)
		{
			if (J > 0)
			{
				O[n] = TDU[0][Iso][vs][0] - TDU[0][Iso][vs][J] 
						- TDL[0][Iso][vss][0] + TDL[0][Iso][vss][J-1];
				F[n++] = TPop[vss][J-1];
				printf("J-1=%d, n=%d\n", J-1, n);
			}
			if (J < MJL)
			{
				O[n] = TDU[0][Iso][vs][0] - TDU[0][Iso][vs][J]
						- TDL[0][Iso][vss][0] + TDL[0][Iso][vss][J+1];
				F[n++] = TPop[vss][J+1];
				printf("J+1=%d, n=%d\n", J+1, n);
			}
		}
	}
	printf("JS=%d, MJL=%d\n", JS, MJL);
	bool S = add(O, F, NL, TDU[0][Iso][vs][0] - TDL[0][Iso][vss][0] - SR, 
		         TDU[0][Iso][vs][0] - TDL[0][Iso][vss][0] + SR, Res, Result);
	printf("Vor Destroy, NL=%d, Nv=%d, NJ=%d\n", NL, Nv, NJ);
	Destroy(TPop, Nv);
	printf("Nach Destroy\n");
	delete[] O;
	delete[] F;
	if (S) printf("S=true\n");
	printf("Vor Ende\n");
	return S;
}

bool Spektrum::add(double *O, double *F, int N, double Ri, double Ra, double Res, double **nData)
{
	//printf("Spektrum::add\n");
	if (MW == 0) return false;
	int n, m, r, ND = int((Ra - Ri) / Res) + 1, nd = Daten->GetDSL();
	//for (n=0; n<N; n++) printf("O[%d]=%f, F[%d]=%f\n", n, O[n], n, F[n]);
	//printf("Ri=%f,, Ra=%f\n", Ri, Ra);
	//printf("nData=%d, ND=%d\n", nData, ND);
	//double rm = Daten.GetValue(0, 0), RM = Daten.GetValue(nd -1, 0);
	double R, rs = 0.0, rb = 0.0, vs, vb;
	for (m=0, R = Ri; m < ND; m++, R += Res)
	{
		nData[m][0] = R;
		nData[m][1] = 0.0;
	}
	for (n=r=0; n<N; n++) for (m=0, R = Ri - O[n]; m < ND; m++, R += Res)
	{
		while (r < nd ? (rb = Daten->GetValue(r, 0)) < R : false) r++;
		if (r == nd) r--;
		vb = Daten->GetValue(r, 1);
		while (r >= 0 ? (rs = Daten->GetValue(r, 0)) > R : false) r--;
		if (r < 0) r=0;
		vs = Daten->GetValue(r, 1);
		nData[m][1] += F[n] * (rs < R ? (rb > R ? vs + (R - rs) / (rb - rs) * (vb - vs) : vb) : vs);
	}
	//printf("Ende Spektrum::add\n");
	return true;
}

bool Spektrum::addSimulation(const double i_TUpLevel, const double *const i_OffSet, const double *const i_ScaleFakt, const int i_N,
                             const double i_Ei, const double i_Ea, const double i_Resolution, double **const io_Result) const
{
    if (m_simulationProfile.GetNumPoints() == 0) return false;
    double Re = (-m_simulationProfile.GetEnergy(0) > m_simulationProfile.GetEMax() ? -m_simulationProfile.GetEnergy(0) :
                                                                                     m_simulationProfile.GetEMax());
    double dResolution = 1.0 / i_Resolution;
    for (int n = 0; n < i_N; ++n)
    {
        double Ec = i_TUpLevel - i_OffSet[n], Es = Ec - Re, Ee = Ec + Re;
        if (Ee < i_Ei || Es > i_Ea) continue;
        int i = (Es > i_Ei ? static_cast<int>(ceil(dResolution * (Es - i_Ei))) : 0);
        double E = i_Ei + i_Resolution * static_cast<double>(i);
        int j = (E < m_simulationProfile.GetEnergy(0) ? m_simulationProfile.GetNumPoints() - 1 : 0);
        for ( ; E <= Ee && E <= i_Ea; i++, E += i_Resolution)
        {
            double Ep = E - Ec;
            if (Ep < m_simulationProfile.GetEnergy(0)) while (-m_simulationProfile.GetEnergy(j) < Ep) --j;
            else if (Ep <= m_simulationProfile.GetEMax())
            {
                if (m_simulationProfile.GetEnergy(j) > Ep) j=0;
                while (m_simulationProfile.GetEnergy(j) < Ep) ++j;
            }
            else
            {
                if (-m_simulationProfile.GetEnergy(j) > Ep) j=0;
                while (-m_simulationProfile.GetEnergy(j) > Ep) j--;
            }
            io_Result[i][1] -= i_ScaleFakt[n] * m_simulationProfile.GetIntensity(j);
        }
    }
    return true;
}

void Spektrum::optimizeWeighting(double** OffSet, double* TE, int NOff, int NTE, double* WF)
{
	int n, m, r=0, md = Daten->GetDSL() - 1;
	double rb, rs, is, ib, RS, M1, M2, M3, M;
	if (Type == Absorption) for (n=0; n  < NTE; n++)
	{
		for (m=0, WF[n] = 0.0; m < NOff; m++)
		{
			for (RS = TE[n] - OffSet[n][m] - 5.0; r > 0 && Daten->GetValue(r, 0) > RS; r--) ;
			while (r < md && Daten->GetValue(r, 0) < RS) r++;
			for (M1 = Daten->GetValue(r, 1), RS += 4.0; Daten->GetValue(r, 0) < RS && r < md; r++) 
				if ((M = Daten->GetValue(r, 1)) > M1) M1 = M;
			for (M2 = Daten->GetValue(r, 1), RS += 2.0; Daten->GetValue(r, 0) < RS && r < md; r++)
				if ((M = Daten->GetValue(r, 1)) > M2) M2 = M;
			for (M3 = Daten->GetValue(r, 1), RS += 4.0; Daten->GetValue(r, 0) < RS && r < md; r++)
				if ((M = Daten->GetValue(r, 1)) > M3) M3 = M;
			for (RS -= 5.0; r > 0 && Daten->GetValue(r, 0) > RS; r--) ;
			rb = Daten->GetValue(r, 0);
			ib = Daten->GetValue(r, 1);
			if (r < md) r++;
			rs = Daten->GetValue(r, 0);
			is = Daten->GetValue(r, 1);
			WF[n] += ((M3 <= M2 && M2 <= M1) || (M3 >= M2 && M2 >= M1) ? M2 
					 : ((M1 <= M3 && M3 <= M2) || (M1 >= M3 && M3 >= M2) ? M3 : M1))
				   - (rs < RS ? (rb > RS ? is + (RS - rs) / (rb - rs) * (ib - is) : ib) : is);
		}
		WF[n] /= NOff;
	}
	else for (n=0; n < NTE; n++)
	{
		for (m=0, WF[n] = 0.0; m < NOff; m++)
		{
			for (RS = TE[n] - OffSet[n][m] - 5.0; r > 0 && Daten->GetValue(r, 0) > RS; r--) ;
			while (r < md && Daten->GetValue(r, 0) < RS) r++;
			for (M1 = Daten->GetValue(r, 1), RS += 4.0; Daten->GetValue(r, 0) < RS && r < md; r++)
				if ((M = Daten->GetValue(r, 1)) < M1) M1 = M;
			for (M2 = Daten->GetValue(r, 1), RS += 2.0; Daten->GetValue(r, 0) < RS && r < md; r++)
				if ((M = Daten->GetValue(r, 1)) < M2) M2 = M;
			for (M3 = Daten->GetValue(r, 1), RS += 4.0; Daten->GetValue(r, 0) < RS && r < md; r++)
				if ((M = Daten->GetValue(r, 1)) < M3) M3 = M;
			for (RS -= 5.0; r > 0 && Daten->GetValue(r, 0) > RS; r--) ;
			rb = Daten->GetValue(r, 0);
			ib = Daten->GetValue(r, 1);
			if (r < md) r++;
			rs = Daten->GetValue(r, 0);
			is = Daten->GetValue(r, 1);
			WF[n] += (rs < RS ? (rb > RS ? is + (RS - rs) / (rb - rs) * (ib - is) : ib) : is)
				   - ((M3 <= M2 && M2 <= M1) || (M3 >= M2 && M2 >= M1) ? M2
				     : ((M1 <= M3 && M3 <= M2) || (M1 >= M3 && M3 >= M2) ? M3 : M1));
		}
		WF[n] /= NOff;
	}
}

void Spektrum::cut(double *start, double *end, int nR)
{
	if (MW == 0) return;
	int n, N=0, r=0, ND = Daten->GetDSL();
	for(n=0; n < ND; n++)
	{
		if (r < nR ? Daten->GetValue(n, 0) <= start[r] : true) N++;
		else
		{
			while (n < ND ? Daten->GetValue(n, 0) < end[r] : false) n++;
			if (n == ND) break;
			n--;
			r++;
		}
	}
	double **nData = Create(N, 2);
	for (N=n=r=0; n < ND; n++)
	{
		if (r < nR ? Daten->GetValue(n, 0) <= start[r] : true)
		{
			nData[N][0] = Daten->GetValue(n, 0);
			nData[N++][1] = Daten->GetValue(n, 1);
		}
		else
		{
			while (n < ND ? Daten->GetValue(n, 0) < end[r] : false) n++;
			if (n == ND) break;
			n--;
			r++;
		}
	}
	Spektrum *nSpektrum = MW->CreateSpectrum();
	if (nSpektrum != 0)
	{
		nSpektrum->setData(nData, N);
		QString FileName = getFileName();
		if ((n = FileName.indexOf(".")) == -1) n = FileName.length();
		nSpektrum->setFileName(FileName.left(n) + "_cut" + FileName.right(FileName.length() - n));
		nSpektrum->setName(getName() + "_cut");
		nSpektrum->setType(getType());
		nSpektrum->show();
	}
	Destroy(nData, N);
}

void Spektrum::cutAssignedLines()
{
	int n, m=0, N = Daten->GetDSL(), NR = 0, p;
	for (n=0; n < AnzahlMarker; n++) if (marker[n].DisplayData) NR++;
	double Start[NR], Stop[N], v1, v2;
	for (n = NR = 0; m < AnzahlMarker; m++) if (marker[m].DisplayData)
	{
		while (Daten->GetValue(n, 0) < marker[m].Line[0]) n++;
		for (p=n, v1 = Daten->GetValue((--n)--, 1); (n>=0 ? (v2 = Daten->GetValue(n, 1)) < v1 : false); 
				   n--, v1 = v2) ;
		Start[NR] = Daten->GetValue(n+1, 0);
		for (n=p+1, v1 = Daten->GetValue(p, 1); (n<N ? (v2 = Daten->GetValue(n, 1)) < v1 : false); 
				   n++, v1 = v2) ;
		Stop[NR++] = Daten->GetValue(n-1, 0);
		//printf("n=%d, Start[%d]=%f, Stop[%d]=%f\n", n, NR - 1, Start[NR - 1], NR - 1, Stop[NR - 1]);
	}
	cut(Start, Stop, NR);
}

void Spektrum::cutStrongLines()
{
	QDialog *D = new QDialog(this);
	QGridLayout *L = new QGridLayout(D);
	QRadioButton *Min = new QRadioButton("Minima", D), *Max = new QRadioButton("Maxima", D);
	QLineEdit *IntE = new QLineEdit("1", D);
	QPushButton *OK = new QPushButton("OK", D), *Cancel = new QPushButton("Cancel", D);
	D->setWindowTitle("Cut...");
	Min->setChecked(true);
	L->addWidget(Min, 0, 0);
	L->addWidget(Max, 0, 1);
	L->addWidget(new QLabel("stronger than:", D), 1, 0);
	L->addWidget(IntE, 1, 1);
	L->setRowMinimumHeight(2, 20);
	L->addWidget(OK, 3, 0);
	L->addWidget(Cancel, 3, 1);
	connect(OK, SIGNAL(clicked()), D, SLOT(accept()));
	connect(Cancel, SIGNAL(clicked()), D, SLOT(reject()));
	if (D->exec() == QDialog::Accepted)
	{
		int n, N, Mi, Si, DL = Daten->GetDSL();
		double *Start, *Stop, M, i1, i2 , i3, MD = IntE->text().toDouble();
		bool S = false;
		if (Min->isChecked())
		{
			for (n=2, N=0, i2 = Daten->GetValue(0, 1), i3 = Daten->GetValue(1, 1), M = 1e99; n < DL; n++)
			{
				i1 = i2;
				i2 = i3;
				i3 = Daten->GetValue(n, 1);
				if (i1 < i2 && i2 >= i3) M = i2;
				else if (i1 > i2 && i2 <= i3 && M - i2 > MD && M != 1e99) N++;
			}
			Start = new double[N];
			Stop = new double[N];
			for (n=2, N=0, i2 = Daten->GetValue(0, 1), i3 = Daten->GetValue(1, 1), M = 1e99; n < DL; n++)
			{
				i1 = i2;
				i2 = i3;
				i3 = Daten->GetValue(n, 1);
				if (i1 < i2 && i2 >= i3)
				{
					if (S) Si = n-1;
					else 
					{
						Mi = n-1;
						M = i2;
					}
				}
				else if (i1 > i2 && i2 <= i3)
				{
					if (M - i2 > MD && M != 1e99) S = true;
					else if (S)
					{
						S = false;
						Start[N] = Daten->GetValue(Mi, 0);
						Stop[N++] = Daten->GetValue(Si, 0);
					}
				}
			}
		}
		else
		{
			for (n=2, N=0, i2 = Daten->GetValue(0, 1), i3 = Daten->GetValue(1, 1), M = 1e99; n < DL; n++)
			{
				i1 = i2;
				i2 = i3;
				i3 = Daten->GetValue(n, 1);
				if (i1 > i2 && i2 <= i3) M = i2;
				else if (i1 < i2 && i2 >= i3 && i2 - M > MD && M != 1e99) N++;
			}
			Start = new double[N];
			Stop = new double[N];
			for (n=2, N=0, i2 = Daten->GetValue(0, 1), i3 = Daten->GetValue(1, 1), M = 1e99; n < DL; n++)
			{
				i1 = i2;
				i2 = i3;
				i3 = Daten->GetValue(n, 1);
				if (i1 > i2 && i2 <= i3)
				{
					if (S) Si = n-1;
					else 
					{
						Mi = n-1;
						M = i2;
					}
				}
				else if (i1 < i2 && i2 >= i3)
				{
					if (i2 - M > MD && M != 1e99) S = true;
					else if (S)
					{
						S = false;
						Start[N] = Daten->GetValue(Mi, 0);
						Stop[N++] = Daten->GetValue(Si, 0);
					}
				}
			}
		}
		cut (Start, Stop, N);
		delete[] Start;
		delete[] Stop;
	}
	delete D;
}

double Spektrum::getMaxR()
{
	int N = Daten->GetDSL();
	if (N == 0) return 0.0;
	return Daten->GetValue(N-1, 0);
}

void Spektrum::FindProgressions()
{
	if (GMarker == 0) return;
	int n, m, k, l, N, i;
	double MinHeight = 1e2 * Rauschen, adiff, diff, mpos;
	for (n=N=0; GMarker[n] != &MStop; n++) if (GMarker[n]->HFLM > MinHeight) N++;
	Marker *M[N+1], **A;
	int p[N], j, b, bi = 0, pm;
	for (n=0, M[0]=&MBeginn, N=1; GMarker[n] != &MStop; n++) if (GMarker[n]->HFLM > MinHeight) 
			M[N++] = GMarker[n];
	for (n=N-1; n > 1; n--) for (m=n-1; m>0; m--)
	{
		A = new Marker*[N];
		printf("N=%d, n=%d, m=%d\n", N, n, m);
		diff = M[l=n]->Line[0] - M[k=m]->Line[0];
		mpos = M[k]->Line[0] - diff;
		A[0] = M[n];
		A[0]->DD = 0.0;
		for (j=1; j<N; j++) p[j]=1;
		for (j=i=1, pm = -1, b=0; pm != 0; i++)
		{
			A[i]=M[k];
			A[i]->DD = mpos - A[i]->Line[0];
			printf("l=%d, A[%d]=%f, A[%d]=%f, mpos=%f\n", l, i-1, A[i-1]->Line[0], i, A[i]->Line[0], mpos);
			adiff = diff;
			diff = A[i-1]->Line[0] - A[i]->Line[0];
			mpos = M[l=k]->Line[0] - diff * diff / adiff;
			//printf("mpos=%f, M[%d]=%f\n", mpos, l-1, M[l-1]->Line[0]);
			for (k=l-1; M[k]->Line[0] > mpos; k--) ;//printf("M[%d]=%f\n", k, M[k]->Line[0]);
			if (pm == -1) 
			{
				pm = l-k-1;
				//kj = k;
			}
			k += p[i];
			printf("k=%d, l=%d, i=%d\n", k, l, i);
			if (k==l)
			{
				if (i > bi) 
				{
					bi = i;
					b = p[j];
				}
				if (p[j] < pm)
				{
					k = ++p[j];
					i=j+1;
				}
				else
				{
					p[j] = b;
					i=++j;
					pm = -1;
				}
			}
		}
		printf("i=%d\n", i);
		for (k=0; k<i; k++)
		{
			A[k]->LState = A[k]->UState = 0;
			A[k]->Iso = 0;
			A[k]->oc = 0.0;
			A[k]->vs = A[k]->Js = A[k]->Jss = -1;
			A[k]->vss = k;
		}
		Prog.Insert(A, i, i);
	}
	ShowFound();
}

void Spektrum::FindLinesFromTable(LineTable *LTab)
{
	
	int N, n, m;
	double SNR = 0.0;
	TableLine *Lines;
	if (LTab == 0) return;
	if (AnzahlMarker == 0) editFind();
	if (AnzahlMarker == 0) return;
	LTab->getSortedLines(Lines, N, 0);
	if (N==0) return;
	ClearMarked();
	Molecule *Mol = LTab->getMolecule();
	IsoTab *Iso = Mol->getIso();
	ElState *LS = LTab->getTransition()->getLowerState(), *US = LTab->getTransition()->getUpperState();
	QString LSN = LS->getName(), USN = US->getName();
	for (n=m=0; n<N && m < AnzahlMarker; n++)
	{
		while (m < AnzahlMarker ? marker[m].Line[0] < Lines[n].WN : false) m++;
		if (m > 0 ? (m < AnzahlMarker ? Lines[n].WN - marker[m-1].Line[0] < marker[m].Line[0] - Lines[n].WN : true) : false) m--;
		if (fabs(Lines[n].WN - marker[m].Line[0]) < ST && (!marker[m].Marked || Lines[n].SNR > SNR) 
			 && (!marker[m].DisplayData || marker[m].Iso != Lines[n].Iso || marker[m].FC != Lines[n].FC 
				|| marker[m].Js != Lines[n].Js || marker[m].Jss != Lines[n].Jss || marker[m].vs != Lines[n].vs
				|| marker[m].vss != Lines[n].vss || marker[m].UState != US || marker[m].LState != LS))
		{
			SNR = Lines[n].SNR;
			if (marker[m].Marked) marker[m].overlap = true;
			else marker[m].overlap = false;
			marker[m].Marked = true;
			marker[m].DisplayData = true;
			marker[m].FC = Lines[n].FC;
			marker[m].Iso = Lines[n].Iso;
			marker[m].IsoName = Iso->getIsoName(Lines[n].Iso);
			marker[m].Js = Lines[n].Js;
			marker[m].Jss = Lines[n].Jss;
			marker[m].LState = LS;
			marker[m].Mol = Mol;
			marker[m].UState = US;
			marker[m].lState = LSN;
			marker[m].satellite = false;
			marker[m].uState = USN;
			marker[m].uncertainty = 0.05;
			marker[m].vs = Lines[n].vs;
			marker[m].vss = Lines[n].vss;
		}
		else if (marker[m].Marked) marker[m].overlap = true;
	}
	delete Iso;
	delete[] Lines;
	addMarker(false);
	if (MW != 0) MW->setActive(this);
    Paint();
}

void Spektrum::normalize(Spektrum *RefSpekt)
{
	//printf("Spektrum::normalize\n");
	int n, m, N = Daten->GetDSL(), M;
	double V, S=0.0, R = Daten->GetValue(0, 0), **Data, **nData = Create(N, 2);
	RefSpekt->getData(Data, M);
	for (m=0; (m<M ? Data[m][0] < R : false); m++) ;
	for (n=0; n<N; n++)
	{
		nData[n][0] = R = Daten->GetValue(n, 0);
		if (R > Data[m][0] && m<M-1)
		{
			m++;
			S = (Data[m][1] - Data[m-1][1]) / (Data[m][0] - Data[m-1][0]);
		}
		else if (R >= Data[m][0]) S=0.0;
		V = Data[m][1] + (R - Data[m][0]) * S;
		nData[n][1] = (V - Daten->GetValue(n, 1)) / V;
	}
	Destroy(Data, M);
	if (MW == 0) setData(nData, N);
	else
	{
		Spektrum *nSpektrum = MW->CreateSpectrum();
		if (nSpektrum != 0)
		{
			QString Name = getName();
			nSpektrum->setData(nData, N);
			nSpektrum->setName(Name += "_norm");
			nSpektrum->setFileName(Name + ".spect");
			nSpektrum->setType(NormalizedAbsorption);
			nSpektrum->show();
		}
	}
	Destroy(nData, N);
	//printf("Normalize Ende\n");
}

void Spektrum::ShowMarker()
{
	if (AnzahlMarker == 0) editFind();
	else if (!DMarkers) DisplayMarkers();
}

int Spektrum::getNPoints()
{
	return Daten->GetDSL();
}
