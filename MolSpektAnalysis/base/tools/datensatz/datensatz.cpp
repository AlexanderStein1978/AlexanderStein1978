//
// C++ Implementation: Datensatz
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2006 - 2020
//
// Copyright: See README file that comes with this source code
//
//


#include "datensatz.h"
#include "element.h"
#include "marker.h"

#include <cassert>


Datensatz::Datensatz()
{
    //printf("Datensatz::Datensatz()\n");
    AE = new Element;
    AE->First = AE->Next = NULL;
    AE->minIndex = 0;
    AE->maxIndex = 999;
    AIndex = 0;
    G = 1000;
    return;
}
    
Datensatz::~Datensatz()
{
    //printf("Datensatz::~Datensatz()\n");
    while (AE->Next != NULL) AE = AE->Next;
    //printf("Start delete, G=%d\n", G);
	//int i=0;
	while (AE->First != NULL)
    {
		//printf("delete Element %d\n", i++);
		AE = AE->First;
		delete AE->Next;
    }
	//printf("delete last element\n");
    delete AE;
	//printf("Ende ~Datensatz\n");
    return;
}

void Datensatz::reinit()
{
    while (AE->Next != NULL) AE = AE ->Next;
    while (AE->First != NULL)
    {
	AE = AE->First;
	delete AE->Next;
    }
    AE->Next = NULL;
    G = 1000;
    AIndex = 0;
}

int Datensatz::GetDSL() const
{
    return G - 1000 + AIndex;
}
    
void Datensatz::AddValue(const double &x, const double &y, const bool marked)
{
    if (AIndex == 1000)
    {
		AIndex = 0;
		G+=1000;
		AE->Next = new Element;
		AE->Next->First = AE;
		AE->Next->Next = NULL;
		AE->Next->minIndex = AE->minIndex + 1000;
		AE->Next->maxIndex = AE->maxIndex + 1000;
		AE = AE->Next;
    }
    AE->Daten[AIndex][0] = x;
    AE->Daten[AIndex][1] = y;
    AE->marker[AIndex] = NULL;
    if (marked) AE->marker[AIndex] += 1;
    AIndex++;
}
    
void Datensatz::InsertML(Datensatz &ML)
{
    //printf("Insert erreicht\n");
    int i, EIndex=0, MLIndex=0, NML = ML.GetDSL(), 
    N = GetDSL();
    int GN = NML + N;
    Element *New = NULL, *Next = NULL; 
    while (AE->First != NULL) AE = AE->First;
    AIndex = 999;
    G = 0;
    double Nx = ML.GetValue(0, 0), Ax = AE->Daten[0][0];
    //printf("NML=%d, N=%d\n", NML, N);
    //printf("Nx=%f; Ax=%f\n", Nx, Ax);
    for (i=0; i<GN; i++)
    {
		if (AE == NULL)
		{
	    	printf("Fehler: zu große X-Werte der Marker\n");
	    	break;
		}
		Ax = AE->Daten[EIndex][0];
		if (++AIndex == 1000)
		{
	    	//printf("A1000 i=%d\n", i);
	    	if (Next == NULL)
	    	{
				if (New == NULL)
				{
		    		New = new Element;
		    		New->First = NULL;
				}
				else
				{
		    		New->Next = new Element;
		    		New->Next->First = New;
		    		New = New->Next;
				}
	    	}
	    	else
	    	{
				Next->First = New;
				New->Next = Next;
				New = Next;
				Next = NULL;
	    	}
	    	New->minIndex = G;
	    	New->maxIndex = (G += 1000) - 1;
	    	AIndex = 0;
		}
		if (EIndex < N ? (MLIndex < NML ? Nx < Ax : false) : true)
		{
	    	New->Daten[AIndex][0] = Nx;
	    	New->Daten[AIndex][1] = ML.GetValue(MLIndex, 1);
	    	//printf("Daten[%d][0]=%f; Daten[%d][1]=%f\n", AIndex, New->Daten[AIndex][0],
	    	//     AIndex, New->Daten[AIndex][1]);
	    	New->marker[AIndex] = ML.GetMarker(MLIndex++);
	    	//if (New->marker[AIndex] == NULL) 
			//	printf("Ungemarkerte Linie bei AIndex=%d, MLIndex=%d\n", AIndex, MLIndex - 1);
	    	New->marker[AIndex]->Line = New->Daten[AIndex];
	    	if (MLIndex < NML) Nx = ML.GetValue(MLIndex, 0);
		}
		else
		{
	    	New->Daten[AIndex][0] = Ax;
	    	New->Daten[AIndex][1] = AE->Daten[EIndex][1];
	    	New->marker[AIndex] = AE->marker[EIndex++];
	    	if (EIndex == 1000)
	    	{
				//printf("EIndex1000 i=%d, APos=%d\n", i, AE->minIndex+EIndex);
				Next = AE;
				AE = AE->Next;
				EIndex = 0;
	    	}
		}
    }
    //printf("MLIndex=%d, EIndex=%d, AIndex=%d\n", MLIndex, EIndex, AIndex);
	//printf("Nach Schleife\n");
    AIndex++;
	if (Next != NULL) delete Next;
    AE = New;
    AE->Next = NULL;
    //printf("Am Ende von InserML\n");
}

void Datensatz::ReverseOrder()
{
    //printf("ReverseOrder\n");
    Element *NF = new Element;
    NF->First = NULL;
    NF->minIndex = 0;
    NF->maxIndex = 999;
    int i, N = GetDSL(), NIndex = -1, NG = 1000;
    while (AE->Next != NULL) AE = AE->Next;
    for (i=0; i<N; i++)
    {
	if (--AIndex < 0)
	{
	    NF->Next = AE;
	    AE->minIndex = NG;
	    //printf("NG=%d\n", NG);
	    AE->maxIndex = (NG += 1000) - 1;
	    AE = AE->First;
	    NF->Next->First = NF;
	    AIndex = 999;
	}
	if (++NIndex == 1000)
	{
	    NIndex = 0;
	    NF = NF->Next;
	}
	NF->Daten[NIndex][0] = AE->Daten[AIndex][0];
	NF->Daten[NIndex][1] = AE->Daten[AIndex][1];
	NF->marker[NIndex] = AE->marker[AIndex];
    }
    NF->Next = NULL;
    AIndex = NIndex;
    AE->Next = NULL;
    AE->First = NULL;
    delete AE;
    AE = NF;
}

void Datensatz::SetMarker(int i, Marker *nmarker)
{
    //if (i>333000) printf("Daten.SetMarker, i=%d\n", i);
    while (AE->minIndex > i) AE = AE->First;
    while (AE->maxIndex < i) AE = AE->Next;
    i-=AE->minIndex;
    AE->marker[i] = nmarker;
    nmarker->Line = AE->Daten[i];
    nmarker->Visible = false;
    nmarker->Marked = false;
    nmarker->DisplayData = false;
	nmarker->satellite = false;
	nmarker->overlap = false;
	nmarker->oc = 0.0;
	nmarker->HFLM = -1;
}
    
Marker *Datensatz::GetMarker(int i)
{
	if (i < 0 || i >= G)
	{
		printf("Datensatz::GetMarker error: i=%d, G=%d\n", i, G);
		return 0;
	}
    while (AE->minIndex > i) AE = AE->First;
    while (AE->maxIndex < i) AE = AE->Next;
    i-=AE->minIndex;
    return AE->marker[i];
}

double Datensatz::GetValue(const int &i0, const int &i1)
{
    if (i0 < 0 || i0 >= G || G==0)
	{
		printf("Datensatz::GetValue error: i0=%d, G=%d\n", i0, G);
		return 0.0;
	}
	int i=i0;
    //if (DebugPrint) printf("GetValue, AE->minIndex=%d, i0=%d\n", AE->minIndex, i0);
    while (AE->minIndex > i0) AE = AE->First;
    /*if (i0==285000) 
{
    printf("GetValue: beim Beginn, i0=%d\n", i0);
    AE = AE->Next;
    printf("minIndex=%d, maxIndex=%d\n", AE->minIndex, AE->maxIndex);
}*/
    while (AE->maxIndex < i0) AE = AE->Next;
    i-=AE->minIndex;
    //if (DebugPrint) printf("GetValue: vor return,i=%d, i1=%d\n",i,i1);
    return AE->Daten[i][i1];
}

double *Datensatz::GetPoint(int i)
{
	if (i < 0 || i >= G)
	{
		printf("Datensatz::GetPoint error: i=%d, G=%d\n", i, G);
		return 0;
	}
	while (AE->minIndex > i) AE = AE->First;
	while (AE->maxIndex < i) AE = AE->Next;
	i -= AE->minIndex;
	return AE->Daten[i];
}

void Datensatz::SubtractLine(const double Xstart, const double Xend, const int Ndata, const double * const Ydiff, const bool subtract)
{
    while (AE->Daten[0][0] > Xstart + 1e-6) AE = AE->First;
    while (AE->Daten[999][0] < Xstart - 1e-6) AE = AE->Next;
    if (abs(AE->First->Daten[999][0] - Xstart) < abs(AE->Daten[0][0] - Xstart)) AE = AE->First;
    double min = abs(AE->Daten[0][0] - Xstart), d;
    int i=1;
    while (i <= 999 && (d = abs(AE->Daten[i][0])) < min)
    {
        min = d;
        ++i;
    }
    --i;
    assert(min < 1e-6);
    for (int n=0; n < Ndata; ++n, ++i)
    {
        if (i==1000)
        {
            AE = AE->Next;
            i=0;
        }
        if (subtract) AE->Daten[i][1] -= Ydiff[n];
        else AE->Daten[i][1] += Ydiff[n];
    }
    assert(abs(AE->Daten[i-1][0] - Xend) < 1e-6);
}

bool Datensatz::GetMarked(int i)
{
	if (i < 0 || i >= G)
	{
		printf("Datensatz::GetMarked error: i=%d, G=%d\n", i, G);
		return false;
	}
    while (AE->minIndex > i) AE = AE->First;
    while (AE->maxIndex < i) AE = AE->Next;
    i -= AE->minIndex;
    if (AE->marker[i] == NULL) return false;
    return true;
}
