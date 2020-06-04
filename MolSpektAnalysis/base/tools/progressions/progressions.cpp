//
// C++ Implementation: Progressions
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2006 - 2019
//
// Copyright: See README file that comes with this source code
//
//


#include "progressions.h"
#include "marker.h"
#include "fli.h"
#include "intprog.h"


Progressions::Progressions() 
{ 
    P = NULL; 
}

Progressions::~Progressions() 
{ 
    Clear(); 
}



void Progressions::Insert(Marker **marker, const int &N, const int &G)
{
    //printf("Progressions::Insert N=%d, G=%d\n", N, G);
	FLI *NFLI = new FLI(N, G);
    NFLI->marker = marker;
	int i;
    for (i=0; i<N; i++)
    {
		NFLI->vss[i] = marker[i]->vss;
		NFLI->Jss[i] = marker[i]->Jss;
		NFLI->DD[i] = marker[i]->DD;
		NFLI->SNR[i] = marker[i]->SNR;
		NFLI->Iso[i] = marker[i]->Iso;
		NFLI->vs[i] = marker[i]->vs;
		NFLI->Js[i] = marker[i]->Js;
		NFLI->oc[i] = marker[i]->oc;
		NFLI->lState[i] = marker[i]->LState;
		NFLI->uState[i] = marker[i]->UState;
		NFLI->satellite[i] = marker[i]->satellite;
		NFLI->overlap[i] = marker[i]->overlap;
    }
    Insert(NFLI);
}

void Progressions::Insert(const IntProg* Prog, const int& G)
{
	FLI *NFLI = new FLI(Prog->N, G);
	int i;
	NFLI->marker = new Marker*[Prog->N];
	for (i=0; i < Prog->N; i++)
	{
		NFLI->marker[i] = Prog->marker[i];
		NFLI->vss[i] = Prog->vss[i];
		NFLI->Jss[i] = Prog->Jss[i];
		NFLI->DD[i] = Prog->DD[i];
		NFLI->SNR[i] = Prog->SNR[i];
		NFLI->Iso[i] = Prog->Iso;
		NFLI->vs[i] = Prog->vs;
		NFLI->Js[i] = Prog->Js;
		NFLI->oc[i] = Prog->oc[i];
		NFLI->lState[i] = Prog->lState;
		NFLI->uState[i] = Prog->uState;
		NFLI->satellite[i] = Prog->satellite;
		NFLI->overlap[i] = Prog->overlap[i];
	}
	Insert(NFLI);
}

void Progressions::Insert(FLI* NFLI)
{
	if (P == NULL) 
    {
		P = NFLI;
		P->First = P->Next = NULL;
    }
    else 
    {
		//printf("Vor While\n");
		while(P->First != NULL) P = P->First;
	   	 //printf("Vor For\n");
		for (; P->G > NFLI->G && P->Next != NULL; P = P->Next) ;
		if (P->G > NFLI->G)
		{
		    //printf("Pos1\n");
	    	P->Next = NFLI;
	    	NFLI->Next = NULL;
	    	NFLI->First = P;
		}
		else
		{
		    //printf("Pos2\n");
		    NFLI->Next = P;
	    	NFLI->First = P->First;
	    	P->First = NFLI;
	    	if (NFLI->First!=NULL) NFLI->First->Next = NFLI;
	    	P = NFLI;
		}
		while(P->First != NULL) P = P->First;
    }
}
    
void Progressions::RotateB()
{
    if (P!=NULL) if (P->First != NULL) P = P->First;
}

void Progressions::RotateF()
{
    if (P!=NULL) if (P->Next != NULL) P = P->Next;
}

void Progressions::GetMarker(Marker **&marker, int &N)
{
    //printf("Progressions::GetMarker\n");
	if (P==NULL) 
    {
		marker = NULL;
		N = 0;
		return;
    }
    N = P->N;
    marker = P->marker;
    int i;
	IsoTab *Iso = 0;
	Molecule *Mol = 0;
    for (i=0; i<N; i++)
    {
		marker[i]->Marked = true;
		marker[i]->DisplayData = true;
		marker[i]->added = false;
		marker[i]->vss = P->vss[i];
		marker[i]->Jss = P->Jss[i];
		marker[i]->DD = P->DD[i];
		marker[i]->SNR = P->SNR[i];
		marker[i]->Iso = P->Iso[i];
		marker[i]->vs = P->vs[i];
		marker[i]->Js = P->Js[i];
		marker[i]->satellite = P->satellite[i];
		marker[i]->oc = P->oc[i];
		marker[i]->lState = ((marker[i]->LState = P->lState[i]) != 0 ? marker[i]->LState->getName() : "");
		marker[i]->uState = ((marker[i]->UState = P->uState[i]) != 0 ? marker[i]->UState->getName() : "");
		marker[i]->Mol = (marker[i]->LState != 0 ? marker[i]->LState->getMolecule() 
			                 : (marker[i]->UState != 0 ? marker[i]->UState->getMolecule() : 0));
		if (marker[i]->Mol != 0 && marker[i]->Mol != Mol) Iso = (Mol = marker[i]->Mol)->getIso();
		marker[i]->IsoName = (Iso != 0 ? Iso->texName[marker[i]->Iso] : "");
		marker[i]->overlap = P->overlap[i];
		marker[i]->FC = -1;
		marker[i]->uncertainty = (marker[i]->overlap ? 0.01 : 0.005);
    }
	//printf("Ende GetMarker\n");
}

void Progressions::Clear()
{
    if (P==NULL) return;
    while (P->First != NULL) P = P->First;
    FLI *F;
    while (P != NULL)
    {
		delete[] P->marker;
		delete[] P->vss;
		delete[] P->Jss;
		delete[] P->DD;
		delete[] P->Iso;
		delete[] P->Js;
		delete[] P->vs;
		delete[] P->SNR;
		delete[] P->lState;
		delete[] P->uState;
		delete[] P->satellite;
		delete[] P->overlap;
		delete[] P->oc;
		F = P;
		P = P->Next;
		delete F;
    }
}
