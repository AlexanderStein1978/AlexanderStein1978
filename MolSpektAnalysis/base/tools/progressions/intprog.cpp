//
// C++ Implementation: IntProg
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2006 - 2020
//
// Copyright: See README file that comes with this source code
//
//


#include "intprog.h"
#include "elstate.h"
#include "molecule.h"
#include "marker.h"


IntProg::IntProg()
{
	clear();
	marker = new Marker*[1000];
}

void IntProg::copy(const IntProg &C)
{
	//printf("Intprog::CConstructor\n");
	int n;
	for (n=0; n < C.N; n++)
	{
		marker[n] = C.marker[n];
		DD[n] = C.DD[n];
		Jss[n] = C.Jss[n];
		oc[n] = C.oc[n];
		overlap[n] = C.overlap[n];
		SNR[n] = C.SNR[n];
		vss[n] = C.vss[n];
	}
	FQS = C.FQS;
	G = C.G;
	Iso = C.Iso;
	Js = C.Js;
	N = C.N;
	NGD = C.NGD;
	Np = C.Np;
	vs = C.vs;
	satellite = C.satellite;
	lState = C.lState;
	uState = C.uState;
	Mol = C.Mol;
	QF = C.QF;
	FC = C.FC;
}

bool IntProg::isbetter(IntProg *C)
{
	if (uState != 0 && C->uState == 0) return true;
	if (uState == 0 && C->uState != 0) return false;
	if (satellite && !C->satellite) return true;
	if (!satellite && C->satellite) return false;
	if (NGD > C->NGD) return true;
	if (NGD < C->NGD) return false;
	if (G > C->G) return true;
	if (G < C->G) return false;
	if (N > C->N) return true;
	if (N < C->N) return false;
	if (FQS < C->FQS) return true;
	return false;
}

bool IntProg::isIncluded(Marker *M)
{
	int i;
	for (i=0; i<N; i++) if (M == marker[i]) return true;
	return false;
}

IntProg::~IntProg()
{
	//printf("Vor delete marker\n");
	if (marker != NULL) delete[] marker;
	//printf("Nach delete marker\n");
}

void IntProg::clear()
{
	N = G = NGD = vs = Js = Iso = 0;
	FC = -1;
	satellite = false;
	FQS = 0.0;
	lState = uState = 0;
	Mol = 0;
}
