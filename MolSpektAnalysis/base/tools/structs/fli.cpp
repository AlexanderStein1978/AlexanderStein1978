//
// C++ Implementation: FLI
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2006 - 2019
//
// Copyright: See README file that comes with this source code
//
//


#include "fli.h"


FLI::FLI(int n, int g)
{
    vs = new int[N=n];
    vss = new int[N];
    Jss = new int[N];
    Js = new int[N];
    DD = new double[N];
    Iso = new int[N];
	SNR = new double[N];
	oc = new double[N];
	lState = new ElState*[N];
	uState = new ElState*[N];
	satellite = new bool[N];
	overlap = new bool[N];
    G = g;
}