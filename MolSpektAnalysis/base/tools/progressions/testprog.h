//
// C++ Interface: TestProg
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2006 - 2019
//
// Copyright: See README file that comes with this source code
//
//


#ifndef TESTPROG_H
#define TESTPROG_H


class ElState;


struct TestProg
{
    int *vss, *Jss, *Iso, *vs, *Js;
    double *DD, SBD, SDiff, oc;
	ElState **lState, **uState;
};


#endif
