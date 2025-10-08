//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
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
