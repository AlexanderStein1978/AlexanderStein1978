//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//

#ifndef INTPROG_H
#define INTPROG_H


class ElState;
class Molecule;

struct Marker;


struct IntProg
{
	IntProg();
	~IntProg();
	void clear();
	void copy(const IntProg&);
	bool isbetter(IntProg*);
	bool isIncluded(Marker*);
	
	int N, G, NGD, Iso, Np;
	Marker **marker;
	bool satellite, overlap[1000];
	int vss[1000], Jss[1000], vs, Js, FC;
	double FQS, uT, oc[1000], DD[1000], SNR[1000], QF;
	ElState *lState, *uState;
	Molecule *Mol;
};

#endif
