//
// C++ Interface: SortProg
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2006 - 2019
//
// Copyright: See README file that comes with this source code
//
//


#ifndef SORTPROG_H
#define SORTPROG_H


struct Marker;


struct SortProg
{
	int G, N, Iso;
	double SDiff;
	Marker **marker;
	SortProg *next;
	int *vss, *Jss, *vs, *Js;
	double *Diff;
};


#endif
