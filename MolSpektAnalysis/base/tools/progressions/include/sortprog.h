//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
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
