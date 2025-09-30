//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//


#ifndef SIMLINE_H
#define SIMLINE_H


class ElState;


struct SimLine
{
	double R, I;
	ElState *Su, *Sl;
	int Iso, vs, vss, Js, Jss;
};

#endif
