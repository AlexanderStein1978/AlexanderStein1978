//
// C++ Interface: SimLine
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2009 - 2019
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
