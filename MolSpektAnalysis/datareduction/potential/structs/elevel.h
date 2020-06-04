//
// C++ Interface: ELevel
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2007 - 2019
//
// Copyright: See README file that comes with this source code
//
//


#ifndef ELEVEL_H
#define ELEVEL_H


struct ELevel
{
	int J, Iso, v;
	double E, *K;
	bool isGT(ELevel L);
	bool isGT(int Iso, int J, int v);
};

#endif
