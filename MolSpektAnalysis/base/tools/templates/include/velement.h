//
// C++ Interface: VElement
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2006 - 2019
//
// Copyright: See README file that comes with this source code
//
//

#ifndef VELEMENT
#define VELEMENT


template <class T> struct VElement
{
	VElement *next, *prev;
	int startI, endI;
	T Data[1000];
};


#endif
