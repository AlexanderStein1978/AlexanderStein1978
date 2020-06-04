//
// C++ Interface: vsOListElement
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2008 - 2019
//
// Copyright: See README file that comes with this source code
//
//


#ifndef VSOLISTELEMENT_H
#define VSOLISTELEMENT_H


struct vsOListElement
{
	int vs, Js, Iso, curMaxOffset;
	vsOListElement *next;
};

#endif
