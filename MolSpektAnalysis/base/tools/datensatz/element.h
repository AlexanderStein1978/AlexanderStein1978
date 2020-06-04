//
// C++ Interface: Element
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2006 - 2019
//
// Copyright: See README file that comes with this source code
//
//


#ifndef ELEMENT_H
#define ELEMENT_H


struct Marker;


struct Element
{
    Element *First;
    Element *Next;
    double Daten[1000][2];
    Marker *marker[1000];
    int minIndex;
    int maxIndex;
};


#endif
