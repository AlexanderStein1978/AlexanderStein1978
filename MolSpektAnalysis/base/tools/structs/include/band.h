//
// C++ Interface: Band
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2006 - 2019
//
// Copyright: See README file that comes with this source code
//
//


#ifndef BAND_H
#define BAND_H


struct Band
{
	Marker *lines[1000];
	int MinI, MaxI, CL1, CL2;
};

#endif
