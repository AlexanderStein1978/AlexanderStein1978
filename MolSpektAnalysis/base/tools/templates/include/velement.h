//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
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
