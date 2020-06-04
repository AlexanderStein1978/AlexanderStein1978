//
// C++ Interface: MarkedPeak
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2009 - 2019
//
// Copyright: See README file that comes with this source code
//
//


#ifndef MARKEDPEAK_H
#define MARKEDPEAK_H


struct MarkedPeak
{
	int JPos, EPos, LevelState, Row, StateNum, SplineNum;
	bool CalcFromProgression, CorrectVibLevel;
};

#endif
