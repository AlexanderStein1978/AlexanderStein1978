//
// C++ Interface: FieldWindow
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2011 - 2019
//
// Copyright: See README file that comes with this source code
//
//


#ifndef FIELDWINDOW_H
#define FIELDWINDOW_H


#include <QWidget>

struct FitResult;


class FieldWindow : public QWidget
{
public:
	FieldWindow(int ***Field, int N1, int N2);
	void update(FitResult *Results, int NResults, int *Pos1, int *Pos2, int MaxParFits);
	
protected:
	void paintEvent(QPaintEvent *P);
	
private:
	int ***Field, NxPos, NyPos, Step, NResults, *xPos, *yPos, MaxParFits;
	FitResult *Results;
};

#endif
