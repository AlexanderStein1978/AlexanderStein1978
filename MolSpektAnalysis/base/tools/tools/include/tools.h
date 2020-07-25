//
// C++ Interface: tools
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2006 - 2019
//
// Copyright: See README file that comes with this source code
//
//


#ifndef TOOLS_H
#define TOOLS_H


class QLineEdit;


void ParabInterpol(double x1, double y1, double x2, double y2, double x3, double y3, double &x, double &y);
bool SolvLinEqS(double **EqnS, int n);
bool SolvLinEqSwItImp(double **EqnS, int n, double *mFQS = 0);
int Runden(const double &d);
void Test(double &d1, double &d2, QLineEdit *L1, QLineEdit *L2);

#endif
