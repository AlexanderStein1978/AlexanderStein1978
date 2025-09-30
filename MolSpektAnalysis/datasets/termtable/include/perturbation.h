//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//


#ifndef PERTURBATION_H
#define PERTURBATION_H


#include <QString>

class TermTable;


struct Perturbation
{
	int Comp, Iso, v, J, Pv, PComp;
	TermTable *Perturber;
	QString PName;
};

#endif
