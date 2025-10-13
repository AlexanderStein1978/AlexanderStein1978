//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//


#ifndef BASEDATA_H
#define BASEDATA_H


struct BaseData
{
    char isotope = 0;
	ushort v = 0;
	ushort J = 0;
	std::string vs;
	ushort Js = 0;
	std::string source;
	int prog = 0;
	std::string file;
	double energy = 0.0;
	double uncert = 0.0;
	double obs_calc = 0.0;
	float devR = 0.0f;
	std::string secondState;
};

#endif
