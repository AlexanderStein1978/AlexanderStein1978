//
// C++ Interface: Js
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2011 - 2019
//
// Copyright: See README file that comes with this source code
//
//


#ifndef JS_H
#define JS_H


struct Js
{
    Js() : J(0), E_calc(0.0), E_obs(0.0), unc(0.0), isUp(true)
    {
    }

    Js(const int i_J, const double i_E_calc, const double i_E_obs, const double i_unc) : J(i_J), E_calc(i_E_calc), E_obs(i_E_obs), unc(i_unc), isUp(true)
    {
    }

    int J;
    double E_calc, E_obs, unc;
    bool isUp;
};

#endif
