//
// C++ Interface: CouplingFuncs
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2007 - 2019
//
// Copyright: See README file that comes with this source code
//
//


#ifndef COUPLINGFUNCS_H
#define COUPLINGFUNCS_H


#include "potential.h"


class CouplingFuncs : public Potential
{
public:

    enum Type{SimpleCouplingFuncs, MLRCouplingFuncs};

    CouplingFuncs(MainWindow* MW);

    inline Type getFunctionType() const
    {
        return m_type;
    }

    inline void SetPotentials(Potential** potentials)
    {
        m_potentials = potentials;
    }

    void setData(int Nb, double* b, double binf, double Rs, double Rc, double epsilon, double xi_Rx_s, double xi_eps);
    void getData(int& Nb, double*& b, double& binf, double& Rs, double &Rc, double& epsilon);
    void getTexTable();

private:

    Type m_type;
    int m_Nb;
    Potential** m_potentials;
};

#endif
