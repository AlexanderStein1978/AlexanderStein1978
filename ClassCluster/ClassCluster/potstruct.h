#ifndef POTSTRUCT_H
#define POTSTRUCT_H


class Potential;


struct PotStruct
{
    PotStruct(Potential* const Pot)
        : pot(Pot), VZoom(1.0), RZoom(1.0)
    {
    }

    Potential *pot;
    double VZoom, RZoom;
};

#endif // POTSTRUCT_H
