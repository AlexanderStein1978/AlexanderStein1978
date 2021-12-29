#ifndef POTSTRUCT_H
#define POTSTRUCT_H


class Potential;


struct PotStruct
{
    PotStruct()
        : pot(nullptr), VZoom(1.0), RZoom(1.0)
    {
    }

    Potential *pot;
    double VZoom, RZoom;
};

#endif // POTSTRUCT_H
