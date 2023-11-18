#include "potstruct.h"
#include "potential.h"


PotStruct::PotStruct() : mIOwner(false), mPot(nullptr), mVZoom(1.0), mRZoom(1.0)
{
}


PotStruct::~PotStruct()
{
    if (mIOwner) delete mPot;
}

Potential & PotStruct::getPotential()
{
    return *mPot;
}

void PotStruct::set(Potential *const pot, const double vZoom, const double rZoom)
{
    mPot = pot;
    mVZoom = vZoom;
    mRZoom = rZoom;
}

void PotStruct::InitAsOwner(const QString fileName)
{
    mIOwner = true;
    mPot = new Potential;
    mPot->readData(fileName);
}
