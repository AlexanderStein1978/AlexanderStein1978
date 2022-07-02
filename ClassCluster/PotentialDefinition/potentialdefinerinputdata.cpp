#include "potentialdefinerinputdata.h"

PotentialDefinerInputData::PotentialDefinerInputData() : mNumPoints(0), mParticleIndex(0), mFirstBound(nullptr), mSecondBound(nullptr), mThirdBound(nullptr), mFourthBound(nullptr),
    mSecondOrderBound(nullptr), mUnbound(nullptr), mEnd1(), mEnd2(), mBoundChanges()
{
}

PotentialDefinerInputData::~PotentialDefinerInputData()
{
    Destroy();
}

void PotentialDefinerInputData::Destroy()
{
    if (0 < mNumPoints)
    {
        delete[] mFirstBound;
        delete[] mSecondBound;
        delete[] mThirdBound;
        delete[] mFourthBound;
        delete[] mSecondOrderBound;
        delete[] mUnbound;
    }
}

void PotentialDefinerInputData::Rescale(const Vector &end1, const Vector &end2, const int numPoints)
{
    Destroy();
    mNumPoints = numPoints;
    if (0 < numPoints)
    {
        mFirstBound = new double[numPoints];
        mSecondBound = new double[numPoints];
        mThirdBound = new double[numPoints];
        mFourthBound = new double[numPoints];
        mSecondOrderBound = new double[numPoints];
        mUnbound = new double[numPoints];
    }
    mEnd1 = end1;
    mEnd2 = end2;
}
