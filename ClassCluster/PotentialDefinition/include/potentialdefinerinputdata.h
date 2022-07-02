#ifndef POTENTIALDEFINERINPUTDATA_H
#define POTENTIALDEFINERINPUTDATA_H


#include <vector>

#include "vector.h"


class PotentialDefinerInputData
{
public:
    PotentialDefinerInputData();
    ~PotentialDefinerInputData();

    void Rescale(const Vector& end1, const Vector& end2, const int numPoints);

    inline void addBoundChange(const int index)
    {
        mBoundChanges.push_back(index);
    }

    inline double GetFirstBound(int index) const
    {
        return mFirstBound[index];
    }

    inline double GetSecondBound(int index) const
    {
        return mSecondBound[index];
    }

    inline double GetThirdBound(int index) const
    {
        return mThirdBound[index];
    }

    inline double GetFourthBound(int index) const
    {
        return mFourthBound[index];
    }

    inline double GetSecondOrderBound(int index) const
    {
        return mSecondOrderBound[index];
    }

    inline double GetUnbound(int index) const
    {
        return mUnbound[index];
    }

    inline void SetFirstBound(int index, double energy)
    {
        mFirstBound[index] = energy;
    }

    inline void SetSecondBound(int index, double energy)
    {
        mSecondBound[index] = energy;
    }

    inline void SetThirdBound(int index, double energy)
    {
        mThirdBound[index] = energy;
    }

    inline void SetFourthBound(int index, double energy)
    {
        mFourthBound[index] = energy;
    }

    inline void SetSecondOrderBound(int index, double energy)
    {
        mSecondOrderBound[index] = energy;
    }

    inline void SetUnbound(int index, double energy)
    {
        mUnbound[index] = energy;
    }

    inline const Vector& getEnd1() const
    {
        return mEnd1;
    }

    inline const Vector& getEnd2() const
    {
        return mEnd2;
    }

    inline int getNumnPoints() const
    {
        return mNumPoints;
    }

    inline void setParticleIndex(const int index)
    {
        mParticleIndex = index;
    }

    inline int getParticleIndex() const
    {
        return mParticleIndex;
    }

private:
    void Destroy();

    int mNumPoints, mParticleIndex;
    double *mFirstBound, *mSecondBound, *mThirdBound, *mFourthBound, *mSecondOrderBound, *mUnbound;
    Vector mEnd1, mEnd2;
    std::vector<int> mBoundChanges;
};

#endif // POTENTIALDEFINERINPUTDATA_H
