#ifndef POTSTRUCT_H
#define POTSTRUCT_H


#include <QString>


class Potential;


class PotStruct
{
public:
    PotStruct();
    ~PotStruct();
    PotStruct(const PotStruct&) = delete;
    PotStruct operator=(const PotStruct&) = delete;

    Potential& getPotential();
    void set(Potential* const pot, const double vZoom, const double rZoom);
    void InitAsOwner(const QString fileName);

    inline double getVZoom() const
    {
        return mVZoom;
    }

    inline double getRZoom() const
    {
        return mRZoom;
    }

private:
    bool mIOwner;
    Potential *mPot;
    double mVZoom, mRZoom;
};

#endif // POTSTRUCT_H
