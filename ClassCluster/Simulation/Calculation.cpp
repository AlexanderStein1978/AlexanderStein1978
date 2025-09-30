//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//

#define _USE_MATH_DEFINES

#include "Calculation.h"
#include "heapsort.h"
#include "random_sortfunctor.h"
#include "yz_sortfunctor.h"
#include "vsortfunctor.h"
#include "utils.h"
#include "potstruct.h"
#include "potential.h"
#include "watchpoint.h"
#include "vector.h"
#include "potentialdefinerinputdata.h"

#include <QGridLayout>
#include <QFile>
#include <QLabel>
#include <QTextStream>
#include <QPainter>
#include <QPaintEvent>
#include <QDebug>

#include <math.h>
#include <cstdlib>


const double UaMax = 1e6;


Calculation::Calculation(PotStruct* PotSs, QObject* parent): QThread(parent), Error_Double(0.0/0.0), N(898), NPot(30000), watchParticle(-1), particleWatchStep(-1), mInstanceId(-1), mMaxIt(-1), PS(1e3),
    Pot(new double*[NumPot]), dPdR(new double*[NumPot]), U(0.0), T(0.0), E(0.0), waveStep(10.0), waveAmp(10.0), potRangeScale(PS), mMaxCalcResult(0.0), mCurEDevUB(0.0), mAbsEDevUB(0.0), mLastEnergy(0.0),
    mRandPOF1(static_cast<double>(N) / (static_cast<double>(RAND_MAX) + 1.0)), mRandPOF2(static_cast<double>(N-1) / (static_cast<double>(RAND_MAX) + 1.0)), FixedWallPos(new Particle*[86]),
    writeSnapShot(false), mRotationChanged(false), mEnergyCsvLogFile(nullptr), mEnergyCsvLog(nullptr)
{
	//printf("Calculation::Calculation\n");
    mEnergyCsvLogFilename = "test.csv";
    double IntDist = 20.0, st;
    int n, x, y, z;
    st = 1/PS;
	Re = 4.0;
    nx = ny = nz = 10;
	Rm = st;
	RM = st * double(NPot);
	h = 0.001; //DefaultStepSize 
	Energy = 0.0; //Default energy per particle
	XS = YS = ZS = 20;
	ScF = 10.0; //Scaling factor
    //QString PotFile = "Potential.dat";
	GridSizeDiv = 5;
	rotated = false;
	Move = false;
	Speed = 1e3;
	
	MaxX = IntDist * double(XS / GridSizeDiv);
	MaxY = IntDist * double(YS / GridSizeDiv);
	MaxZ = IntDist * double(ZS / GridSizeDiv);
	
	YMid = 0.5 * MaxY;
	
    for (n=0; n < NumPot; ++n)
    {
        Pot[n] = nullptr;
        dPdR[n] = nullptr;
        potentialOK[n] = false;
    }
    if (nullptr != PotSs) for (n=0; n < NumPot; ++n) setPotential(static_cast<PotRole>(n), PotSs[n]);
	PZS = 23;
	PYS = int(round(MaxY / Re)); 
	PXS = int(round(MaxX / Re));
	P = new Particle[N];
	D = new Particle*[N];
	G = new Particle***[XS];
    MAR = new MARStruct*[N];
	for (x=0; x < XS; x++)
	{
		G[x] = new Particle**[YS];
		for (y=0; y < YS; y++)
        {
            G[x][y] = new Particle*[ZS];
            for (z=0; z < ZS; ++z) G[x][y][z] = nullptr;
        }
	}
    for (n=0; n < N; n++) MAR[n] = new MARStruct[Particle::BoundAL];
    Pos = new Vector[N];

	initialize();

    /*DebugLogFile = new QFile("DebugLog.txt");
    DebugLogFile->open(QIODevice::WriteOnly);
    DebugLog = new QTextStream(DebugLogFile);
    *DebugLog << "It\t";
    for (n=0; n<N; ++n) *DebugLog << "E[" << n << "]    \tT    \tdeltaE\tnew E  \tM\t";
    *DebugLog << "\n";*/
}

Calculation::~Calculation()
{
	int x, y, n;
    for (n=0; n < NumPot; ++n)
    {
        if (nullptr != Pot[n]) delete[] Pot[n];
        if (nullptr != dPdR[n]) delete[] dPdR[n];
    }
    delete[] Pot;
    delete[] dPdR;
	delete[] P;
	delete[] D;
    delete[] Pos;
	for (n=0; n<N; n++) delete[] MAR[n];
	delete[] MAR;
	for (x=0; x < XS; x++)
	{
		for (y=0; y < YS; y++) delete[] G[x][y];
		delete[] G[x];
	}
	delete[] G;
    delete[] FixedWallPos;
}

void Calculation::calcMAR()
{
    int n, p, mx, my, mz, lx, ly, lz, i1, i2;
    double r;
    Particle *PP1, *PP2;
    for (n=0; n<N; n++) for (p=0; p < Particle::BoundAL; p++) MAR[n][p].R = RM;
    for (mz = 0; mz < ZS; mz++) for (my = 0; my < YS; my++) for (mx = 0; mx < XS; mx++)
    {
        if (G[mx][my][mz] != 0)
        {
            for (PP1 = G[mx][my][mz]; PP1->next != 0; PP1 = PP1->next)
                for (PP2 = PP1->next; PP2 != 0; PP2 = PP2->next)
            {
                i1 = PP1 - P;
                i2 = PP2 - P;
                r = (PP2->R - PP1->R).length();
                for (n=0; (n < Particle::BoundAL ? r > MAR[i1][n].R : false); n++) ;
                for (p = Particle::BoundAL - 1; p>n; p--) MAR[i1][p] = MAR[i1][p-1];
                if (n < Particle::BoundAL)
                {
                    MAR[i1][n].R = r;
                    MAR[i1][n].index = i2;
                }
                for (n=0; (n < Particle::BoundAL ? r > MAR[i2][n].R : false); n++) ;
                for (p = Particle::BoundAL - 1; p>n; p--) MAR[i2][p] = MAR[i2][p-1];
                if (n < Particle::BoundAL)
                {
                    MAR[i2][n].R = r;
                    MAR[i2][n].index = i1;
                }
            }
            for (lz = mz; lz < ZS && lz <= mz + GridSizeDiv; lz++)
                for (ly = (lz > mz ? ((ly = my - GridSizeDiv) >= 0 ? ly : 0) : my);
                        ly < YS && ly <= my + GridSizeDiv; ly++)
                    for (lx = (ly > my || lz > mz ? ((lx = mx - GridSizeDiv) >= 0 ? lx : 0) : mx + 1);
                            lx < XS && lx <= mx + GridSizeDiv; lx++)
                        for (PP2 = G[lx][ly][lz]; PP2 != 0; PP2 = PP2->next)
                            for (PP1 = G[mx][my][mz]; PP1 != 0; PP1 = PP1->next)
            {
                i1 = PP1 - P;
                i2 = PP2 - P;
                r = (PP2->R - PP1->R).length();
                for (n=0; (n < Particle::BoundAL ? r > MAR[i1][n].R : false); n++) ;
                for (p = Particle::BoundAL - 1; p>n; p--) MAR[i1][p] = MAR[i1][p-1];
                if (n < Particle::BoundAL)
                {
                    MAR[i1][n].R = r;
                    MAR[i1][n].index = i2;
                }
                for (n=0; (n < Particle::BoundAL ? r > MAR[i2][n].R : false); n++) ;
                for (p = Particle::BoundAL - 1; p>n; p--) MAR[i2][p] = MAR[i2][p-1];
                if (n < Particle::BoundAL)
                {
                    MAR[i2][n].R = r;
                    MAR[i2][n].index = i1;
                }
            }
        }
    }
}

Calculation::Result Calculation::geta(Vector* t0, Vector *a, const bool collectCandidates)
{
	//printf("geta\n");
    int mx, my, mz, lx, ly, lz;// i1, i2, p, n;
    //double r, dx, dy, dz;
	Particle *PP1, *PP2;
    for (mx = 0; mx < N; mx++) a[mx].clear();
	for (mz = 0; mz < ZS; mz++) for (my = 0; my < YS; my++) for (mx = 0; mx < XS; mx++)
	{
		if (G[mx][my][mz] != 0)
		{
			//printf("l0\n");
			for (PP1 = G[mx][my][mz]; PP1->next != 0; PP1 = PP1->next)
                for (PP2 = PP1->next; PP2 != 0; PP2 = PP2->next) if (getU(PP1, PP2, U, t0, temporaryPos, a, collectCandidates) == Error) return Error;
			//printf("l1\n");
			for (lz = mz; lz < ZS && lz <= mz + GridSizeDiv; lz++)
				for (ly = (lz > mz ? ((ly = my - GridSizeDiv) >= 0 ? ly : 0) : my); 
						ly < YS && ly <= my + GridSizeDiv; ly++)
					for (lx = (ly > my || lz > mz ? ((lx = mx - GridSizeDiv) >= 0 ? lx : 0) : mx + 1);
							lx < XS && lx <= mx + GridSizeDiv; lx++)
						for (PP2 = G[lx][ly][lz]; PP2 != 0; PP2 = PP2->next)
                            for (PP1 = G[mx][my][mz]; PP1 != 0; PP1 = PP1->next) if (getU(PP1, PP2, U, t0, temporaryPos, a, collectCandidates) == Error) return Error;
		}
	}
    /*for (n=0; n<N; ++n) if (isnan(ax[n]) || isnan(ay[n]) || isnan(az[n]))
    {
		printf("After calc a: Particel %d is nan!\n", n);
		Run = false;
    }*/
	//printf("U=%f\n", U);
	T *= 0.5;
	//printf("U=%f, T=%f, U+T=%f, E=%f\n", U, T, U+T, E);
	//printf("End geta\n");
    return Success;
}

bool Calculation::wasStepOK() const
{
    int mx, my, mz, lx, ly, lz;
	Particle *PP1, *PP2;
	for (mz = 0; mz < ZS; mz++) for (my = 0; my < YS; my++) for (mx = 0; mx < XS; mx++)
	{
		if (G[mx][my][mz] != 0)
		{
			for (PP1 = G[mx][my][mz]; PP1->next != 0; PP1 = PP1->next)
                for (PP2 = PP1->next; PP2 != 0; PP2 = PP2->next) if ((PP1->R - PP2->R).dot(PP1->lR - PP2->lR) < 0.0) return false;
			for (lz = mz; lz < ZS && lz <= mz + GridSizeDiv; lz++)
				for (ly = (lz > mz ? ((ly = my - GridSizeDiv) >= 0 ? ly : 0) : my); 
						ly < YS && ly <= my + GridSizeDiv; ly++)
					for (lx = (ly > my || lz > mz ? ((lx = mx - GridSizeDiv) >= 0 ? lx : 0) : mx + 1);
							lx < XS && lx <= mx + GridSizeDiv; lx++)
						for (PP2 = G[lx][ly][lz]; PP2 != 0; PP2 = PP2->next)
                            for (PP1 = G[mx][my][mz]; PP1 != 0; PP1 = PP1->next) if ((PP1->R - PP2->R).dot(PP1->lR - PP2->lR) < 0.0) return false;
		}
	}
	return true;
}

Calculation::Result Calculation::getU(Particle *const P1, Particle *const P2, double& U, const Vector *const t0, Positions pos, Vector* a, const bool collectCandidates) const
{
    // int *debugnullptr = nullptr;
    double r, amp(0.0);
    Vector d, b, d1, d2;
    bool calcA = (NULL != a);
    int i1 = P1 - P, i2 = P2 - P, p, bi1=N, bi2=N;
    static Particle *errorP1(nullptr), *errorP2(nullptr);
    if (errorP1 == P1 && errorP2 == P2) printf("i1=%d, i2=%d, P1->x=%g, P1->y=%g, P1->z=%g, P2->x=%g, P2->y=%g, P2->z=%g\n", i1, i2, P1->R.X(), P1->R.Y(), P1->R.Z(), P2->R.X(), P2->R.Y(), P2->R.Z());
    switch (pos)
    {
    case temporaryPos:
        d = t0[i2] - t0[i1];
        break;
    case lastPos:
        d = P2->lR - *t0;
        break;
    case currentPos:
        d = P2->R - *t0;
        break;
    case particles:
        d = P2->R - P1->R;
        break;
    }
    //printf("tx1=%f, tx2=%f, ty1=%f, ty2=%f, tz1=%f, tz2=%f\n",
        //   tx[i1], tx[i2], ty[i1], ty[i2], tz[i1], tz[i2]);
    r = d.length();
    p = int((r - Rm) * potRangeScale);
    const double dR = 1.0 / r;
    if (p < 0)
    {
        //*debugnullptr = 5;
        return Error;
    }
    for (int n=0; n < Particle::BoundAL; ++n)
    {
        if (n < P1->NB && P1->bound[n].p == P2) bi1=n;
        if (n < P2->NB && P2->bound[n].p == P1) bi2=n;
    }
    if (p >= NPot)
    {
        if (bi1 <= P1->NB && bi2 <= P2->NB)
        {
            removeBinding(P1, bi1);
            removeBinding(P2, bi2);
            if (isBindingDoubled(P1 - P)) return Error; // *debugnullptr = 5;
            if (isBindingDoubled(P2 - P)) return Error;  // *debugnullptr = 5;
        }
        else if (bi1 <= P1->NB || bi2 <= P2->NB) return Error; // *debugnullptr = 5;
        return Success;
    }
    //printf("p=%d, r=%f, Rm=%f, PS=%f\n", p, r, Rm, PS);
    if (bi1 < P1->NB && bi2 < P2->NB)
    {
        if (bi1 <= 1 || bi2 <= 1)
        {
            //if (Pot[ClosestTwo][p] > UaMax)
              //  return Error;
            if (calcA) amp = dPdR[ClosestTwo][p] * dR;
            U += Pot[ClosestTwo][p];
            //printf("Amp=%g, dR=%g, dPdR[ClosestTwo][%d]=%g\n", amp, dR, p, dPdR[ClosestTwo][p]);
        }
        else 
        {
            //if (Pot[NextTwo][p] > UaMax)
              //  return Error;
            if (calcA) amp = dPdR[NextTwo][p] * dR;
            U += Pot[NextTwo][p];
            //printf("Amp=%g, dR=%g, dPdR[NextTwo][%d]=%g\n", amp, dR, p, dPdR[NextTwo][p]);
        }
        if (collectCandidates)
        {
            P1->bound[bi1].lastDist = r;
            P2->bound[bi2].lastDist = r;
        }
    }
    else
    {
        if (collectCandidates)
        {
            addCandidate(P1, P2, r);
            addCandidate(P2, P1, r);
        }
        //if (Pot[Remaining][p] > UaMax)
              //  return Error;
        if (calcA) amp = dPdR[Remaining][p] * dR;
        U += Pot[Remaining][p];
        //printf("Amp=%g, dR=%g, dPdR[Remaining][%d]=%g\n", amp, dR, p, dPdR[Remaining][p]);
        for (int n=0; n < P1->NB; ++n) for (int m=0; m < P1->bound[n].p->NB; ++m) if (P1->bound[n].p->bound[m].p == P1)
        {
            int i = getPartnerBindingIndex(P1->bound[n].p, m);
            if (i != -1 && P1->bound[n].p->bound[i].p == P2)
            {
                int im = P1->bound[n].p - P;
                switch (pos)
                {
                case temporaryPos:
                    d1 = t0[i1] - t0[im];
                    d2 = t0[i2] - t0[im];
                    break;
                case lastPos:
                    d1 = *t0 - P[im].lR;
                    d2 = P2->lR - P[im].lR;
                    break;
                case currentPos:
                    d1 = *t0 - P[im].R;
                    d2 = P2->R - P[im].R;
                    break;
                case particles:
                    d1 = P1->R - P1->bound[n].p->R;
                    d2 = P2->R - P1->bound[n].p->R;
                    break;
                }
                double dr1 = 1.0 / d1.length(), dr2 = 1.0 / d2.length();         
                int ap = static_cast<int>((1.0 + d1.dot(d2)*dr1*dr2) * 0.5 * (NPot - 1));
                U += Pot[Angular][ap];
                if (calcA)
                {
                    Vector a1 = d1.cross(d1.cross(d2)).unit() * dr1 * dPdR[Angular][ap];
                    Vector a2 = d2.cross(d2.cross(d1)).unit() * dr2 * dPdR[Angular][ap];
                    a[i1] += a1;
                    a[i2] += a2;
                    a[im] -= a1 + a2;
                }
            }
        }
    }
    /*if (abs(amp * r) > UaMax)
    {
        //printf("i1=%d, i2=%d, a=%f: Stepsize gets reduced!\n", i1, i2, amp*r);
        errorP1 = P1;
        errorP2 = P2;
        return Error;
    }*/
    if (calcA)
    {
        a[i1] += (b = amp * d);
        a[i2] -= b;
        if (isnan(a[i1].X()) || isnan(a[i1].Y()) || isnan(a[i1].Z()) || isnan(a[i2].X()) || isnan(a[i2].Y()) || isnan(a[i2].Z())) return Error; // *debugnullptr = 5;
        if (particleWatchStep >= 0)
        {
            if (watchParticle == i1) ParticleWatchPoint->set(particleWatchStep, i2, amp * d);
            else if (watchParticle == i2) ParticleWatchPoint->set(particleWatchStep, i1, -amp * d);
        }
    }
    return Success;
}

/*Vector Calculation::calcF2(const Vector& d1, const Vector& d2, Vector& F1)
{
    int sortOrder[Vector::dimension];
    d2.getSortOrder(sortOrder);

}*/

int Calculation::getPartnerBindingIndex(const Particle* const P, const int index)
{
    switch (index)
    {
        case 0:
            if (P->NB > 1) return 1;
            break;
        case 1:
            return 0;
        case 2:
            if (P->NB > 3) return 3;
            break;
        case 3:
            return 2;
        case 4:
            if (P->NB > 5) return 5;
            break;
        case 5:
            return 4;
        default:
            return -1;
    }
    return -1;
}

void Calculation::removeBinding(Particle *const part, const int index) const
{
    for (int n = index + 1; n < part->NB; ++n) part->bound[n-1] = part->bound[n];
    --(part->NB);
    part->bound[part->NB].p = nullptr;
}

void Calculation::addCandidate(Particle *const currPart, Particle *const candidate, const double dist) const
{
    int n, m;
    if (currPart->NB == currPart->MNB)
    {
        for (n=0; n < currPart->MNB; ++n) if (currPart->bound[n].lastDist > dist) break;
        if (n == currPart->MNB) return;
    }
    for (n = currPart->NC - 1; n >= 0; n--) if (currPart->candidates[n].lastDist < dist) break;
    if (n+1 == Particle::NCandidates) return;
    for (m=n+1; m < currPart->NC && currPart->candidates[m].lastDist == dist; ++m) if (currPart->candidates[m].p == candidate) return;
    if (currPart->NC < Particle::NCandidates) ++(currPart->NC);
    for (m = currPart->NC - 2; m > n; --m) currPart->candidates[m+1] = currPart->candidates[m];
    currPart->candidates[n+1].p = candidate;
    currPart->candidates[n+1].lastDist = dist;
}

void Calculation::correctEnergy()
{
    double T = getKineticEnergy(), V = getPotentialEnergy(), curEnergy = T + V, delta = Energy - curEnergy;
    // printf("Energy=%g, delta=%g\n", curEnergy, delta);
    //printf("OverallDeltaE=%g, lastDeltaE=%g, overallByUpdateBinding=%g, lastByUpdateBinding=%g\n", curEnergy - Energy, curEnergy - mLastEnergy, mAbsEDevUB += mCurEDevUB, mCurEDevUB);
    if (nullptr != mEnergyCsvLog) *mEnergyCsvLog << "\t" << (curEnergy - Energy) << "\t" << (curEnergy - mLastEnergy) << "\t" << (mAbsEDevUB += mCurEDevUB) << "\t" << mCurEDevUB << "\n";
    mLastEnergy = curEnergy;
    /*if (abs(delta / Energy) > 0.01)
    {
        setEnergy(T, V, delta);
        T = getKineticEnergy();
    }*/
    emit EnergiesChanged(T, T+V);
}

double Calculation::getE(Particle * const cP, const Vector &R, const bool useLastPos, const bool collectCandidates) const
{
    Particle* P2;
    int lx, ly, lz;
    double E = 0.0;
    Positions pos(collectCandidates ? particles : (useLastPos ? lastPos : currentPos));
    for (lz = ((lz = cP->zp - GridSizeDiv) >= 0 ? lz : 0); lz < ZS && lz <= cP->zp + GridSizeDiv; lz++)
        for (ly = ((ly = cP->yp - GridSizeDiv) >= 0 ? ly : 0); ly < YS && ly <= cP->yp + GridSizeDiv; ly++)
            for (lx = ((lx = cP->xp - GridSizeDiv) >= 0 ? lx : 0); lx < XS && lx <= cP->xp + GridSizeDiv; lx++)
                for (P2 = G[lx][ly][lz]; P2 != 0; P2 = P2->next) if (P2 != cP) if (getU(cP, P2, E, &R, pos, nullptr, collectCandidates) == Error) return Error_Double;
    return 0.5 * E;
}

double Calculation::getKineticEnergy() const
{
    int n;
    double T;
    for (n = 0, T = 0.0; n < N; n++)
        T += P[n].v.lengthSquared();
    return 0.5 * T;
}

double Calculation::getPotentialEnergy() const
{
    int mx, my, mz, lx, ly, lz, n;
    double U = 0.0;
    Particle *PP1, *PP2;
    for (n=0; n < NumPot; ++n) if (Pot[n] == nullptr || dPdR[n] == nullptr) return -1.0;
    for (mz = 0; mz < ZS; mz++) for (my = 0; my < YS; my++) for (mx = 0; mx < XS; mx++)
    {
        if (G[mx][my][mz] != 0)
        {
            //printf("l0\n");
            for (PP1 = G[mx][my][mz]; PP1->next != 0; PP1 = PP1->next)
                for (PP2 = PP1->next; PP2 != 0; PP2 = PP2->next) if (getU(PP1, PP2, U, nullptr, particles, nullptr, false) == Error) return Error_Double;
            //printf("l1\n");
            for (lz = mz; lz < ZS && lz <= mz + GridSizeDiv; lz++)
                for (ly = (lz > mz ? ((ly = my - GridSizeDiv) >= 0 ? ly : 0) : my);
                        ly < YS && ly <= my + GridSizeDiv; ly++)
                    for (lx = (ly > my || lz > mz ? ((lx = mx - GridSizeDiv) >= 0 ? lx : 0) : mx + 1);
                            lx < XS && lx <= mx + GridSizeDiv; lx++)
                        for (PP2 = G[lx][ly][lz]; PP2 != 0; PP2 = PP2->next)
                            for (PP1 = G[mx][my][mz]; PP1 != 0; PP1 = PP1->next) if (getU(PP1, PP2, U, nullptr, particles, nullptr, false) == Error) return Error_Double;
        }
    }
    return U;
}

void Calculation::getScales(double& rScF, int& rMaxX, int&rMaxY, int& rMaxZ)
{
	rScF = ScF;
    rMaxX = MaxX;
    rMaxY = MaxY;
	rMaxZ = MaxZ;
}

void Calculation::getSize(int& width, int& height)
{
	width = int(MaxX * ScF) + 20;
	height = int(MaxY * ScF) + 50;
}

double Calculation::getSpeed()
{
	return Speed;
}

double Calculation::getStepSize()
{
	return h;
}

void Calculation::initialize()
{
	//printf("Calculation::initialize\n");
	int n, x, y, z, nwp=0;
	double rx, rz, rxs = MaxX / double(PXS), rzs = rxs * sin(M_PI / 3.0);
    mLayerDistance = 2 * MaxY / double(PYS);
    mWavePhase = 0.0;
    mLSH = h;
	double y1 = 0.5 * (MaxY - mLayerDistance);
	double y2 = 0.5 * (MaxY + mLayerDistance);
    Vector F(double(XS) / MaxX, double(YS) / MaxY, double(ZS) / MaxZ * 1.001);
    bool even = false;
	for (x=0; x < XS; x++) for (y=0; y < YS; y++) for (z=0; z < ZS; z++) G[x][y][z] = 0;
	for (n=z=0, rz = 0.5 * rzs; z < PZS; z++, rz += rzs) 
	{
        even = (even ? false : true);
		for (x=0, rx = (even ? 0.5 * rxs : rxs); x < (even ? PXS: PXS - 1); x++, rx += rxs)
		{
            Vector r1(rx, y1, rz), r2(rx, y2, rz);
            initializeParticle(P[n++], x, z, r1, F, even);
            initializeParticle(P[n++], x, z, r2, F, even);
            if (P[n-1].Fixed)
            {
                P[n-2].WallPosIndex = nwp;
                FixedWallPos[nwp++] = P + n - 1;
            }
		}
	}
    calcMAR();
    for (n=0; n<N; ++n) for (x=0; x<2 && MAR[n][x].R != RM; ++x)
    {
        z = MAR[n][x].index;
        if (z <= n) continue;
        for (y=0; y<2 && MAR[z][y].R != RM && MAR[z][y].index != n; ++y) ;
        if (y<2 && MAR[z][y].index == n)
        {
            P[n].bound[P[n].NB].p = P+z;
            P[n].bound[P[n].NB++].lastDist = MAR[n][x].R;
            P[z].bound[P[z].NB].p = P+n;
            P[z].bound[P[z].NB++].lastDist = MAR[n][x].R;
        }
    }
    for (n=0; n<N; ++n) for (x=0; x < Particle::BoundAL && MAR[n][x].R != RM; ++x)
    {
        z = MAR[n][x].index;
        if (z <= n || MAR[n][x].R > 1.01 * rxs) continue;
        for (y=0; y < Particle::BoundAL && MAR[z][y].R != RM && MAR[z][y].index != n; ++y) ;
        if (y < Particle::BoundAL && MAR[z][y].index == n && P[n].bound[0].p != P+z && P[n].bound[1].p != P+z)
        {
            P[n].bound[P[n].NB].p = P+z;
            P[n].bound[P[n].NB++].lastDist = MAR[n][x].R;
            P[z].bound[P[z].NB].p = P+n;
            P[z].bound[P[z].NB++].lastDist = MAR[n][x].R;
        }
    }
    updateBindingPairs();
}

void Calculation::initializeParticle(Particle &cP, const int x, const int z, const Vector &iR, const Vector &Fact, const bool even) const
{
    int n = &cP - P, b;
    cP.lR = cP.R = iR;
    cP.lv = cP.v = Vector();
    cP.xp = ((b = int(Fact.X() * iR.X())) >= 0 ? (b < XS ? b : XS - 1) : 0);
    cP.yp = ((b = int(Fact.Y() * iR.Y())) >= 0 ? (b < YS ? b : YS - 1) : 0);
    cP.zp = ((b = int(Fact.Z() * iR.Z())) >= 0 ? (b < ZS ? b : ZS - 1) : 0);
	cP.prev = 0;
	cP.next = G[cP.xp][cP.yp][cP.zp];
    cP.NB = 0;
    memset(cP.bound, 0, Particle::BoundAL * sizeof(Particle*));
	G[cP.xp][cP.yp][cP.zp] = D[n] = &cP;
	if (cP.next != 0) cP.next->prev = D[n];
	if (z == 0 || x == 0 || z == PZS - 1 || x == (even ? PXS - 1 : PXS - 2))
	{
		cP.Fixed = true;
        cP.MNB = (z==0 || z == PZS - 1 ? (x==0 || x == PXS - 1 ? Particle::BoundAL - 4 : Particle::BoundAL - 2) : (even ? Particle::BoundAL - 3 : Particle::BoundAL - 1));
        if (x==0) cP.WaveParticle = true;
	}
    else
    {
        cP.Fixed = false;
        cP.MNB = Particle::BoundAL;
    }
}

void Calculation::setLayerDistance(double newDistance)
{
    const Vector diff(0.0, 0.5 * (newDistance - mLayerDistance), 0.0);
    for (int n=0; n < N-1; n+=2)
    {
        P[n].R -= diff;
        P[n+1].R += diff;
    }
    mLayerDistance = newDistance;
}

void Calculation::move()
{
	if (Move) Move = false;
	else Move = true;
}

void Calculation::run()
{
    // Contains the rk4 algorithm from Numerical Recipes, Third Edition
    mErrorCode = ECSuccess;
    for (int n=0; n < NumPot; ++n)
    {
        if (Pot[n] == nullptr)
        {
            mErrorCode = ECPotentialNotAvailable;
            qCritical() << "Pot[" << n << "] == NULL!";
            return;
        }
        if (dPdR[n] == nullptr)
        {
            mErrorCode = ECGradientNotAvailable;
            qCritical() << "dPdR[" << n << "] == NULL!";
            return;
        }
        if (!potentialOK[n])
        {
            mErrorCode =  ECPotentialNotOK;
            qCritical() << "Pot[" << n << "] is not OK!";
            return;
        }
    }
    int n, i;
    Vector *a = new Vector[N], *dm = new Vector[N], *dvm = new Vector[N], *dt = new Vector[N], *t0 = new Vector[N];
    Vector *dvt = new Vector[N];
	Particle *PB;
    if (!mEnergyCsvLogFilename.isEmpty())
    {
        mEnergyCsvLogFile = new QFile(mEnergyCsvLogFilename);
        mEnergyCsvLogFile->open(QIODevice::WriteOnly);
        mEnergyCsvLog = new QTextStream(mEnergyCsvLogFile);
        *mEnergyCsvLog << "Iteration\tOverallDeltaE\tlastDeltaE\toverallByUpdateBinding\tlastByUpdateBinding\n";
    }
    Energy = getKineticEnergy() + getPotentialEnergy();
	for (i=0, Run = true; Run && i != mMaxIt; i++)
	{
		/*if (i==30178)
        {
            printf("Break!");
        }*/

        rk4(t0, dvt, a, dt, dm, dvm, h);
        if (mErrorCode != ECSuccess) return;

        if (writeSnapShot)
        {
            WriteSnapshot();
            writeSnapShot = false;
        }
        if (watchParticle >= 0)
            ParticleWatchPoint->setSum(particleWatchStep, dvt[watchParticle] - Vector(ParticleWatchPoint->getSumX(1), ParticleWatchPoint->getSumY(1), ParticleWatchPoint->getSumZ(1)));
        /*if (i==147)
        {
            printf("Break!");
        }*/
        // printf("iteration=%d, ", i);
        updateBindings();
        //*DebugLog << i;
        if (nullptr != mEnergyCsvLog) *mEnergyCsvLog << i;
        correctEnergy();
        if (mInstanceId == -1)
        {
            if (mRotationChanged)
            {
                int *Sort = utils::heapSort(yz_SortFunctor(P, (rotated ? yz_SortFunctor::yzY : yz_SortFunctor::yzZ)), N);
                for (n=0; n<N; ++n) D[Sort[n]] = P + n;
                mRotationChanged = false;
            }
            else for (PB = P; PB != 0; ) for (n=1, PB = 0; n<N; n++) if (rotated ? D[n]->R.Y() < D[n-1]->R.Y() : D[n]->R.Z() < D[n-1]->R.Z())
            {
                PB = D[n];
                D[n] = D[n-1];
                D[n-1] = PB;
            }
            mutex.lock();
            if (rotated) for (n=0; n < N; n++) Pos[n] = Vector(D[n]->R.X(), D[n]->R.Z(), D[n]->R.Y());
            else for (n=0; n < N; n++) Pos[n] = D[n]->R;
            mutex.unlock();
            emit PictureChanged(Pos, N);
        }
        else sendCalcResult(i);
        if (watchParticle >= 0)
        {
            particleWatchStep = -1;
            break;
        }
	}
    delete[] a;
    delete[] dm;
    delete[] dvm;
    delete[] dt;
    delete[] dvt;
    delete[] t0;
}

void Calculation::sendCalcResult(const int iteration)
{
    static const int CPI[8] = {378, 379, 380, 381, 418, 419, 420, 421};
    static const double y1 = 0.5 * (MaxY - mLayerDistance), y2 = 0.5 * (MaxY + mLayerDistance);
    double currentDev = 0.0;
    for (int i=0; i<8; i+=2) currentDev += (P[CPI[i]].R.Y() - y1);
    for (int i=1; i<8; i+=2) currentDev += (P[CPI[i]].R.Y() - y2);
    currentDev *= 0.125;
    if (abs(currentDev) > abs(mMaxCalcResult)) mMaxCalcResult = currentDev;
    emit CalcState(mInstanceId, iteration, currentDev, mMaxCalcResult);
}

void Calculation::updateBlock(int n)
{
    int x, y, z;
    getGridAtPos(P[n].R, x, y, z);
    if (P[n].xp != x || P[n].yp != y || P[n].zp != z)
    {
        if (P[n].next != 0) P[n].next->prev = P[n].prev;
        if (P[n].prev != 0) P[n].prev->next = P[n].next;
        else G[P[n].xp][P[n].yp][P[n].zp] = P[n].next;
        P[n].prev = 0;
        P[n].next = G[x][y][z];

        if (G[x][y][z] != 0) G[x][y][z]->prev = P + n;
        G[x][y][z] = P + n;
        P[n].xp = x;
        P[n].yp = y;
        P[n].zp = z;
    }
}

void Calculation::getGridAtPos(const Vector& Pos, int& x, int& y, int& z) const
{
    static const Vector F(double(XS) / MaxX, double(YS) / MaxY, double(ZS) / MaxZ);
    x = ((x = int(F.X() * Pos.X())) >= 0 ? (x < XS ? x : XS - 1) : 0);
    y = ((y = int(F.Y() * Pos.Y())) >= 0 ? (y < YS ? y : YS - 1) : 0);
    z = ((z = int(F.Z() * Pos.Z())) >= 0 ? (z < ZS ? z : ZS - 1) : 0);
}

void Calculation::rk4(Vector *t0, Vector *dvt, Vector *a, Vector *dt, Vector* dm, Vector* dvm, const double lh)
{
    double hh = 0.5 * lh, h6 = lh / 6.0;
    int n;
    // int *debugNullPtr = nullptr;
    Result result = Success;
    for (int i=1; i==1; ++i)
    {
        if (Move)
        {
            //applyMove(lh);
            applyWave(lh);
        }
        for (n=0; n<N; n++)
        {
            updateBlock(n);
            t0[n] = P[n].R;
        }
        for (n=0, U = T = 0.0; n<N; n++) T += P[n].v.lengthSquared();
        if (watchParticle >= 0) particleWatchStep = 0;
        result = geta(t0, a, false);
        if (result == Error)
        {
            for (n=0; n<N; ++n)
            {
                P[n].v = P[n].lv;
                P[n].R = P[n].lR;
            }
            break;
        }
        if (E == 0.0 || T == 0.0 || Move) E = T + U;
        for (n=0; n<N; n++) if (!P[n].Fixed)
        {
            t0[n] = P[n].R + hh * P[n].v;
            dt[n] = P[n].v + hh * a[n];
        }
        for (n=0, U = T = 0.0; n<N; n++) T += dt[n].lengthSquared();
        if (watchParticle >= 0) ParticleWatchPoint->setSum(particleWatchStep++, a[watchParticle]);
        result = geta(t0, dvt, false);
        if (result == Error) break;
        for (n=0; n<N; n++) if (!P[n].Fixed)
        {
            t0[n] = P[n].R + hh * dt[n];
            dm[n] = P[n].v + hh * dvt[n];
        }
        for (n=0, U = T = 0.0; n<N; n++) T += dm[n].lengthSquared();
        if (watchParticle >= 0) ParticleWatchPoint->setSum(particleWatchStep++, dvt[watchParticle]);
        result = geta(t0, dvm, false);
        if (result == Error) break;
        for (n=0; n<N; n++) if (!P[n].Fixed)
        {
            t0[n] = P[n].R + lh * dm[n];
            dm[n] += dt[n];
            dt[n] = P[n].v + lh * dvm[n];
            dvm[n] += dvt[n];
        }
        for (n=0, U = T = 0.0; n<N; n++)
        {
            T += dt[n].lengthSquared();
            P[n].NC = 0;
        }
        if (watchParticle >= 0) ParticleWatchPoint->setSum(particleWatchStep++, dvm[watchParticle]);
        result = geta(t0, dvt, true);
        if (result == Error) break;
        for (n=0; n<N; n++)
        {
            if (P[n].Fixed) P[n].lv = P[n].lR;
            P[n].lR = P[n].R;
            if (!P[n].Fixed)
            {
                P[n].lv = P[n].v;
                P[n].R += h6 * (P[n].v + dt[n] + 2.0 * dm[n]);
                P[n].v += h6 * (a[n] + dvt[n] + 2.0 * dvm[n]);
                if ((P[n].R.X() < 0.0 && P[n].v.X() < 0.0) || (P[n].R.X() > MaxX && P[n].v.X() > 0.0)) P[n].v *= Vector(-1.0, 1.0, 1.0);
                if ((P[n].R.Y() < 0.0 && P[n].v.Y() < 0.0) || (P[n].R.Y() > MaxY && P[n].v.Y() > 0.0)) P[n].v *= Vector(1.0, -1.0, 1.0);
                if ((P[n].R.Z() < 0.0 && P[n].v.Z() < 0.0) || (P[n].R.Z() > MaxZ && P[n].v.Z() > 0.0)) P[n].v *= Vector(1.0, 1.0, -1.0);
                if (isnan(P[n].R.X()) || isnan(P[n].R.Y()) || isnan(P[n].R.Z()) || isnan(P[n].v.X()) || isnan(P[n].v.Y()) || isnan(P[n].v.Z()))
                {
                    //*debugNullPtr = 5;
                    printf("After calculation of new position and v: Particel %d is nan!\n", n);
                    Run = false;
                }
            }
        }
        mLSH = lh;
    }
    if (result == Error)
    {
        Run = false;
        /* const double nh = 0.5 * lh;
        if (nh < 1e-10)
        {
            mErrorCode = ECParticlesTooClose;
            //printf("break!");
        }
        else
        {
            for (int n=0; n<N; ++n) if (P[n].Fixed) P[n].R = P[n].lR;
            mWavePhase = mLastWavePhase;
            rk4(t0, dvt, a, dt, dm, dvm, nh);
            if (mErrorCode == ECSuccess) rk4(t0, dvt, a, dt, dm, dvm, nh);
            if (mErrorCode != ECSuccess && mLSH == lh)
            {
                qInfo() << "Second stage of error correction!";
                printf("Second stage\n");
                mErrorCode = ECSuccess;
                mWavePhase = mSecondToLastWavePhase;
                for (int n=0; n<N; ++n)
                {
                    if (P[n].Fixed) P[n].R = P[n].lv;
                    else
                    {
                        P[n].R = P[n].lR;
                        P[n].v = P[n].lv;
                    }

                }
                rk4(t0, dvt, a, dt, dm, dvm, nh);
                if (mErrorCode == ECSuccess) rk4(t0, dvt, a, dt, dm, dvm, nh);
                if (mErrorCode == ECSuccess) rk4(t0, dvt, a, dt, dm, dvm, nh);
                if (mErrorCode == ECSuccess) rk4(t0, dvt, a, dt, dm, dvm, nh);
            }
        }*/
    }
}

void Calculation::applyMove(const double lh)
{

    double XMid = 0.5 * MaxX, R, dX, dZ, ZMid = 0.5 * MaxZ, dist = MaxX / double(PXS);
    for (int n=0; n<N; n++) if (P[n].Fixed)
    {
        P[n].lR = P[n].R;
        dX = XMid - P[n].R.X();
        dZ = ZMid - P[n].R.Z();
        R = lh * Speed / (dX * dX + dZ * dZ);
        Vector step(dX * R, 0.0, dZ * R);
        if (P[n].R.Y() > YMid) P[n].R += step;
        else P[n].R -= step;
        if (P[n].R.X() < 0.0) doParticleLayerSwitch(P+n, step * (dist / step.X()));
        else if (P[n].R.Z() < 0.0) doParticleLayerSwitch(P+n, step * (dist / step.Z()));
        else if (P[n].R.X() > MaxX) doParticleLayerSwitch(P+n, step * (-dist / step.X()));
        else if (P[n].R.Z() > MaxZ) doParticleLayerSwitch(P+n, step * (-dist / step.Z()));
    }
}

void Calculation::applyWave(const double lh)
{
    const double newWavePhase = mWavePhase + waveStep * lh, sinNWP = sin(newWavePhase), sinWP = sin(mWavePhase);
    if (newWavePhase > M_PI)
    {
        Move = false;
        return;
    }
    const Vector cAmp(0.0, waveAmp * (sinNWP * sinNWP - sinWP * sinWP), 0.0 /* waveAmp * sinNWP * sinNWP /  M_PI */);
    for (int n=0; n<N; ++n) if (P[n].WaveParticle) P[n].R += cAmp;
    mSecondToLastWavePhase = mLastWavePhase;
    mLastWavePhase = mWavePhase;
    mWavePhase = newWavePhase;
}

bool Calculation::isBoundToWaveParticle(const Particle *const P)
{
    if (P->Fixed) return false;
    for (int n=0; n < P->NB; ++n) if (P->bound[n].p->WaveParticle) return true;
    return false;
}

void Calculation::doParticleLayerSwitch(Particle *const cP, const Vector& Dist)
{
    int *debugNullPtr = nullptr;
    Vector Pos = cP->R + Dist;
    double dist = Dist.length(), cDist;
    Particle* newF = getClosestBound(cP, Pos);
    if (nullptr != newF)
    {
        newF->R = Pos;
        newF->Fixed = true;
        newF->WallPosIndex = cP->WallPosIndex;
    }
    for (int m = cP->NB - 1; m >= 0; --m)
    {
        if (cP->bound[m].p->Fixed) --(cP->bound[m].p->MNB);
        removeBinding(cP->bound[m].p, getBindingIndexAtBound(cP, m));
        removeBinding(cP, m);
    }
    if (-1 == cP->WallPosIndex) *debugNullPtr = 5;
    Particle* newFree = FixedWallPos[cP->WallPosIndex];
    newFree->Fixed = false;
    if (Particle::BoundAL == newFree->NB) bindToRadical(cP, newFree, dist, true);
    else
    {
        newFree->MNB = Particle::BoundAL;
        newFree->bound[newFree->NB++].p = cP;
        cP->bound[0].p = newFree;
        newFree->WallPosIndex = -1;
    }
    cP->R = newFree->R - Dist;
    cP->NB = cP->MNB = 1;
    if (isBindingDoubled(newFree - P)) *debugNullPtr = 5;
    const double MaxDev = 2.0 * dist;
    for (int i=0; i < newFree->NB; ++i) if (cP != newFree->bound[i].p) for (int j=0; j < newFree->bound[i].p->NB; ++j)
        if (newFree->bound[i].p->bound[j].p->Fixed && cP != newFree->bound[i].p->bound[j].p && isNotBound(cP, newFree->bound[i].p->bound[j].p)
            && (cDist = (cP->R - newFree->bound[i].p->bound[j].p->R).length()) < MaxDev)
    {
        if (Particle::BoundAL == cP->MNB) break;
        Particle* newBound = newFree->bound[i].p->bound[j].p;
        if (Particle::BoundAL > newBound->NB)
        {
            cP->bound[cP->NB].p = newBound;
            newBound->bound[newBound->NB].p = cP;
            newBound->MNB = ++(newBound->NB);
            cP->MNB = ++(cP->NB);
        }
        else if (bindToRadical(cP, newBound, cDist, false)) ++(cP->MNB);
        if (isBindingDoubled(newBound - P)) *debugNullPtr = 5;
    }
    if (isBindingDoubled(cP - P)) *debugNullPtr = 5;
}

int Calculation::getBindingIndexAtBound(const Particle *const P1, const int index)
{
    for (int n=0; n < P1->bound[index].p->NB; ++n) if (P1->bound[index].p->bound[n].p == P1) return n;
    return -1;
}

Particle* Calculation::getClosestBound(const Particle *const P1, const Vector& Pos)
{
    Particle* Res = nullptr;
    double minD = -1.0;
    for (int n=0; n < P1->NB; ++n) if (!P1->bound[n].p->Fixed)
    {
        if (Res == nullptr) Res = P1->bound[n].p;
        else
        {
            if (minD < 0.0) minD = (Res->R - Pos).lengthSquared();
            double D = (P1->bound[n].p->R - Pos).lengthSquared();
            if (D < minD)
            {
                minD = D;
                Res = P1->bound[n].p;
            }
        }
    }
    if (Res == nullptr) for (int n=0; n < P1->NB; ++n) for (int m=0; m < P1->bound[n].p->NB; ++m) if (!P1->bound[n].p->bound[m].p->Fixed)
    {
        if (Res == nullptr) Res = P1->bound[n].p->bound[m].p;
        else
        {
            if (minD < 0.0) minD = (Res->R - Pos).lengthSquared();
            double D = (P1->bound[n].p->bound[m].p->R - Pos).lengthSquared();
            if (D < minD)
            {
                minD = D;
                Res = P1->bound[n].p->bound[m].p;
            }
        }
    }
    return Res;
}


bool Calculation::UpdateBindings()
{
    for (int n=0; n<N; ++n) P[n].NC = 0;
    for (int mz = 0; mz < ZS; ++mz) for (int my = 0; my < YS; ++my) for (int mx = 0; mx < XS; ++mx)
        for (Particle *PP1 = G[mx][my][mz]; PP1 != 0; PP1 = PP1->next) getE(PP1, PP1->R, false, true);
    return updateBindings();
}

bool Calculation::updateBindings()
{
    int *debugNullPtr = nullptr;
    bool rValue = true;
    int* randomOrder = createRandomParticleOrder();
    for (int n=0; n<N; ++n)
    {
        Particle* CP = P + randomOrder[n];
        for (int m=0; m < CP->NC; ++m)
        {
            Particle* CanP = CP->candidates[m].p;
            if (isNotBound(CP, CanP))
            {
                //qInfo() << "CP=" << CP-P << ", CanP=" << CanP-P << ", CP->NB=" << CP->NB << ", CanP->NB=" << CanP->NB << ", CP->bound[0]=" << CP->bound[0].p-P << ", CP->bound[1]=" << CP->bound[1].p-P
                  //      << "CP->bound[2]=" << CP->bound[2].p-P << ", CP->bound[3]=" << CP->bound[3].p-P << ", CanP->bound[0]=" << CanP->bound[0].p-P << ", CanP->bound[1]=" << CanP->bound[1].p-P
                    //    << "CanP->bound[2]=" << CanP->bound[2].p-P << "CanP->bound[3]=" << CanP->bound[3].p-P;
                if (CP->NB < CP->MNB)
                {
                    if (CanP->NB < CanP->MNB)
                    {
                        double oldBE = 0.0, newBE = 0.0;
                        getU(CP, CanP, oldBE, nullptr, particles, nullptr, false);
                        CP->bound[CP->NB++] = CP->candidates[m];
                        CanP->bound[CanP->NB].p = CP;
                        CanP->bound[CanP->NB++].lastDist = CP->candidates[m].lastDist;
                        getU(CP, CanP, newBE, nullptr, particles, nullptr, false);
                        mCurEDevUB += newBE - oldBE;
                        if (isBindingDoubled(CP - P)) *debugNullPtr = 5;
                        if (isBindingDoubled(CanP - P)) *debugNullPtr = 5;
                    }
                    else bindToRadical(CP, CanP, CP->candidates[m].lastDist, false);
                }
                else
                {
                    //qInfo() << "1_CP->bound[0]=" << CP->bound[0].p-P << ", CP->bound[1]=" << CP->bound[1].p-P << "CP->bound[2]=" << CP->bound[2].p-P << ", CP->bound[3]=" <<
                      //                    CP->bound[3].p-P;
                    if (CanP->NB < CanP->MNB) bindToRadical(CanP, CP, CP->candidates[m].lastDist, false);
                    else
                    {
                        //qInfo() << "2_CP->bound[0]=" << CP->bound[0].p-P << ", CP->bound[1]=" << CP->bound[1].p-P << "CP->bound[2]=" << CP->bound[2].p-P << ", CP->bound[3]=" <<
                          //                CP->bound[3].p-P;
                        for (int i=0; i < CanP->NB; ++i) if (goCentral(CP, CanP, m, i)) break;
                    }
                }
            }
        }
    }
    verifyNoBindingDoubled();
    updateBindingPairs();
    return rValue;
}

bool Calculation::goCentral(Particle *const CP, Particle *const CanP, const int CPCanIndex, const int CanPBIndex)
{
    int *debugNullPtr = nullptr, n = CP->NB, CPBoundOrder[Particle::BoundAL], b=1;
    //printf("goCentral: CP=%ld, CanP=%ld, CPCanIndex=%d, CanP->bound[CanPBIndex]=%ld\n", CP-P, CanP-P, CPCanIndex, CanP->bound[CanPBIndex].p-P);
    for (int i=0; i < CP->NB; ++i)
    {
        if (CP->bound[i].p == CanP->bound[CanPBIndex].p) n=i;
        CPBoundOrder[i] = i;
    }
    while (b>=0)
    {
        b=-1;
        for (int i=1; i < CP->NB; ++i) if (CP->bound[CPBoundOrder[i]].lastDist > CP->bound[CPBoundOrder[i-1]].lastDist)
        {
            b = CPBoundOrder[i];
            CPBoundOrder[i] = CPBoundOrder[i-1];
            CPBoundOrder[i-1] = b;
        }
    }
    if (n < CP->NB)
    {
        Particle *P3 = CP->bound[n].p; 
        //printf("P3=%ld\n", P3-P);
        int m, j, k=-1;
        double candidatesDist;
        for (m=0; m < CP->NB; ++m)
        {
            int i = CPBoundOrder[m];
            Particle *P4 = nullptr;
            for (j=0; j < P3->NC; ++j) if (CP->bound[i].p == P3->candidates[j].p && isNotBound(P3, P3->candidates[j].p))
            {
                P4 = CP->bound[i].p;
                candidatesDist = CP->candidates[CPCanIndex].lastDist + P3->candidates[j].lastDist;
                break;
            }
            if (P4 == nullptr) for (k=0; k < CP->bound[i].p->NC; ++k) if (CP->bound[i].p->candidates[k].p == P3 && isNotBound(CP->bound[i].p, P3))
            {
                P4 = CP->bound[i].p;
                candidatesDist = CP->candidates[CPCanIndex].lastDist + P4->candidates[k].lastDist;
                break;
            }
            if (P4 != nullptr && candidatesDist < CP->bound[n].lastDist + CanP->bound[CanPBIndex].lastDist)
            {
                int l, o;
                double oldBE = 0.0, newBE = 0.0;
                getU(CP, CanP, oldBE, nullptr, particles, nullptr, false);
                getU(CP, P4, oldBE, nullptr, particles, nullptr, false);
                getU(CanP, P3, oldBE, nullptr, particles, nullptr, false);
                getU(P3, P4, oldBE, nullptr, particles, nullptr, false);
                for (l=0; l < P3->NB; ++l) if (P3->bound[l].p == CanP) break;
                for (o=0; o < P4->NB; ++o) if (P4->bound[o].p == CP) break;
                CP->bound[i].lastDist = CanP->bound[CanPBIndex].lastDist = CP->candidates[CPCanIndex].lastDist;
                P3->bound[l].lastDist = P4->bound[o].lastDist = (j < P3->NC ? P3->candidates[j].lastDist : P4->candidates[k].lastDist);
                //qInfo() << "7_CP->bound[0]=" << CP->bound[0].p-P << ", CP->bound[1]=" << CP->bound[1].p-P << "CP->bound[2]=" << CP->bound[2].p-P << ", CP->bound[3]=" <<
                  //                      CP->bound[3].p-P << ", j=" << j << ", NC=" << CanP->bound[CanPBIndex].p->NC;
                //printf("goCentral: %ld-%ld, %ld-%ld, %ld-%ld and %ld-%ld to %ld-%ld, %ld-%ld, %ld-%ld and %ld-%ld\n",
                  //     CP-P, CP->bound[i].p-P, CanP-P, CanP->bound[CanPBIndex].p-P, P3-P, P3->bound[l].p-P, P4-P, P4->bound[o].p-P,
                    //   CP-P, P3->bound[l].p-P, CanP-P, P4->bound[o].p-P, P3-P, CP->bound[i].p-P, P4-P, CanP->bound[CanPBIndex].p-P);
                std::swap(CP->bound[i].p, P3->bound[l].p);
                std::swap(CanP->bound[CanPBIndex].p, P4->bound[o].p);
                getU(CP, CanP, newBE, nullptr, particles, nullptr, false);
                getU(CP, P4, newBE, nullptr, particles, nullptr, false);
                getU(CanP, P3, newBE, nullptr, particles, nullptr, false);
                getU(P3, P4, newBE, nullptr, particles, nullptr, false);
                mCurEDevUB += newBE - oldBE;
                if (isBindingDoubled(CP - P)) *debugNullPtr = 5;
                if (isBindingDoubled(CanP - P)) *debugNullPtr = 5;
                if (isBindingDoubled(P3 - P)) *debugNullPtr = 5;
                if (isBindingDoubled(P4 - P)) *debugNullPtr = 5;
                //qInfo() << "Break k, j=" << j << ", NC=" << CanP->bound[CanPBIndex].p->NC;
                return true;
            }
            // else if (P4 != nullptr) printf("P4=%ld, candidatesDist=%g, boundDist=%g\n", P4-P, candidatesDist, CP->bound[n].lastDist + CanP->bound[CanPBIndex].lastDist);
        }
    }
    return false;
}

bool Calculation::bindToRadical(Particle *const CP, Particle *const CanP, const double lastDist, const bool force)
{
    int leastBound = 0;
    int *debugNullPtr = nullptr;
    for (int i=1; i < CanP->NB; ++i) if (CanP->bound[i].lastDist > CanP->bound[leastBound].lastDist) leastBound = i;
    if (CanP->bound[leastBound].lastDist > lastDist || force)
    {
        Particle* LBP = CanP->bound[leastBound].p;
        int i;
        double oldBE = 0.0, newBE = 0.0;
        getU(CP, CanP, oldBE, nullptr, particles, nullptr, false);
        getU(CanP, LBP, oldBE, nullptr, particles, nullptr, false);
        for (i=0; i < LBP->NB; ++i) if (LBP->bound[i].p == CanP) break;
        if (i < LBP->NB) removeBinding(LBP, i);
        else leastBound = CanP->NB++;
        CanP->bound[leastBound].p = CP;
        CanP->bound[leastBound].lastDist = lastDist;
        CP->bound[CP->NB].p = CanP;
        CP->bound[CP->NB++].lastDist = lastDist;
        getU(CP, CanP, newBE, nullptr, particles, nullptr, false);
        getU(CanP, LBP, newBE, nullptr, particles, nullptr, false);
        mCurEDevUB += newBE - oldBE;
        if (isBindingDoubled(CP - P)) *debugNullPtr = 5;
        if (isBindingDoubled(CanP - P)) *debugNullPtr = 5;
        if (isBindingDoubled(LBP - P)) *debugNullPtr = 5;
        return true;
    }
    return false;
}

void Calculation::updateBindingPairs()
{
    Vector unit[Particle::BoundAL];
    for (int n=0; n<N; ++n) if (P[n].NB > 2)
    {
        for (int i=0; i < P[n].NB; ++i) unit[i] = (P[n].bound[i].p->R - P[n].R) / P[n].bound[i].lastDist;
        int bp1, bp2;
        double minValue = 1.0, value;
        for (int i=0; i < P[n].NB - 1; ++i) for (int j=i+1; j < P[n].NB; ++j) if ((value = unit[i].dot(unit[j])) < minValue)
        {
            value = minValue;
            bp1 = i;
            bp2 = j;
        }
        if (bp1 != 0 && bp2 != 0)
        {
            if (bp1 > 2)
            {
                std::swap(P[n].bound[0], P[n].bound[bp1]);
                bp1 = 0;
            }
            else
            {
                std::swap(P[n].bound[0], P[n].bound[bp2]);
                bp2 = 0;
            }
        }
        if (bp1 != 1 && bp2 != 1)
        {
            if (bp1 > 2) std::swap(P[n].bound[1], P[n].bound[bp1]);
            else std::swap(P[n].bound[1], P[n].bound[bp2]);
        }
    }
}

int* Calculation::createRandomParticleOrder()
{
    int *randArray = new int[N];
    for (int n=0; n<N; ++n) randArray[n] = n;
    for (int n=0; n<N; ++n)
    {
        int n1 = static_cast<int>(mRandPOF1 * rand()), n2 = static_cast<int>(mRandPOF2 * rand());
        if (n1 == n2) n2 = N-1;
        std::swap(randArray[n1], randArray[n2]);
    }
    return randArray;
}

double Calculation::dist(const Particle * const P1, const Particle * const P2)
{
    return (P1->R - P2->R).length();
}

bool Calculation::isNotBound(const Particle *const P1, const Particle *const P2)
{
    for (int i=0; i < P1->NB; ++i) if (P1->bound[i].p == P2) return false;
    return true;
}

void Calculation::WriteSnapshot()
{
    Particle* PC = new Particle[N];
    memcpy(PC, P, sizeof(Particle) * N);
    emit WriteSnapShot(PC, N);
}

void Calculation::rotate()
{
	if (rotated) rotated = false;
	else rotated = true;
    mRotationChanged = true;
}

double Calculation::setEnergy(const double deltaEnergy)
{
    return Energy = setEnergy(getKineticEnergy(), getPotentialEnergy(), deltaEnergy);
}

double Calculation::setEnergy(const double T, const double V, const double deltaEnergy)
{
	int n;
    double ParE, ParV;
	double VC, A1, A2, RD = M_PI / RAND_MAX, EnDiff = 2.0 * (deltaEnergy) / double(N);
    for (n=0; n < NumPot; ++n) if (Pot[n] == nullptr || dPdR[n] == nullptr) return -1.0;
    double Energy = E = T + V;
    if (deltaEnergy > 0.0)
	{
		Energy += deltaEnergy;
		VC = sqrt(EnDiff);
		for (n=0; n<N; n++)
		{
            if (P[n].v == Vector())
			{
				A1 = 2.0 * RD * rand();
				A2 = RD * rand();
                P[n].v = Vector(sin(A1) * sin(A2), sin(A1) * cos(A2), cos(A1)) * VC;
			}
			else
			{
                ParE = P[n].v.lengthSquared();
				ParV = sqrt((ParE + EnDiff) / ParE);
                P[n].v *= ParV;
			}
		}
	}
	else if (-deltaEnergy >= T)
    {
        for (int n=0; n<N; ++n) P[n].v.clear();
        Energy = V;
    }
    else
	{
        int *Sort = utils::heapSort(VSortFunctor(P), N), EOrder[N];
        for (n=0; n<N; ++n) EOrder[Sort[n]] = n;
        delete[] Sort;
        for (n=N-1; n>=0; --n)
        {
            double lT = P[EOrder[n]].v.lengthSquared();
            if (lT > -EnDiff) EnDiff = 2.0 * (deltaEnergy - Energy + E) / double(n+1);
            if (lT > -EnDiff) break;
            P[EOrder[n]].v.clear();
            Energy -= 0.5 * lT;
        }
        for (; n>=0; --n)
        {
            ParE = P[EOrder[n]].v.lengthSquared();
			ParV = sqrt((ParE + EnDiff) / ParE);
            P[EOrder[n]].v *= ParV;
        }
        Energy = E - deltaEnergy;
	}
	qInfo() << "SetEnergy: wanted Delta=" << deltaEnergy << ", current E=" << E << ", new Energy=" << Energy;
	E = Energy;
    if (isnan(E))
    {
        int *debugNullPtr = nullptr;
        *debugNullPtr = 5;
    }
    return Energy;
}

bool Calculation::setPotential(const PotRole role, PotStruct &PotS)
{
    if (Pot[role] != nullptr) delete[] Pot[role];
    if (dPdR[role] != nullptr) delete[] dPdR[role];
    double dRScale = 1.0 / PotS.getRZoom(), Rmin = Rm * dRScale, Rmax = RM * dRScale, devF = PotS.getVZoom() * dRScale;
    if (role == Angular)
    {
        Rmin = -1.0;
        Rmax = 1.0;
        devF = PotS.getVZoom();
    }
    Pot[role] = PotS.getPotential().getPoints(Rmin, Rmax, NPot);
    dPdR[role] = PotS.getPotential().get_dVdR(Rmin, Rmax, NPot);
    for (int n=0; n < NPot; ++n)
    {
        Pot[role][n] *= PotS.getVZoom();
        dPdR[role][n] *= devF;
        if (role == Angular && n>0 && n < NPot - 1)
        {
            const double cosP = static_cast<double>(n) * 2.0 / (NPot - 1) - 1.0;
            const double phi = acos(cosP);
            dPdR[role][n] *= sin(phi) *  cos(phi) * -2.0;
        }
    }
    if (role == Angular)
    {
        dPdR[role][0] = 0.0;
        dPdR[role][NPot - 1] = 0.0;
    }
    checkPotential(role);
    return potentialOK[role];
}

void Calculation::checkPotential(const PotRole role)
{
    if (role == ClosestTwo)
    {
        if (Pot[ClosestTwo] != nullptr) potentialOK[role] = true;
        checkPotentials();
    }
    else if (Pot[role] != nullptr)
    {
        potentialOK[role] = true;
        if (Pot[ClosestTwo]  != nullptr && role != Angular) for (int n=0; n < NPot && Pot[ClosestTwo][n] > 0.0; ++n) if (Pot[role][n] < Pot[ClosestTwo][n])
        {
            potentialOK[role] = false;
            break;
        }
    }
}

void Calculation::checkPotentials()
{
    for (int n = 1; n < NumPot; ++n) checkPotential(static_cast<PotRole>(n));
}

bool Calculation::arePotentialsOK()
{
    for (int n = ClosestTwo; n != NumPot; ++n) if (!potentialOK[n]) return false;
    return true;
}

void Calculation::setSpeed(double S)
{
	Speed = S;
}

void Calculation::setStepSize(double nh)
{
	mLSH = h = nh;
}

void Calculation::stop()
{
	Run = false;
}

void Calculation::triggerSnapShot()
{
    if (isRunning()) writeSnapShot = true;
    else WriteSnapshot();
}

Particle* Calculation::getParticles(int &num)
{
    num = N;
    return P;
}

void Calculation::CalcEndpointsOfEnergyDefinitionAxis(const int particeIndex, const Vector &direction, Vector &end1, Vector &end2) const
{
    const Vector& Point1 = P[particeIndex].R;
    if (direction.X() != 0.0)
    {
        end1 = Point1 - (Point1.X() / direction.X()) * direction;
        end2 = Point1 + ((MaxX - Point1.X()) / direction.X()) * direction;
    }
    else if (direction.Y() != 0.0)
    {
        end1 = Point1 - (Point1.Y() / direction.Y()) * direction;
        end2 = Point1 + ((MaxY - Point1.Y()) / direction.Y()) * direction;
    }
    else
    {
        end1 = Point1 - (Point1.Z() / direction.Z()) * direction;
        end2 = Point1 + ((MaxZ - Point1.Z()) / direction.Z()) * direction;
    }
    if (end1.Y() < 0.0) end1 = Point1 - (Point1.Y() / direction.Y()) * direction;
    else if (end1.Y() > MaxY) end1 = Point1 + ((MaxY - Point1.Y()) / direction.Y()) * direction;
    if (end2.Y() < 0.0) end2 = Point1 - (Point1.Y() / direction.Y()) * direction;
    else if (end2.Y() > MaxY) end2 = Point1 + ((MaxY - Point1.Y()) / direction.Y()) * direction;
    if (end1.Z() < 0.0) end1 = Point1 - (Point1.Z() / direction.Z()) * direction;
    else if (end1.Z() > MaxZ) end1 = Point1 + ((MaxZ - Point1.Z()) / direction.Z()) * direction;
    if (end2.Z() < 0.0) end2 = Point1 - (Point1.Z() / direction.Z()) * direction;
    else if (end2.Z() > MaxZ) end2 = Point1 + ((MaxZ - Point1.Z()) / direction.Z()) * direction;
}

void Calculation::GetAxisEnergies(PotentialDefinerInputData &data)
{
    const int particleIndex = data.getParticleIndex(), N = data.getNumnPoints();
    Particle& currentParticle = P[particleIndex];
    const Vector start = data.getEnd1(), end = data.getEnd2(), diff = end - start, particlePos = currentParticle.R;
    const double s = 1.0 / (N-1);
    Vector point = start;
    for (int n=0; n<N; ++n, point += s * diff)
    {
        currentParticle.R = point;
        updateBlock(particleIndex);
        if (UpdateBindings()) data.addBoundChange(n);
        const Particle* P2;
        double UnboundSum = 0.0;
        int lx, ly, lz;
        for (lz = ((lz = currentParticle.zp - GridSizeDiv) >= 0 ? lz : 0); lz < ZS && lz <= currentParticle.zp + GridSizeDiv; lz++)
            for (ly = ((ly = currentParticle.yp - GridSizeDiv) >= 0 ? ly : 0); ly < YS && ly <= currentParticle.yp + GridSizeDiv; ly++)
                for (lx = ((lx = currentParticle.xp - GridSizeDiv) >= 0 ? lx : 0); lx < XS && lx <= currentParticle.xp + GridSizeDiv; lx++)
                    for (P2 = G[lx][ly][lz]; P2 != 0; P2 = P2->next) if (P2 != &currentParticle)
        {
            Vector d = P2->R - currentParticle.R, b;
            int bi1=N, bi2=N;
            double r = d.length();
            int p = int((r - Rm) * potRangeScale);
            if (p < 0 || p >= NPot) continue;
            for (int m=0; m<4; ++m)
            {
                if (currentParticle.bound[m].p == P2) bi1=m;
                if (P2->bound[m].p == &currentParticle) bi2=m;
            }
            if (bi1 < Particle::BoundAL && bi2 < Particle::BoundAL)
            {
                switch (bi1)
                {
                case 0:
                    data.SetFirstBound(n, Pot[ClosestTwo][p]);
                    break;
                case 1:
                    data.SetSecondBound(n, Pot[ClosestTwo][p]);
                    break;
                case 2:
                    data.SetThirdBound(n, Pot[NextTwo][p]);
                    break;
                case 3:
                    data.SetFourthBound(n, Pot[NextTwo][p]);
                    break;
                case 4:
                    data.SetFifthBound(n, Pot[NextTwo][p]);
                    break;
                default:
                    data.SetSixthBound(n, Pot[NextTwo][p]);
                    break;
                }
            }
            else UnboundSum += Pot[Remaining][p];
        }
        data.SetUnbound(n, UnboundSum);
    }
    currentParticle.R = particlePos;
    updateBlock(particleIndex);
    UpdateBindings();
}

int Calculation::TranslateParticleIndex(int index) const
{
    for (int i=0; i < N; ++i) if (D[i] == P + index) return i;
    return -1;
}

void Calculation::verifyNoBindingDoubled() const
{
    for (int n=0; n<N; ++n) for (int i=0; i < P[n].NB - 1; ++i) for (int j=i+1; j < P[n].NB; ++j) if (P[n].bound[i].p == P[n].bound[j].p)
        qCritical() << "For P[" << n << "] is bound[" << i << "]==P[" << (P[n].bound[i].p-P) << "] and bound[" << j << "]==P[" <<  (P[n].bound[j].p-P) << "]!";
}

bool Calculation::isBindingDoubled(const int n) const
{
    if (Particle::BoundAL < P[n].NB || Particle::BoundAL < P[n].MNB) return true;
    for (int i=0; i < P[n].NB; ++i) if (nullptr == P[n].bound[i].p) return true;
    for (int i=0; i < P[n].NB - 1; ++i) for (int j=i+1; j < P[n].NB; ++j) if (P[n].bound[i].p == P[n].bound[j].p) return true;
    for (int i=0; i < P[n].NB; ++i) if (P[n].bound[i].p == P+n) return true;
    for (int i=0; i < P[n].NB; ++i)
    {
        int j;
        for (j=0; j < P[n].bound[i].p->NB; ++j) if (P[n].bound[i].p->bound[j].p == P+n) break;
        if (j == P[n].bound[i].p->NB) return true;
    }
    return false;
}

