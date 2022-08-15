#define _USE_MATH_DEFINES

#include "Calculation.h"
#include "heapsort.h"
#include "deltaesortfunctor.h"
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


Calculation::Calculation(PotStruct* PotSs, QObject* parent): QThread(parent), Error_Double(0.0/0.0), NPot(30000), watchParticle(-1), particleWatchStep(-1), PS(1e3),
    Pot(new double*[NumPot]), dPdR(new double*[NumPot]), potRangeScale(PS), writeSnapShot(false), mRotationChanged(false)
{
	//printf("Calculation::Calculation\n");
    double IntDist = 20.0, st;
    int n, x, y;
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
    if (nullptr != PotSs) for (n=0; n < NumPot; ++n) if (nullptr != PotSs[n].pot) setPotential(static_cast<PotRole>(n), PotSs[n]);
	PZS = int(round(MaxZ / Re));
	PYS = int(round(MaxY / Re)); 
	PXS = int(round(MaxX / Re));
    N = 2 * PXS * PZS;
	P = new Particle[N];
	D = new Particle*[N];
	G = new Particle***[XS];
    MAR = new MARStruct*[N];
	for (x=0; x < XS; x++)
	{
		G[x] = new Particle**[YS];
		for (y=0; y < YS; y++) G[x][y] = new Particle*[ZS];
	}
    for (n=0; n < N; n++) MAR[n] = new MARStruct[4];
    Pos = new Vector[N];
	Fixed = new bool[N];
	
	initialize();

    DebugLogFile = new QFile("DebugLog.txt");
    DebugLogFile->open(QIODevice::WriteOnly);
    DebugLog = new QTextStream(DebugLogFile);
    *DebugLog << "It\t";
    for (n=0; n<N; ++n) *DebugLog << "E[" << n << "]    \tT    \tdeltaE\tnew E  \tM\t";
    *DebugLog << "\n";
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
	delete[] Fixed;
}

void Calculation::calcMAR()
{
    int n, p, mx, my, mz, lx, ly, lz, i1, i2;
    double r;
    Particle *PP1, *PP2;
    for (n=0; n<N; n++) for (p=0; p<4; p++) MAR[n][p].R = RM;
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
                for (n=0; (n<4 ? r > MAR[i1][n].R : false); n++) ;
                for (p=3; p>n; p--) MAR[i1][p] = MAR[i1][p-1];
                if (n<4)
                {
                    MAR[i1][n].R = r;
                    MAR[i1][n].index = i2;
                }
                for (n=0; (n<4 ? r > MAR[i2][n].R : false); n++) ;
                for (p=3; p>n; p--) MAR[i2][p] = MAR[i2][p-1];
                if (n<4)
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
                for (n=0; (n<4 ? r > MAR[i1][n].R : false); n++) ;
                for (p=3; p>n; p--) MAR[i1][p] = MAR[i1][p-1];
                if (n<4)
                {
                    MAR[i1][n].R = r;
                    MAR[i1][n].index = i2;
                }
                for (n=0; (n<4 ? r > MAR[i2][n].R : false); n++) ;
                for (p=3; p>n; p--) MAR[i2][p] = MAR[i2][p-1];
                if (n<4)
                {
                    MAR[i2][n].R = r;
                    MAR[i2][n].index = i1;
                }
            }
        }
    }
}

Calculation::Result Calculation::geta(Vector* t0, Vector *a)
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
                for (PP2 = PP1->next; PP2 != 0; PP2 = PP2->next) if (getU(PP1, PP2, U, t0, temporaryPos, a) == Error) return Error;
			//printf("l1\n");
			for (lz = mz; lz < ZS && lz <= mz + GridSizeDiv; lz++)
				for (ly = (lz > mz ? ((ly = my - GridSizeDiv) >= 0 ? ly : 0) : my); 
						ly < YS && ly <= my + GridSizeDiv; ly++)
					for (lx = (ly > my || lz > mz ? ((lx = mx - GridSizeDiv) >= 0 ? lx : 0) : mx + 1);
							lx < XS && lx <= mx + GridSizeDiv; lx++)
						for (PP2 = G[lx][ly][lz]; PP2 != 0; PP2 = PP2->next)
                            for (PP1 = G[mx][my][mz]; PP1 != 0; PP1 = PP1->next) if (getU(PP1, PP2, U, t0, temporaryPos, a) == Error) return Error;
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

Calculation::Result Calculation::getU(const Particle * const P1, const Particle * const P2, double &U, const Vector* const t0, Positions pos, Vector* a) const
{
    double r, amp(0.0);
    Vector d, b;
    bool calcA = (NULL != a);
    int i1 = P1 - P, i2 = P2 - P, p, bi1=N, bi2=N;
    //printf("i1=%d, i2=%d\n", i1, i2);
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
    if (p < 0 || p >= NPot) return Success;
    for (int n=0; n<4; ++n)
    {
        if (P1->bound[n] == P2) bi1=n;
        if (P2->bound[n] == P1) bi2=n;
    }
    //printf("p=%d, r=%f, Rm=%f, PS=%f\n", p, r, Rm, PS);
    if (bi1 <= 1 && bi2 <= 1)
    {
        if (Pot[ClosestTwo][p] > UaMax)
            return Error;
        if (calcA) amp = dPdR[ClosestTwo][p] / r;
        U += Pot[ClosestTwo][p];
    }
    else if (bi1 <= 3 && bi2 <= 3)
    {
        if (Pot[NextTwo][p] > UaMax)
            return Error;
        if (calcA) amp = dPdR[NextTwo][p] / r;
        U += Pot[NextTwo][p];
    }
    else
    {
        bool SecondOrderBound(false);
        for (int n=0; !SecondOrderBound && n<4; ++n) for (int m=0; m<4; ++m) if (P1->bound[n] == P2->bound[m])
        {
            SecondOrderBound = true;
            break;
        }
        if (SecondOrderBound)
        {
            if (Pot[SecondOrder][p] > UaMax)
                return Error;
            if (calcA) amp = dPdR[SecondOrder][p] / r;
            U += Pot[SecondOrder][p];
        }
        else
        {
            if (Pot[Remaining][p] > UaMax)
                return Error;
            if (calcA) amp = dPdR[Remaining][p] / r;
            U += Pot[Remaining][p];
        }
    }
    if (abs(amp * r) > UaMax)
    {
        //printf("i1=%d, i2=%d, a=%f: Stepsize gets reduced!\n", i1, i2, amp*r);
        return Error;
    }
    if (calcA)
    {
        a[i1] += (b = amp * d);
        a[i2] -= b;
        if (particleWatchStep >= 0)
        {
            if (watchParticle == i1) ParticleWatchPoint->set(particleWatchStep, i2, amp * d);
            else if (watchParticle == i2) ParticleWatchPoint->set(particleWatchStep, i1, -amp * d);
        }
    }
    return Success;
}

void Calculation::correctLocalE()
{
    int mx, my, mz, n;
    Particle *PP1;
    double currSumE = 0.0, currT(0.0), EStart[N], TStart[N];
    int M[N];
    for (mz = 0; mz < ZS; ++mz) for (my = 0; my < YS; ++my) for (mx = 0; mx < XS; ++mx)
    {
        for (PP1 = G[mx][my][mz]; PP1 != 0; PP1 = PP1->next)
        {
            updateDelta(PP1->T, PP1->deltaT, 0.5 * (PP1->v.lengthSquared()));
            updateDelta(PP1->U, PP1->deltaU, getE(PP1, PP1->R, false));
            updateDelta(PP1->E, PP1->deltaE, PP1->T + PP1->U);
            currT += PP1->T;
            currSumE += PP1->E;
            /*if (isnan(currSumE))
			{
				printf("After calc delta and sum energies: Particle %ld is nan!\n", PP1 - P);
				Run = false;
            }*/
        }
    }
    emit EnergiesChanged(currT, currSumE);
    for (n=0; n<N; ++n)
    {
        EStart[n] = P[n].E;
        TStart[n] = P[n].T;
        M[n] = 0;
    }
    printf("current temporary energy = %g\n", currSumE);
    if (currSumE > Energy)
    {
        int *Sort = utils::heapSort(DeltaESortFunctor(P), N), EOrder[N];
        for (n=0; n<N; ++n) EOrder[Sort[n]] = n;
        delete[] Sort;
        for (n=0; n<N && currSumE > Energy; ++n)
        {
            Particle* curPar = P + EOrder[n];
            if (curPar->deltaE <= 0.0) break;
            const double EDelta = (curPar->deltaE > currSumE - Energy ? currSumE - Energy : curPar->deltaE);
            if (EDelta >= curPar->T)
            {
                if (curPar->U > 0.0)
                {
                    const double TMin = 25e4;
                    if (curPar->T > TMin)
                    {
                        curPar->v *= sqrt(TMin / curPar->T);
                        currSumE -= curPar->T + TMin;
                        M[EOrder[n]] = 4;
                    }
                    else M[EOrder[n]] = 5;
                    curPar->E -= curPar->deltaE;
                    curPar->U -= curPar->deltaU;
                    curPar->T -= curPar->deltaT;
                }
                else
                {
                    curPar->v.clear();
                    currSumE -= curPar->T;
                    if (EDelta > curPar->T)
                    {
                        curPar->E -= curPar->deltaE;
                        curPar->U -= curPar->deltaU;
                        curPar->T -= curPar->deltaT;
                        M[EOrder[n]] = 1;
                    }
                    else
                    {
                        curPar->E -= curPar->T;
                        curPar->T = 0.0;
                        M[EOrder[n]] = 2;
                    }
                }

            }
            else
            {
                curPar->v *= sqrt((curPar->T - EDelta) / curPar->T);
                curPar->T -= EDelta;
                curPar->E -= EDelta;
                currSumE -= EDelta;
                M[EOrder[n]] = 3;
                /*if (isnan(vF))
                {
                    printf("After T reduction: Particle %d is nan!\n", n);
                    Run = false;
                }*/
            }
        }
    }
    for (n=0; n<N; ++n)
        *DebugLog << "\t" << QString::number(EStart[n], 'g') << '\t' << QString::number(TStart[n], 'g')
                  << '\t' << QString::number(P[n].deltaE, 'g') << '\t' << QString::number(P[n].E, 'g') << '\t' << M[n];
    *DebugLog << '\n';
}

void Calculation::updateDelta(double &toUpdate, double &delta, const double newValue)
{
    delta = newValue - toUpdate;
    toUpdate = newValue;
}

double Calculation::getE(const Particle * const cP, const Vector &R, const bool useLastPos) const
{
    const Particle* P2;
    int lx, ly, lz;
    double E = 0.0;
    for (lz = ((lz = cP->zp - GridSizeDiv) >= 0 ? lz : 0); lz < ZS && lz <= cP->zp + GridSizeDiv; lz++)
        for (ly = ((ly = cP->yp - GridSizeDiv) >= 0 ? ly : 0); ly < YS && ly <= cP->yp + GridSizeDiv; ly++)
            for (lx = ((lx = cP->xp - GridSizeDiv) >= 0 ? lx : 0); lx < XS && lx <= cP->xp + GridSizeDiv; lx++)
                for (P2 = G[lx][ly][lz]; P2 != 0; P2 = P2->next) if (P2 != cP) if (getU(cP, P2, E, &R, (useLastPos ? lastPos : currentPos), nullptr) == Error) return Error_Double;
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
                for (PP2 = PP1->next; PP2 != 0; PP2 = PP2->next) if (getU(PP1, PP2, U, nullptr, particles, nullptr) == Error) return Error_Double;
            //printf("l1\n");
            for (lz = mz; lz < ZS && lz <= mz + GridSizeDiv; lz++)
                for (ly = (lz > mz ? ((ly = my - GridSizeDiv) >= 0 ? ly : 0) : my);
                        ly < YS && ly <= my + GridSizeDiv; ly++)
                    for (lx = (ly > my || lz > mz ? ((lx = mx - GridSizeDiv) >= 0 ? lx : 0) : mx + 1);
                            lx < XS && lx <= mx + GridSizeDiv; lx++)
                        for (PP2 = G[lx][ly][lz]; PP2 != 0; PP2 = PP2->next)
                            for (PP1 = G[mx][my][mz]; PP1 != 0; PP1 = PP1->next) if (getU(PP1, PP2, U, nullptr, particles, nullptr) == Error) return Error_Double;
        }
    }
    return U;
}

void Calculation::getScales(double& rScF, int& rMaxZ)
{
	rScF = ScF;
	rMaxZ = MaxZ;
}

void Calculation::getSize(int& width, int& height)
{
	width = int(MaxX * ScF) + 10;
	height = int(MaxY * ScF) + 40;
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
	int n, x, y, z;
	double rx, rz, rys = sqrt(0.5) * MaxY / double(PYS), rxs = MaxX / double(PXS), rzs = MaxZ / double(PZS);
	double y1 = 0.5 * MaxY - rys;
	double y2 = 0.5 * MaxY + rys;
    Vector F(double(XS) / MaxX, double(YS) / MaxY, double(ZS) / MaxZ * 1.001);
	for (x=0; x < XS; x++) for (y=0; y < YS; y++) for (z=0; z < ZS; z++) G[x][y][z] = 0;
	for (n=z=0, rz = 0.5 * rzs; z < PZS; z++, rz += rzs) 
	{
		for (x=0, rx = 0.5 * rxs; x < PXS; x++, rx += rxs)
		{
            Vector r1(rx, y1, rz), r2(rx, y2, rz);
            initializeParticle(P[n++], x, z, r1, F);
            initializeParticle(P[n++], x, z, r2, F);
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
            P[n].bound[P[n].NB++] = P+z;
            P[z].bound[P[z].NB++] = P+n;
        }
    }
    for (n=0; n<N; ++n) for (x=0; x<4 && MAR[n][x].R != RM; ++x)
    {
        z = MAR[n][x].index;
        if (z <= n) continue;
        for (y=0; y<4 && MAR[z][y].R != RM && MAR[z][y].index != n; ++y) ;
        if (y<4 && MAR[z][y].index == n && P[n].bound[0] != P+z && P[n].bound[1] != P+z)
        {
            P[n].bound[P[n].NB++] = P+z;
            P[z].bound[P[z].NB++] = P+n;
        }
    }
}

void Calculation::initializeParticle(Particle &cP, const int x, const int z, const Vector &iR, const Vector &Fact) const
{
    int n = &cP - P, b;
    double R;
    cP.lR = cP.R = iR;
    cP.xp = ((b = int(Fact.X() * iR.X())) >= 0 ? (b < XS ? b : XS - 1) : 0);
    cP.yp = ((b = int(Fact.Y() * iR.Y())) >= 0 ? (b < YS ? b : YS - 1) : 0);
    cP.zp = ((b = int(Fact.Z() * iR.Z())) >= 0 ? (b < ZS ? b : ZS - 1) : 0);
	cP.prev = 0;
	cP.next = G[cP.xp][cP.yp][cP.zp];
    cP.NB = 0;
    memset(cP.bound, 0, 4u * sizeof(Particle*));
	G[cP.xp][cP.yp][cP.zp] = D[n] = &cP;
	if (cP.next != 0) cP.next->prev = D[n];
	if (z == 0 || x == 0 || z == PZS - 1 || x == PXS - 1) 
	{
		Fixed[n] = true;
        Vector delta = P[x].R - Vector(0.5 * MaxX, 0.0, 0.5 * MaxZ);
        R = 1.0 / delta.length();
        P[x].v = delta * Vector(R, 0.0, R);
	}
	else Fixed[n] = false;
}

void Calculation::setLayerDistance(double newDistance)
{
    const double cDist(P[1].R.Y() - P[0].R.Y());
    const Vector diff(0.0, 0.5 * (newDistance - cDist), 0.0);
    for (int n=0; n < N-1; n+=2)
    {
        P[n].R -= diff;
        P[n+1].R += diff;
    }
}

void Calculation::move()
{
	if (Move) Move = false;
	else Move = true;
}

void Calculation::run()
{
    // Contains the rk4 algorithm from Numerical Recipes, Third Edition
    for (int n=0; n < NumPot; ++n)
    {
        if (Pot[n] == nullptr)
        {
            qCritical() << "Pot[" << n << "] == NULL!";
            return;
        }
        if (dPdR[n] == nullptr)
        {
            qCritical() << "dPdR[" << n << "] == NULL!";
            return;
        }
        if (!potentialOK[n])
        {
            qCritical() << "Pot[" << n << "] is not OK!";
            return;
        }
    }
    int n, i, x, y, z; // m;
    bool isNotFirstIt(false);
    Vector *a = new Vector[N], *dm = new Vector[N], *dvm = new Vector[N], *dt = new Vector[N], *t0 = new Vector[N];
    Vector *dvt = new Vector[N];
	Particle *PB;
	for (i=0, Run = true; Run; i++)
	{
		if (i==30178)
        {
            printf("Break!");
        }

        rk4(t0, dvt, a, dt, dm, dvm, h);

        if (writeSnapShot)
        {
            WriteSnapshot();
            writeSnapShot = false;
        }
        if (watchParticle >= 0)
            ParticleWatchPoint->setSum(particleWatchStep, dvt[watchParticle] - Vector(ParticleWatchPoint->getSumX(1), ParticleWatchPoint->getSumY(1), ParticleWatchPoint->getSumZ(1)));
        if (i==147)
        {
            printf("Break!");
        }
        printf("iteration=%d, ", i);
        updateBindings();
        *DebugLog << i;
        if (isNotFirstIt) correctLocalE();
        else isNotFirstIt = true;
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

void Calculation::updateBlock(int n)
{
    static const Vector F(double(XS) / MaxX, double(YS) / MaxY, double(ZS) / MaxZ);
    int x, y, z;
    x = ((x = int(F.X() * P[n].R.X())) >= 0 ? (x < XS ? x : XS - 1) : 0);
    y = ((y = int(F.Y() * P[n].R.Y())) >= 0 ? (y < YS ? y : YS - 1) : 0);
    z = ((z = int(F.Z() * P[n].R.Z())) >= 0 ? (z < ZS ? z : ZS - 1) : 0);
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

void Calculation::rk4(Vector *t0, Vector *dvt, Vector *a, Vector *dt, Vector* dm, Vector* dvm, const double lh)
{

    double hh = 0.5 * lh, h6 = lh / 6.0, R, dX, dZ, ZMid = 0.5 * MaxZ;
    double XMid = 0.5 * MaxX;
    int n;
    Result result = Success;
    for (int i=1; i==1; ++i)
    {
        for (n=0; n<N; n++)
        {
            updateBlock(n);
            t0[n] = P[n].R;
        }
        if (Move)
        {
            for (n=0; n<N; n++) if (Fixed[n])
            {
                dX = XMid - P[n].R.X();
                dZ = ZMid - P[n].R.Z();
                R = lh * Speed / (dX * dX + dZ * dZ);
                if (P[n].R.Y() > YMid) P[n].R += Vector(dX * R, 0.0, dZ * R);
                else P[n].R -= Vector(dX * R, 0.0, dZ * R);
            }
        }
        for (n=0, U = T = 0.0; n<N; n++) T += P[n].v.lengthSquared();
        if (watchParticle >= 0) particleWatchStep = 0;
        result = geta(t0, a);
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
        for (n=0; n<N; n++) if (!Fixed[n])
        {
            t0[n] = P[n].R + hh * P[n].v;
            dt[n] = P[n].v + hh * a[n];
        }
        for (n=0, U = T = 0.0; n<N; n++) T += dt[n].lengthSquared();
        if (watchParticle >= 0) ParticleWatchPoint->setSum(particleWatchStep++, a[watchParticle]);
        result = geta(t0, dvt);
        if (result == Error) break;
        for (n=0; n<N; n++) if (!Fixed[n])
        {
            P[n].aa = dvt[n];
            t0[n] = P[n].R + hh * dt[n];
            dm[n] = P[n].v + hh * dvt[n];
        }
        for (n=0, U = T = 0.0; n<N; n++) T += dm[n].lengthSquared();
        if (watchParticle >= 0) ParticleWatchPoint->setSum(particleWatchStep++, dvt[watchParticle]);
        result = geta(t0, dvm);
        if (result == Error) break;
        for (n=0; n<N; n++) if (!Fixed[n])
        {
            P[n].aa += 2.0 * dvm[n];
            t0[n] = P[n].R + lh * dm[n];
            dm[n] += dt[n];
            dt[n] = P[n].v + lh * dvm[n];
            dvm[n] += dvt[n];
        }
        for (n=0, U = T = 0.0; n<N; n++) T += dt[n].lengthSquared();
        if (watchParticle >= 0) ParticleWatchPoint->setSum(particleWatchStep++, dvm[watchParticle]);
        result = geta(t0, dvt);
        if (result == Error) break;
        for (n=0; n<N; n++) if (!Fixed[n])
        {
            P[n].lR = P[n].R;
            P[n].lv = P[n].v;
            P[n].aa += dvt[n];
            P[n].R += h6 * (P[n].v + dt[n] + 2.0 * dm[n]);
            P[n].v += h6 * (a[n] + dvt[n] + 2.0 * dvm[n]);
            if ((P[n].R.X() < 0.0 && P[n].v.X() < 0.0) || (P[n].R.X() > MaxX && P[n].v.X() > 0.0)) P[n].v *= Vector(-1.0, 0.0, 0.0);
            if ((P[n].R.Y() < 0.0 && P[n].v.Y() < 0.0) || (P[n].R.Y() > MaxY && P[n].v.Y() > 0.0)) P[n].v *= Vector(0.0, -1.0, 0.0);
            if ((P[n].R.Z() < 0.0 && P[n].v.Z() < 0.0) || (P[n].R.Z() > MaxZ && P[n].v.Z() > 0.0)) P[n].v *= Vector(0.0, 0.0, -1.0);
            /*if (isnan(P[n].X) || isnan(P[n].Y) || isnan(P[n].Z) || isnan(P[n].vX) || isnan(P[n].vY) || isnan(P[n].vZ))
            {
                printf("After calculation of new position and v: Particel %d is nan!\n", n);
                Run = false;
            }*/
        }
    }
    if (result == Error)
    {
        const double nh = 0.5 * lh;
        if (nh < 1e-10)
        {
            printf("Break!");
        }
        else
        {
            rk4(t0, dvt, a, dt, dm, dvm, nh);
            rk4(t0, dvt, a, dt, dm, dvm, nh);
        }
    }
}

bool Calculation::updateBindings()
{
    bool rValue = false;
    double randF = static_cast<double>(static_cast<int>(Particle::NBound) * N) / RAND_MAX;
    for (int n=0; n < static_cast<int>(Particle::NBound)*N/2; ++n)
    {
        int random(static_cast<int>(static_cast<double>(rand()) * randF));
        int i0 = random % N, i1(random % static_cast<int>(Particle::NBound));
        if (i0 == N || i1 == static_cast<int>(Particle::NBound) || nullptr == P[i0].bound[i1]) continue;
        std::map<double,int> map1, map2;
        for (int i2 = 0; i2 < static_cast<int>(Particle::NBound); ++i2)
        {
            if (i1 != i2 && nullptr != P[i0].bound[i2] && isNotBound(P[i0].bound[i1], P[i0].bound[i2]))
                map1.insert(std::make_pair(dist(P + i0, P[i0].bound[i2]) - dist(P[i0].bound[i1], P[i0].bound[i2]), i2));
            if (P + i0 != P[i0].bound[i1]->bound[i2] && nullptr != P[i0].bound[i1]->bound[i2] && isNotBound(P + i0, P[i0].bound[i1]->bound[i2]))
                map2.insert(std::make_pair(dist(P[i0].bound[i1], P[i0].bound[i1]->bound[i2]) - dist(P + i0, P[i0].bound[i1]->bound[i2]), i2));
        }
        for (std::map<double,int>::const_reverse_iterator it1 = map1.rbegin(), it2 = map2.rbegin(); it1 != map1.rend() && it2 != map2.rend() && it1->first + it2->first > 0.0; ++it1, ++it2)
        {
            rValue = true;
            for (int i4 = 0; i4 < static_cast<int>(Particle::NBound); ++i4)
            {
                if (P[i0].bound[it1->second]->bound[i4] == P + i0) P[i0].bound[it1->second]->bound[i4] = P[i0].bound[i1];
                if (P[i0].bound[i1]->bound[it2->second]->bound[i4] == P[i0].bound[i1]) P[i0].bound[i1]->bound[it2->second]->bound[i4] = P + i0;
            }
            std::swap(P[i0].bound[it1->second], P[i0].bound[i1]->bound[it2->second]);
        }
    }
    return rValue;
}

double Calculation::dist(const Particle * const P1, const Particle * const P2)
{
    return (P1->R - P2->R).length();
}

bool Calculation::isNotBound(const Particle *const P1, const Particle *const P2)
{
    for (int i=0; i<4; ++i) if (P1->bound[i] == P2) return false;
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

double Calculation::setKineticEnergy(const double newT)
{
	int n;
    double nE = getPotentialEnergy() + newT, ParE, ParV;
	double VC, A1, A2, RD = M_PI / RAND_MAX, EnDiff = 2.0 * (nE - Energy) / double(N);
    for (n=0; n < NumPot; ++n) if (Pot[n] == nullptr || dPdR[n] == nullptr) return -1.0;
    qInfo() << "SetKineticEnergy: newT=" << newT << " nE=" << nE << " Energy=" << Energy;
    if (nE > Energy)
	{
		Energy = nE;
		VC = sqrt(EnDiff);
		for (n=0; n<N; n++)
		{
            if (P[n].v == Vector(0.0, 0.0, 0.0))
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
	else
	{
        int *Sort = utils::heapSort(VSortFunctor(P), N), EOrder[N];
        for (n=0; n<N; ++n) EOrder[Sort[n]] = n;
        delete[] Sort;
        for (n=0; Energy > nE && n<N; ++n)
        {
            double T = 0.5 * P[EOrder[n]].v.lengthSquared();
            if (T <= Energy - nE)
            {
                P[EOrder[n]].v.clear();
                Energy -= T;
            }
            else
            {
                const double vF = sqrt((T - Energy + nE) / T);
                P[EOrder[n]].v *= vF;
                Energy = nE;
            }
        }
        Energy = nE;
	}
	E = Energy;
	return Energy;
}

bool Calculation::setPotential(const PotRole role, PotStruct &PotS)
{
    if (Pot[role] != nullptr) delete[] Pot[role];
    if (dPdR[role] != nullptr) delete[] dPdR[role];
    const double dRScale = 1.0 / PotS.RZoom, Rmin = Rm * dRScale, Rmax = RM * dRScale, devF = PotS.VZoom * dRScale;
    Pot[role] = PotS.pot->getPoints(Rmin, Rmax, NPot);
    dPdR[role] = PotS.pot->get_dVdR(Rmin, Rmax, NPot);
    if (PotS.RZoom != 1.0 || PotS.VZoom != 1.0) for (int n=0; n < NPot; ++n)
    {
        Pot[role][n] *= PotS.VZoom;
        dPdR[role][n] *= devF;
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
        if (Pot[ClosestTwo]  != nullptr) for (int n=0; n < NPot && Pot[ClosestTwo][n] > 0.0; ++n) if (Pot[role][n] < Pot[ClosestTwo][n])
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
	h = nh;
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
        if (updateBindings()) data.addBoundChange(n);
        const Particle* P2;
        double SecondOrderSum = 0.0, UnboundSum = 0.0;
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
                if (currentParticle.bound[m] == P2) bi1=m;
                if (P2->bound[m] == &currentParticle) bi2=m;
            }
            if (bi1 <= 1 && bi2 <= 1)
            {
                if (bi1 == 0) data.SetFirstBound(n, Pot[ClosestTwo][p]);
                else data.SetSecondBound(n, Pot[ClosestTwo][p]);
            }
            else if (bi1 <= 3 && bi2 <= 3)
            {
                if (bi1 == 2) data.SetThirdBound(n, Pot[NextTwo][p]);
                else data.SetFourthBound(n, Pot[NextTwo][p]);
            }
            else
            {
                bool SecondOrderBound(false);
                for (int n=0; !SecondOrderBound && n<4; ++n) for (int m=0; m<4; ++m) if (currentParticle.bound[n] == P2->bound[m])
                {
                    SecondOrderBound = true;
                    break;
                }
                if (SecondOrderBound) SecondOrderSum += Pot[SecondOrder][p];
                else UnboundSum += Pot[Remaining][p];
            }
        }
        data.SetSecondOrderBound(n, SecondOrderSum);
        data.SetUnbound(n, UnboundSum);
    }
    currentParticle.R = particlePos;
    updateBlock(particleIndex);
    updateBindings();
}
