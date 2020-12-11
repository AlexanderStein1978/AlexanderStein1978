#define _USE_MATH_DEFINES

#include "Calculation.h"
#include "heapsort.h"
#include "deltaesortfunctor.h"
#include "utils.h"
#include "potstruct.h"
#include "potential.h"
#include "watchpoint.h"

#include <QGridLayout>
#include <QFile>
#include <QLabel>
#include <QTextStream>
#include <QPainter>
#include <QPaintEvent>

#include <math.h>
#include <cstdlib>


const double UaMax = 1e6;


Calculation::Calculation(PotStruct* PotSs, QObject* parent): QThread(parent), NPot(30000), watchParticle(-1), particleWatchStep(-1), PS(1e3),
    Pot(new double*[NumPot]), dPdR(new double*[NumPot]), potRangeScale(PS), writeSnapShot(false)
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
        if (nullptr != PotSs && nullptr != PotSs[n].pot) setPotential(static_cast<PotRole>(n), PotSs[n]);
    }
	
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
	XP = new double[N];
	YP = new double[N];
	ZP = new double[N];
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
	delete[] XP;
	delete[] YP;
	delete[] ZP;
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
    double r, dx, dy, dz;
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
                dx = PP2->X - PP1->X;
                dy = PP2->Y - PP1->Y;
                dz = PP2->Z - PP1->Z;
                r = sqrt(dx * dx + dy * dy + dz * dz);
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
                dx = PP2->X - PP1->X;
                dy = PP2->Y - PP1->Y;
                dz = PP2->Z - PP1->Z;
                r = sqrt(dx * dx + dy * dy + dz * dz);
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

void Calculation::geta(double *tx, double *ty, double *tz, double *ax, double *ay, double *az)
{
	//printf("geta\n");
    int mx, my, mz, lx, ly, lz;// i1, i2, p, n;
    //double r, dx, dy, dz;
	Particle *PP1, *PP2;
	for (mx = 0; mx < N; mx++) ax[mx] = ay[mx] = az[mx] = 0.0;
	for (mz = 0; mz < ZS; mz++) for (my = 0; my < YS; my++) for (mx = 0; mx < XS; mx++)
	{
		if (G[mx][my][mz] != 0)
		{
			//printf("l0\n");
			for (PP1 = G[mx][my][mz]; PP1->next != 0; PP1 = PP1->next)
                for (PP2 = PP1->next; PP2 != 0; PP2 = PP2->next) getU(PP1, PP2, U, tx, ty, tz, temporaryPos, ax, ay, az);
			//printf("l1\n");
			for (lz = mz; lz < ZS && lz <= mz + GridSizeDiv; lz++)
				for (ly = (lz > mz ? ((ly = my - GridSizeDiv) >= 0 ? ly : 0) : my); 
						ly < YS && ly <= my + GridSizeDiv; ly++)
					for (lx = (ly > my || lz > mz ? ((lx = mx - GridSizeDiv) >= 0 ? lx : 0) : mx + 1);
							lx < XS && lx <= mx + GridSizeDiv; lx++)
						for (PP2 = G[lx][ly][lz]; PP2 != 0; PP2 = PP2->next)
                            for (PP1 = G[mx][my][mz]; PP1 != 0; PP1 = PP1->next) getU(PP1, PP2, U, tx, ty, tz, temporaryPos, ax, ay, az);
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
}

void Calculation::getU(const Particle * const P1, const Particle * const P2, double &U, const double* tx, const double * const ty,
                       const double* const tz, Positions pos, double* ax, double* ay, double* az) const
{
    double dx, dy, dz, r, a, b;
    bool calcA = (NULL != ax && NULL != ay && NULL != az);
    int i1 = P1 - P, i2 = P2 - P, p, bi1=N, bi2=N;
    //printf("i1=%d, i2=%d\n", i1, i2);
    switch (pos)
    {
    case temporaryPos:
        dx = tx[i2] - tx[i1];
        dy = ty[i2] - ty[i1];
        dz = tz[i2] - tz[i1];
        break;
    case lastPos:
        dx = P2->lX - *tx;
        dy = P2->lY - *ty;
        dz = P2->lZ - *tz;
        break;
    case currentPos:
        dx = P2->X - *tx;
        dy = P2->Y - *ty;
        dz = P2->Z - *tz;
        break;
    case particles:
        dx = P2->X - P1->X;
        dy = P2->Y - P1->Y;
        dz = P2->Z - P1->Z;
        break;
    }
    //printf("tx1=%f, tx2=%f, ty1=%f, ty2=%f, tz1=%f, tz2=%f\n",
        //   tx[i1], tx[i2], ty[i1], ty[i2], tz[i1], tz[i2]);
    r = sqrt(dx * dx + dy * dy + dz * dz);
    p = int((r - Rm) * potRangeScale);
    if (p < 0 || p >= NPot) return;
    for (int n=0; n<4; ++n)
    {
        if (P1->bound[n] == P2) bi1=n;
        if (P2->bound[n] == P1) bi2=n;
    }
    //printf("p=%d, r=%f, Rm=%f, PS=%f\n", p, r, Rm, PS);
    if (bi1 <= 1 && bi2 <= 1)
    {
        if (calcA) a = dPdR[ClosestTwo][p] / r;
        U += Pot[ClosestTwo][p];
    }
    else if (bi1 <= 3 && bi2 <= 3)
    {
        if (calcA) a = dPdR[NextTwo][p] / r;
        U += Pot[NextTwo][p];
    }
    else
    {
        if (calcA) a = dPdR[Remaining][p] / r;
        U += Pot[Remaining][p];
    }
    if (abs(a) > UaMax)
    {
        printf("i1=%d, i2=%d, a=%f gets reduced to %f\n", i1, i2, a, UaMax);
        a = UaMax;
    }
    if (calcA)
    {
        ax[i1] += (b = a * dx);
        ax[i2] -= b;
        ay[i1] += (b = a * dy);
        ay[i2] -= b;
        az[i1] += (b = a * dz);
        az[i2] -= b;
        if (particleWatchStep >= 0)
        {
            if (watchParticle == i1) ParticleWatchPoint->set(particleWatchStep, i2, a * dx, a * dy, a * dz);
            else if (watchParticle == i2) ParticleWatchPoint->set(particleWatchStep, i1, -a * dx, -a * dy, -a * dz);
        }
    }
}

void Calculation::correctLocalE()
{
    int mx, my, mz, n;
    Particle *PP1;
    double currSumE = 0.0, EStart[N], TStart[N];
    int M[N];
    for (mz = 0; mz < ZS; ++mz) for (my = 0; my < YS; ++my) for (mx = 0; mx < XS; ++mx)
    {
        for (PP1 = G[mx][my][mz]; PP1 != 0; PP1 = PP1->next)
        {
            updateDelta(PP1->T, PP1->deltaT, 0.5 * (PP1->vX * PP1->vX + PP1->vY * PP1->vY + PP1->vZ * PP1->vZ));
            updateDelta(PP1->U, PP1->deltaU, getE(PP1, PP1->X, PP1->Y, PP1->Z, false));
            updateDelta(PP1->E, PP1->deltaE, PP1->T + PP1->U);
            currSumE += PP1->E;
            /*if (isnan(currSumE))
			{
				printf("After calc delta and sum energies: Particle %ld is nan!\n", PP1 - P);
				Run = false;
            }*/
        }
    }
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
                        const double vF = sqrt(TMin / curPar->T);
                        curPar->vX *= vF;
                        curPar->vY *= vF;
                        curPar->vZ *= vF;
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
                    curPar->vX = curPar->vY = curPar->vZ = 0.0;
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
                const double vF = sqrt((curPar->T - EDelta) / curPar->T);
                curPar->vX *= vF;
                curPar->vY *= vF;
                curPar->vZ *= vF;
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

double Calculation::getE(const Particle * const cP, const double X, const double Y, const double Z, const bool useLastPos) const
{
    const Particle* P2;
    int lx, ly, lz;
    double E = 0.0;
    for (lz = ((lz = cP->zp - GridSizeDiv) >= 0 ? lz : 0); lz < ZS && lz <= cP->zp + GridSizeDiv; lz++)
        for (ly = ((ly = cP->yp - GridSizeDiv) >= 0 ? ly : 0); ly < YS && ly <= cP->yp + GridSizeDiv; ly++)
            for (lx = ((lx = cP->xp - GridSizeDiv) >= 0 ? lx : 0); lx < XS && lx <= cP->xp + GridSizeDiv; lx++)
                for (P2 = G[lx][ly][lz]; P2 != 0; P2 = P2->next) if (P2 != cP) getU(cP, P2, E, &X, &Y, &Z, (useLastPos ? lastPos : currentPos));
    return 0.5 * E;
}

double Calculation::getEnergy()
{
    int mx, my, mz, lx, ly, lz, n;
    double T, U;
	Particle *PP1, *PP2;
    for (n=0; n < NumPot; ++n) if (Pot[n] == nullptr || dPdR[n] == nullptr) return -1.0;
	for (mx = 0, T = U = 0.0; mx < N; mx++) 
		T += P[mx].vX * P[mx].vX + P[mx].vY * P[mx].vY + P[mx].vZ * P[mx].vZ;
	for (mz = 0; mz < ZS; mz++) for (my = 0; my < YS; my++) for (mx = 0; mx < XS; mx++)
	{
		if (G[mx][my][mz] != 0)
		{
			//printf("l0\n");
			for (PP1 = G[mx][my][mz]; PP1->next != 0; PP1 = PP1->next)
                for (PP2 = PP1->next; PP2 != 0; PP2 = PP2->next) getU(PP1, PP2, U, NULL, NULL, NULL, particles);
			//printf("l1\n");
			for (lz = mz; lz < ZS && lz <= mz + GridSizeDiv; lz++)
				for (ly = (lz > mz ? ((ly = my - GridSizeDiv) >= 0 ? ly : 0) : my); 
						ly < YS && ly <= my + GridSizeDiv; ly++)
					for (lx = (ly > my || lz > mz ? ((lx = mx - GridSizeDiv) >= 0 ? lx : 0) : mx + 1);
							lx < XS && lx <= mx + GridSizeDiv; lx++)
						for (PP2 = G[lx][ly][lz]; PP2 != 0; PP2 = PP2->next)
                            for (PP1 = G[mx][my][mz]; PP1 != 0; PP1 = PP1->next) getU(PP1, PP2, U, NULL, NULL, NULL, particles);
		}
	}
	return Energy = 0.5*T+U;
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
    double XF = double(XS) / MaxX, YF = double(YS) / MaxY, ZF = double(ZS) / MaxZ * 1.001;
	for (x=0; x < XS; x++) for (y=0; y < YS; y++) for (z=0; z < ZS; z++) G[x][y][z] = 0;
	for (n=z=0, rz = 0.5 * rzs; z < PZS; z++, rz += rzs) 
	{
		for (x=0, rx = 0.5 * rxs; x < PXS; x++, rx += rxs)
		{
            initializeParticle(P[n++], x, z, rx, y1, rz, XF, YF, ZF);
            initializeParticle(P[n++], x, z, rx, y2, rz, XF, YF, ZF);
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

void Calculation::initializeParticle(Particle &cP, const int x, const int z, const double X, const double Y, const double Z,
    const double XF, const double YF, const double ZF) const
{
    int n = &cP - P, b;
    double dX, dZ, R;
    cP.lvX = cP.lvY = cP.lvZ = cP.vX = cP.vY = cP.vZ = 0.0;
	cP.lX = cP.X = X;
	cP.lY = cP.Y = Y;
	cP.lZ = cP.Z = Z;
	cP.xp = ((b = int(XF * X)) >= 0 ? (b < XS ? b : XS - 1) : 0);
	cP.yp = ((b = int(YF * Y)) >= 0 ? (b < YS ? b : YS - 1) : 0);
	cP.zp = ((b = int(ZF * Z)) >= 0 ? (b < ZS ? b : ZS - 1) : 0);
	cP.prev = 0;
	cP.next = G[cP.xp][cP.yp][cP.zp];
    cP.NB = 0;
    memset(cP.bound, 0, 4u * sizeof(Particle*));
	G[cP.xp][cP.yp][cP.zp] = D[n] = &cP;
	if (cP.next != 0) cP.next->prev = D[n];
	if (z == 0 || x == 0 || z == PZS - 1 || x == PXS - 1) 
	{
		Fixed[n] = true;
		dX = P[x].X - 0.5 * MaxX;
		dZ = P[x].Z - 0.5 * MaxZ;
		R = 1.0 / sqrt(dX * dX + dZ * dZ);
		P[x].vX = dX * R;
		P[x].vZ = dZ * R;
	}
	else Fixed[n] = false;
}

void Calculation::move()
{
	if (Move) Move = false;
	else Move = true;
}

void Calculation::run()
{
    // Contains the rk4 algorithm from Numerical Recipes, Third Edition
    for (int n=0; n < NumPot; ++n) if (Pot[n] == nullptr || dPdR[n] == nullptr) return;
    int n, i, x, y, z; // m;
    bool isNotFirstIt(false);
	double hh = 0.5 * h, h6 = h / 6.0, *ax = new double[N], *ay = new double[N];
	double *az = new double[N], *dxm = new double[N], *dym = new double[N], *dzm = new double[N];
	double *dvxm = new double[N], *dvym = new double[N], *dvzm = new double[N];
	double *dxt = new double[N], *dyt = new double[N], *dzt = new double[N];
	double *dvxt = new double[N], *dvyt = new double[N], *dvzt = new double[N];
    double *xt = new double[N], *yt = new double[N], *zt = new double[N], R, dX, dZ;
    double XF = double(XS) / MaxX, YF = double(YS) / MaxY, ZF = double(ZS) / MaxZ, ZMid = 0.5 * MaxZ;
    double *Angle = new double[N], XMid = 0.5 * MaxX;
	Particle *PB;
	for (n=0; n<N; n++) if (Fixed[n]) Angle[n] = tan((XMid - P[n].X) / (ZMid - P[n].Z));
	for (i=0, Run = true; Run; i++)
	{
		if (i==30178)
        {
            printf("Break!");
        }
		for (n=0; n<N; n++)
		{
			x = ((x = int(XF * P[n].X)) >= 0 ? (x < XS ? x : XS - 1) : 0);
			y = ((y = int(YF * P[n].Y)) >= 0 ? (y < YS ? y : YS - 1) : 0);
			z = ((z = int(ZF * P[n].Z)) >= 0 ? (z < ZS ? z : ZS - 1) : 0);
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
			xt[n] = P[n].X;
			yt[n] = P[n].Y;
			zt[n] = P[n].Z;
		}
		if (Move)
		{
			for (n=0; n<N; n++) if (Fixed[n])
			{
				dX = XMid - P[n].X;
				dZ = ZMid - P[n].Z;
				R = h * Speed / (dX * dX + dZ * dZ);
				if (P[n].Y > YMid)
				{
					P[n].X += dX * R;
					P[n].Z += dZ * R;
				}
				else
				{
					P[n].X -= dX * R;
					P[n].Z -= dZ * R;
				}
			}
		}
		for (n=0, U = T = 0.0; n<N; n++) 
			T += P[n].vX * P[n].vX + P[n].vY * P[n].vY + P[n].vZ * P[n].vZ;
        if (watchParticle >= 0) particleWatchStep = 0;
		geta(xt, yt, zt, ax, ay, az);
		if (E == 0.0 || T == 0.0 || Move) E = T + U;
		for (n=0; n<N; n++) if (!Fixed[n])
		{
			xt[n] = P[n].X + hh * P[n].vX;
			yt[n] = P[n].Y + hh * P[n].vY;
			zt[n] = P[n].Z + hh * P[n].vZ;
			dxt[n] = P[n].vX + hh * ax[n];
			dyt[n] = P[n].vY + hh * ay[n];
			dzt[n] = P[n].vZ + hh * az[n];
		}
		for (n=0, U = T = 0.0; n<N; n++) 
			T += dxt[n] * dxt[n] + dyt[n] * dyt[n] + dzt[n] * dzt[n];
        if (watchParticle >= 0) ParticleWatchPoint->setSum(particleWatchStep++, ax[watchParticle], ay[watchParticle], az[watchParticle]);
		geta(xt, yt, zt, dvxt, dvyt, dvzt);
		for (n=0; n<N; n++) if (!Fixed[n])
		{
            P[n].aaX = dvxt[n];
            P[n].aaY = dvyt[n];
            P[n].aaZ = dvzt[n];
			xt[n] = P[n].X + hh * dxt[n];
			yt[n] = P[n].Y + hh * dyt[n];
			zt[n] = P[n].Z + hh * dzt[n];
			dxm[n] = P[n].vX + hh * dvxt[n];
			dym[n] = P[n].vY + hh * dvyt[n];
			dzm[n] = P[n].vZ + hh * dvzt[n];
		}
		for (n=0, U = T = 0.0; n<N; n++) 
			T += dxm[n] * dxm[n] + dym[n] * dym[n] + dzm[n] * dzm[n];
        if (watchParticle >= 0) ParticleWatchPoint->setSum(particleWatchStep++, dvxt[watchParticle], dvyt[watchParticle], dvzt[watchParticle]);
		geta(xt, yt, zt, dvxm, dvym, dvzm);
		for (n=0; n<N; n++) if (!Fixed[n])
		{
            P[n].aaX += 2.0 * dvxm[n];
            P[n].aaY += 2.0 * dvym[n];
            P[n].aaZ += 2.0 * dvzm[n];
			xt[n] = P[n].X + h * dxm[n];
			yt[n] = P[n].Y + h * dym[n];
			zt[n] = P[n].Z + h * dzm[n];
			dxm[n] += dxt[n];
			dym[n] += dyt[n];
			dzm[n] += dzt[n];
			dxt[n] = P[n].vX + h * dvxm[n];
			dyt[n] = P[n].vY + h * dvym[n];
			dzt[n] = P[n].vZ + h * dvzm[n];
			dvxm[n] += dvxt[n];
			dvym[n] += dvyt[n];
			dvzm[n] += dvzt[n];
		}
		for (n=0, U = T = 0.0; n<N; n++) 
			T += dxt[n] * dxt[n] + dyt[n] * dyt[n] + dzt[n] * dzt[n];
        if (watchParticle >= 0) ParticleWatchPoint->setSum(particleWatchStep++, dvxm[watchParticle], dvym[watchParticle], dvzm[watchParticle]);
		geta(xt, yt, zt, dvxt, dvyt, dvzt);
		for (n=0; n<N; n++) if (!Fixed[n])
		{
            P[n].lX = P[n].X;
            P[n].lY = P[n].Y;
            P[n].lZ = P[n].Z;
            P[n].lvX = P[n].vX;
            P[n].lvY = P[n].vY;
            P[n].lvZ = P[n].vZ;
            P[n].aaX += dvxt[n];
            P[n].aaY += dvyt[n];
            P[n].aaZ += dvzt[n];
			P[n].X += h6 * (P[n].vX + dxt[n] + 2.0 * dxm[n]);
			P[n].Y += h6 * (P[n].vY + dyt[n] + 2.0 * dym[n]);
			P[n].Z += h6 * (P[n].vZ + dzt[n] + 2.0 * dzm[n]);
			P[n].vX += h6 * (ax[n] + dvxt[n] + 2.0 * dvxm[n]);
			P[n].vY += h6 * (ay[n] + dvyt[n] + 2.0 * dvym[n]);
			P[n].vZ += h6 * (az[n] + dvzt[n] + 2.0 * dvzm[n]);
			if ((P[n].X < 0.0 && P[n].vX < 0.0) || (P[n].X > MaxX && P[n].vX > 0.0)) P[n].vX *= -1.0;
			if ((P[n].Y < 0.0 && P[n].vY < 0.0) || (P[n].Y > MaxY && P[n].vY > 0.0)) P[n].vY *= -1.0;
			if ((P[n].Z < 0.0 && P[n].vZ < 0.0) || (P[n].Z > MaxZ && P[n].vZ > 0.0)) P[n].vZ *= -1.0;
            /*if (isnan(P[n].X) || isnan(P[n].Y) || isnan(P[n].Z) || isnan(P[n].vX) || isnan(P[n].vY) || isnan(P[n].vZ))
			{
				printf("After calculation of new position and v: Particel %d is nan!\n", n);
				Run = false;
            }*/
		}
        if (writeSnapShot)
        {
            WriteSnapshot();
            writeSnapShot = false;
        }
        if (watchParticle >= 0)
            ParticleWatchPoint->setSum(particleWatchStep, dvxt[watchParticle] - ParticleWatchPoint->getSumX(1),
                                       dvyt[watchParticle] - ParticleWatchPoint->getSumY(1), dvzt[watchParticle] - ParticleWatchPoint->getSumZ(1));
        if (i==147)
        {
            printf("Break!");
        }
        printf("iteration=%d, ", i);
        updateBindings();
        *DebugLog << i;
        if (isNotFirstIt) correctLocalE();
        else isNotFirstIt = true;
		for (PB = P; PB != 0; ) for (n=1, PB = 0; n<N; n++) if (D[n]->Z < D[n-1]->Z)
		{
			PB = D[n];
			D[n] = D[n-1];
			D[n-1] = PB;
		}
		mutex.lock();
		if (rotated)
		{
			for (n=x=0, z=N; n < N; n++) 
			{
				y = D[n] - P;
				if (y < N && y >= 0)
				{
					XP[x] = D[n]->X;
					YP[x] = D[n]->Z;
					ZP[x++] = D[n]->Y;
				}
				else
				{
					XP[z] = D[n]->X;
					YP[z] = D[n]->Z;
					ZP[z++] = D[n]->Y;
				}
			}
		}
		else for (n=0; n < N; n++)
		{
			XP[n] = D[n]->X;
			YP[n] = D[n]->Y;
			ZP[n] = D[n]->Z;
		}
		mutex.unlock();
		emit PictureChanged(XP, YP, ZP, N);
        if (watchParticle >= 0)
        {
            particleWatchStep = -1;
            break;
        }
	}
	delete[] ax;
	delete[] ay;
	delete[] az;
	delete[] dxm;
	delete[] dym;
	delete[] dzm;
	delete[] dvxm;
	delete[] dvym;
	delete[] dvzm;
	delete[] dxt;
	delete[] dyt;
	delete[] dzt;
	delete[] dvxt;
	delete[] dvyt;
	delete[] dvzt;
	delete[] xt;
	delete[] yt;
	delete[] zt;
	delete[] Angle;
}

void Calculation::updateBindings()
{
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
            for (int i4 = 0; i4 < static_cast<int>(Particle::NBound); ++i4)
            {
                if (P[i0].bound[it1->second]->bound[i4] == P + i0) P[i0].bound[it1->second]->bound[i4] = P[i0].bound[i1];
                if (P[i0].bound[i1]->bound[it2->second]->bound[i4] == P[i0].bound[i1]) P[i0].bound[i1]->bound[it2->second]->bound[i4] = P + i0;
            }
            std::swap(P[i0].bound[it1->second], P[i0].bound[i1]->bound[it2->second]);
        }
    }
}

double Calculation::dist(const Particle * const P1, const Particle * const P2)
{
    double dx = P1->X - P2->X;
    double dy = P1->Y - P2->Y;
    double dz = P1->Z - P2->Z;
    return sqrt(dx * dx + dy * dy + dz * dz);
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
}

double Calculation::setEnergy(double newE)
{
	int n;
	double nE = newE, ParE, ParV, ERem = 0.0;
	double VC, A1, A2, RD = M_PI / RAND_MAX, EnDiff = 2.0 * (nE - Energy) / double(N);
    for (n=0; n < NumPot; ++n) if (Pot[n] == nullptr || dPdR[n] == nullptr) return -1.0;
	if (nE > Energy)
	{
		Energy = nE;
		VC = sqrt(EnDiff);
		for (n=0; n<N; n++)
		{
			if (P[n].vX == 0.0 && P[n].vY == 0.0 && P[n].vZ == 0.0)
			{
				A1 = 2.0 * RD * rand();
				A2 = RD * rand();
				P[n].vX = sin(A1) * sin(A2) * VC;
				P[n].vY = sin(A1) * cos(A2) * VC;
				P[n].vZ = cos(A1) * VC;
			}
			else
			{
				ParE = P[n].vX * P[n].vX + P[n].vY * P[n].vY + P[n].vZ * P[n].vZ;
				ParV = sqrt((ParE + EnDiff) / ParE);
				P[n].vX *= ParV;
				P[n].vY *= ParV;
				P[n].vZ *= ParV;
			}
		}
	}
	else
	{
		double *E = new double[N], *EA = new double[N];
		int S;
		for (n=0, ERem = Energy - nE; n<N; n++) 
			ERem -= (EA[n] = E[n] = P[n].vX * P[n].vX + P[n].vY * P[n].vY + P[n].vZ * P[n].vZ);
		if (ERem >= 0.0) for (n=0; n<N; n++) P[n].vX = P[n].vY = P[n].vZ = 0.0; 
		else
		{
			for (S=0; ERem != 0.0; ) 
			{
				for (n=0, ERem = 0.0; n<N; n++)
				{
					if (E[n] >= EnDiff) E[n] -= EnDiff;
					else 
					{
						ERem += EnDiff - E[n];
						E[n] = 0.0;
						S++;
					}
				}
				EnDiff = ERem / double(N-S);
			}
			for (n=0; n<N; n++) 
			{
				if (E[n] == 0.0) P[n].vX = P[n].vY = P[n].vZ = 0.0;
				else
				{
					ParV = sqrt(E[n] / EA[n]);
					P[n].vX *= ParV;
					P[n].vY *= ParV;
					P[n].vZ *= ParV;
				}
			}
		}
		delete[] E;
		delete[] EA;
		if (ERem > 0.0) Energy = nE + ERem / double(2);
		else Energy = nE;
	}
	E = Energy;
	return Energy;
}

void Calculation::setPotential(const PotRole role, PotStruct &PotS)
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
        if (Pot[role][n] > UaMax) Pot[role][n] = UaMax;
    }
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
