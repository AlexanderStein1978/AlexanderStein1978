#define _USE_MATH_DEFINES

#include "window.h"

#include <QGridLayout>
#include <QFile>
#include <QLabel>
#include <QTextStream>
#include <QPainter>
#include <QPaintEvent>

#include <math.h>
#include <cstdlib>

Window::Window()
{
	QGridLayout *L = new QGridLayout(this);
	L->addWidget(Start = new QPushButton("Start", this), 0, 0);
	L->addWidget(Restart = new QPushButton("Restart", this), 0, 1);
	L->addWidget(new QLabel("Step size:", this), 0, 2);
	L->addWidget(StepE = new QLineEdit(QString::number(Calc.getStepSize(), 'f', 3), this), 0, 3);
	L->addWidget(new QLabel("Energy:", this), 0, 4);
	L->addWidget(EnE = new QLineEdit(QString::number(Calc.getEnergy(), 'f', 3), this), 0, 5);
	L->addWidget(End = new QPushButton("Quit", this), 0, 6);
	L->addWidget(Pict = new Picture(this), 1, 0, 1, 7);
	L->setRowStretch(1, 1);
	int XSize, YSize;
	Calc.getSize(XSize, YSize);
	setMinimumSize(XSize, YSize);
	setMaximumSize(XSize, YSize);
	connect(Start, SIGNAL(clicked()), this, SLOT(run()));
	connect(Restart, SIGNAL(clicked()), this, SLOT(restart()));
	connect(End, SIGNAL(clicked()), this, SLOT(close()));
	connect(&Calc, SIGNAL(PictureChanged(double*, double*, double*, bool*, int)), 
			this, SLOT(draw(double*, double*, double*, bool*, int)));
}

Window::~Window()
{

}

void Window::closeEvent(QCloseEvent* event)
{
    if (Calc.isRunning())
	{
		Calc.stop();
		Calc.wait();
	}
	event->accept();
}

void Window::draw(double* XP, double* YP, double* ZP, bool *M, int N)
{
	int w = Pict->width(), he = Pict->height(), col, n, z, MaxZ;
	double ScF;
	Calc.getScales(ScF, MaxZ);
	double zSc = 255.0 / MaxZ, xc, yc;
	QPainter Paint(Pict->getPixmap());
	Paint.eraseRect(0, 0, w, he);
	Calc.mutex.lock();
	for (n=0, z=-100; n<N; n++)
	{
		if (int(ZP[n] * zSc) != z)
		{
			col = 254 - (z = int(ZP[n] * zSc));
			if (col < 0) col = 0;
			if (col > 254) col = 254;
			Paint.setPen(QColor(col, col, col));
			Paint.setBrush(QColor(col, col, col));
		}
		xc = XP[n] * ScF;
		yc = YP[n] * ScF;
		if (M[n]) Paint.drawEllipse(QRectF(xc - 2.0, yc - 2.0, 4.0, 4.0));
		else Paint.drawEllipse(QRectF(xc - 5.0, yc - 5.0, 10.0, 10.0));
	}
	Calc.mutex.unlock();
	Pict->repaint();
}

void Window::run()
{
	if (Calc.isRunning())
	{
		Calc.stop();
		Start->setText("Start");
	}
	else
	{
		EnE->setText(QString::number(Calc.setEnergy(EnE->text().toDouble()), 'f', 3));
		Calc.setStepSize(StepE->text().toDouble());
		Calc.start();
		Start->setText("Stop");
	}
}

void Window::restart()
{
	if (Calc.isRunning())
	{
		Calc.stop();
		Calc.wait();
	}
	Calc.initialize();
	run();
}

Calculation::Calculation(QObject* parent): QThread(parent)
{
	printf("Calculation::Calculation\n");
	nx = ny = nz = 10;
	double IntDist = 20.0;
	h = 0.001; //DefaultStepSize 
	Energy = 0.0; //Default energy per particle
	XS = YS = ZS = 20;
	ScF = 5.0; //Scaling factor
	QString PotFile = "Potential.dat";
	GridSizeDiv = 5;
	
	MaxX = IntDist * double(XS / GridSizeDiv);
	MaxY = IntDist * double(YS / GridSizeDiv);
	MaxZ = IntDist * double(ZS / GridSizeDiv);
	
	int n, x, y;
	double ddR;
	QFile Datei(PotFile);
	Datei.open(QIODevice::ReadOnly);
	QTextStream S(&Datei);
	QStringList SL = S.readLine().split(' ', QString::SkipEmptyParts);
	NPot = SL[0].toInt();
	Re = SL[1].toDouble();
	Pot = new double[NPot];
	dPdR = new double[NPot];
	SL = S.readLine().split(' ', QString::SkipEmptyParts);
	Rm = SL[0].toDouble();
	Pot[0] = SL[1].toDouble();
	SL = S.readLine().split(' ', QString::SkipEmptyParts);
	Pot[1] = SL[1].toDouble();
	dPdR[0] = (Pot[1] - Pot[0]) * (ddR = 1.0 / (SL[0].toDouble() - Rm));
	for (n=2; n < NPot; n++)
	{
		Pot[n] = (SL = S.readLine().split(' ', QString::SkipEmptyParts))[1].toDouble();
		dPdR[n-1] = (Pot[n] - Pot[n-1]) * ddR;
	}
	dPdR[NPot - 1] = 0.0;
	for (n=0; n < 100; n++) printf("dPdR[%d]=%f\n", n, dPdR[n]);
	RM = SL[0].toDouble();
	PS = double(NPot - 1) / (RM - Rm);
	PXS = int(round((MaxX - 0.5 * Re) / Re));
	PZS = int(round(MaxZ * 2.0 / (sqrt(3.0) * Re)));
	N = 2 * PXS * PZS;
	P = new Particle[N];
	D = new Particle*[N];
	G = new Particle***[XS];
	for (x=0; x < XS; x++)
	{
		G[x] = new Particle**[YS];
		for (y=0; y < YS; y++) G[x][y] = new Particle*[ZS];
	}
	
	MXS = int(round((MaxX - 0.125 * Re) * 4.0 / Re));
	MZS = int(round(MaxZ * 8.0 / (sqrt(3.0) * Re)));
	MD = new int[NM = MXS * MZS];
	for (n=0, MG = new int*[MXS]; n < MXS; n++) MG[n] = new int[MZS];
	MPD = 40.0;
	MPPD = 2.0;
	MXP = new double[NM];
	MYP = new double[NM];
	MZP = new double[NM];
	MX = new double[NM];
	MY = new double[NM];
	MZ = new double[NM];
	MvX = new double[NM];
	MvY = new double[NM];
	MvZ = new double[NM];
	MRSc = Re / (MaxX / double(2 * MXS + 1) + MaxZ / double(2 * MZS));
	MPRSc = 2.0;
	
	XP = new double[N + NM];
	YP = new double[N + NM];
	ZP = new double[N + NM];
	M = new bool[N + NM];
	
	initialize();
}

Calculation::~Calculation()
{
	int x, y;
	delete[] Pot;
	delete[] dPdR;
	delete[] P;
	delete[] D;
	delete[] XP;
	delete[] YP;
	delete[] ZP;
	delete[] M;
	delete[] MD;
	for (x=0; x < XS; x++)
	{
		for (y=0; y < YS; y++) delete[] G[x][y];
		delete[] G[x];
	}
	delete[] G;
	for (x=0; x < MXS; x++) delete[] MG[x];
	delete[] MG;
	delete[] MXP;
	delete[] MYP;
	delete[] MZP;
	delete[] MX;
	delete[] MY;
	delete[] MZ;
	delete[] MvX;
	delete[] MvY;
	delete[] MvZ;
}

void Calculation::geta(double* tx, double* ty, double* tz, double* ax, double* ay, double* az,
					   double *mtx, double *mty, double *mtz, 
					   double *max, double *may, double *maz)
{
	//printf("geta\n");
	int mx, my, mz, lx, ly, lz, i1, i2, p, n;
	double r, b, dx, dy, dz, a, GIR = RM / MRSc, GPIR = RM / MPRSc, Md = MPD * MRSc * MRSc;
	double ZGF = double(MZS) / MaxZ, XGF = double(MXS) / MaxX, MPd = MPPD * MPRSc * MPRSc;
	int gr = int(2.0 * GIR * MaxZ / double(MZS)), gpr = int(2.0 * GPIR * MaxZ / double(MZS)); 
	Particle *P1, *P2;
	for (mx = 0; mx < N; mx++) ax[mx] = ay[mx] = az[mx] = 0.0;
	for (mx = 0; mx < NM; mx++) max[mx] = may[mx] = maz[mx] = 0.0;
	for (mz = 0; mz < ZS; mz++) for (my = 0; my < YS; my++) for (mx = 0; mx < XS; mx++)
	{
		if (G[mx][my][mz] == 0) continue;
		//printf("l0\n");
		for (P1 = G[mx][my][mz]; P1->next != 0; P1 = P1->next)
			for (P2 = P1->next; P2 != 0; P2 = P2->next)
		{
			i1 = P1 - P;
			i2 = P2 - P;
			dx = tx[i2] - tx[i1];
			dy = ty[i2] - ty[i1];
			dz = tz[i2] - tz[i1];
			r = sqrt(dx * dx + dy * dy + dz * dz);
			p = int((r - Rm) * PS);
			if (p < 0 || p >= NPot) continue;
			//printf("p=%d\n", p);
			a = dPdR[p] / r;
			U += Pot[p];
			ax[i1] += (b = a * dx);
			ax[i2] -= b;
			ay[i1] += (b = a * dy);
			ay[i2] -= b;
			az[i1] += (b = a * dz);
			az[i2] -= b;
		}
		//printf("l1\n");
		for (lz = mz; lz < ZS && lz <= mz + GridSizeDiv; lz++)
			for (ly = (lz > mz ? ((ly = my - GridSizeDiv) >= 0 ? ly : 0) : my); 
					ly < YS && ly <= my + GridSizeDiv; ly++)
				for (lx = (ly > my || lz > mz ? ((lx = mx - GridSizeDiv) >= 0 ? lx : 0) : mx + 1);
						lx < XS && lx <= mx + GridSizeDiv; lx++)
					for (P2 = G[lx][ly][lz]; P2 != 0; P2 = P2->next)
						for (P1 = G[mx][my][mz]; P1 != 0; P1 = P1->next)
		{
			i1 = P1 - P;
			i2 = P2 - P;
			//printf("i1=%d, i2=%d\n", i1, i2);
			dx = tx[i2] - tx[i1];
			dy = ty[i2] - ty[i1];
			dz = tz[i2] - tz[i1];
			//printf("tx1=%f, tx2=%f, ty1=%f, ty2=%f, tz1=%f, tz2=%f\n", 
				//   tx[i1], tx[i2], ty[i1], ty[i2], tz[i1], tz[i2]);
			r = sqrt(dx * dx + dy * dy + dz * dz);
			p = int((r - Rm) * PS);
			if (p < 0 || p >= NPot) continue;
			//printf("p=%d, r=%f, Rm=%f, PS=%f\n", p, r, Rm, PS);
			a = dPdR[p] / r;
			U += Pot[p];
			ax[i1] += (b = a * dx);
			ax[i2] -= b;
			ay[i1] += (b = a * dy);
			ay[i2] -= b;
			az[i1] += (b = a * dz);
			az[i2] -= b;
		}
		//printf("l2\n");
	}
	printf("U=%f\n", U);
	for (n=0; n<N; n++) 
	{
		mx = ((mx = int(P[n].X * XGF)) >= 0 ? (mx < MXS ? mx : MXS - 1) : 0);
		mz = ((mz = int(P[n].Z * ZGF)) >= 0 ? (mz < MZS ? mz : MZS - 1) : 0);
		if (fabs(P[n].Y - MY[MG[mx][mz]]) > 2.0 * GPIR) continue;
		for (lx = ((lx = mx - gpr) >= 0 ? lx : 0); lx <= mx + gpr && lx < MXS; lx++)
			for (lz = ((lz = mz - gpr) >= 0 ? lz : 0); lz <= mz + gpr && lz < MZS; lz++)
		{
			dx = tx[n] - mtx[i1 = MG[lx][lz]];
			dy = ty[n] - mty[i1];
			dz = tz[n] - mtz[i1];
			r = sqrt(dx * dx + dy * dy + dz * dz) * MPRSc;
			p = int((r - Rm) * PS);
			if (p < 0 || p >= NPot) continue;
			a = dPdR[p] * MPd / r;
			U += Pot[p] * MPPD;
			ax[n] -= (dx *= a);
			ay[n] -= (dy *= a);
			az[n] -= (dz *= a);
			if (lz > 0 && lz < MZS - 1 && (lx > 0 || lx % 2 == 0) 
				&& (lx < MXS - 1 || lx % 2 != 0))
			{
				max[i1] += dx;
				may[i1] += dy;
				maz[i1] += dz;
			}
		}
	}
	printf("U=%f\n", U);
	for (mz = 0; mz < MZS; mz++) for (mx = 0; mx < MXS; mx++) 
		for (lz = mz; lz < MZS && lz <= mz + gr; lz++)
			for (lx = (lz > mz ? ((lx = mx - gr) >= 0 ? lx : 0) : mx + 1); 
					lx < MXS && lx <= mx + gr; lx++)
	{
		i1 = MG[lx][lz];
		i2 = MG[mx][mz];
		dx = mtx[i2] - mtx[i1];
		dy = mty[i2] - mty[i1];
		dz = mtz[i2] - mtz[i1];
		r = sqrt(dx * dx + dy * dy + dz * dz) * MRSc;
		p = int((r - Rm) * PS);
		if (p < 0 || p >= NPot) continue;
		a = dPdR[p] * Md / r;
		U += Pot[p] * MPD;
		dx *= a;
		dy *= a;
		dz *= a;
		if (lz > 0 && lz < MZS - 1 && (lx > 0 || lx % 2 == 0) 
				&& (lx < MXS - 1 || lx % 2 != 0))
		{
			max[i1] += dx;
			may[i1] += dy;
			maz[i1] += dz;
		}
		if (mz > 0 && mz < MZS - 1 && (mx > 0 || mx % 2 == 0) && (mx < MXS - 1 || mx % 2 != 0))
		{
			max[i2] -= dx;
			may[i2] -= dy;
			maz[i2] -= dz;
		}
	}
	T *= 0.5;
	printf("U=%f, T=%f, U+T=%f, E=%f\n", U, T, U+T, E);
	//printf("End geta\n");
}

double Calculation::getEnergy()
{
	int mx, my, mz, lx, ly, lz, i1, i2, p, n;
	double r, dx, dy, dz, T, U, GIR = RM / MRSc, GPIR = RM / MPRSc;
	double ZGF = double(MZS) / MaxZ, XGF = double(MXS) / MaxX;
	int gr = int(2.0 * GIR * MaxZ / double(MZS)), gpr = int(2.0 * GPIR * MaxZ / double(MZS)); 
	Particle *P1, *P2;
	for (mx = 0, T = U = 0.0; mx < N; mx++) 
		T += P[mx].vX * P[mx].vX + P[mx].vY * P[mx].vY + P[mx].vZ * P[mx].vZ;
	for (mx = 0; mx < NM; mx++) T += MX[mx] * MX[mx] + MY[mx] * MY[mx] + MZ[mx] * MZ[mx];
	for (mz = 0; mz < ZS; mz++) for (my = 0; my < YS; my++) for (mx = 0; mx < XS; mx++)
	{
		if (G[mx][my][mz] == 0) continue;
		//printf("l0\n");
		for (P1 = G[mx][my][mz]; P1->next != 0; P1 = P1->next)
			for (P2 = P1->next; P2 != 0; P2 = P2->next)
		{
			i1 = P1 - P;
			i2 = P2 - P;
			dx = P[i2].X - P[i1].X;
			dy = P[i2].Y - P[i1].Y;
			dz = P[i2].Z - P[i1].Z;
			r = sqrt(dx * dx + dy * dy + dz * dz);
			if (r >= RM || r <= Rm) continue;
			U += Pot[p = int((r - Rm) * PS)];
		}
		//printf("l1\n");
		for (lz = mz; lz < ZS && lz < mz + GridSizeDiv; lz++)
			for (ly = (lz > mz ? ((ly = my - GridSizeDiv) >= 0 ? ly : 0) : my); 
					ly < YS && ly < my + GridSizeDiv; ly++)
				for (lx = (ly > my || lz > mz ? ((lx = mx - GridSizeDiv) >= 0 ? lx : 0) : mx);
						lx < XS && lx < mx + GridSizeDiv; lx++)
					for (P2 = G[lx][ly][lz]; P2 != 0; P2 = P2->next)
						for (P1 = G[mx][my][mz]; P1 != 0; P1 = P1->next)
		{
			i1 = P1 - P;
			i2 = P2 - P;
			dx = P[i2].X - P[i1].X;
			dy = P[i2].Y - P[i1].Y;
			dz = P[i2].Z - P[i1].Z;
			r = sqrt(dx * dx + dy * dy + dz * dz);
			if (r >= RM || r <= Rm) continue;
			U += Pot[p = int((r - Rm) * PS)];
		}
		//printf("l2\n");
	}
	for (n=0; n<N; n++) 
	{
		mx = ((mx = int(P[n].X * XGF)) >= 0 ? (mx < MXS ? mx : MXS - 1) : 0);
		mz = ((mz = int(P[n].Z * ZGF)) >= 0 ? (mz < MZS ? mz : MZS - 1) : 0);
		if (fabs(P[n].Y - MY[MG[mx][mz]]) > 2.0 * GPIR) continue;
		for (lx = ((lx = mx - gpr) >= 0 ? lx : 0); lx <= mx + gpr && lx < MXS; lx++)
			for (lz = ((lz = mz - gpr) >= 0 ? lz : 0); lz <= mz + gpr && lz < MZS; lz++)
		{
			dx = P[n].X - MX[i1 = MG[lx][lz]];
			dy = P[n].Y - MY[i1];
			dz = P[n].Z - MZ[i1];
			r = sqrt(dx * dx + dy * dy + dz * dz) * MPRSc;
			if (r >= RM || r <= Rm) continue;
			p = int((r - Rm) * PS);
			U += Pot[p] * MPPD;
		}
	}
	for (mz = 0; mz < MZS; mz++) for (mx = 0; mx < MXS; mx++) 
		for (lz = mz; lz < MZS && lz <= mz + gr; lz++)
			for (lx = (lz > mz ? ((lx = mx - gr) >= 0 ? lx : 0) : mx); 
					lx < MXS && lx <= mx + gr; lx++)
	{
		i1 = MG[lx][lz];
		i2 = MG[mx][mz];
		dx = MX[i2] - MX[i1];
		dy = MY[i2] - MY[i1];
		dz = MZ[i2] - MZ[i1];
		r = sqrt(dx * dx + dy * dy + dz * dz) * MRSc;
		if (r >= RM || r <= Rm) continue;
		p = int((r - Rm) * PS);
		U += Pot[p] * MPD;
	}
	return Energy = T+U;
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

double Calculation::getStepSize()
{
	return h;
}

void Calculation::initialize()
{
	printf("Calculation::initialize\n");
	int n, x, y, z, b;
	double XF = double(XS) / MaxX, YF = double(YS) / MaxY, ZF = double(ZS) / MaxZ;
	double rx, ry, rz, rys, rxs; 
	rxs = MaxX / (double(PXS) + 0.5);
	rys = MaxZ / double(PZS);
	for (x=0; x < XS; x++) for (y=0; y < YS; y++) for (z=0; z < ZS; z++) G[x][y][z] = 0;
	for (n=z=0, rz = 0.0, ry = 0.5 * (MaxY - Re); z < PZS; z++, rz += rys) 
		for (x=0, rx = (z % 2 == 0 ? 0.0 : 0.5 * rxs); x < PXS; x++, n+=2, rx += rxs)
	{
		P[n].vX = P[n].vY = P[n].vZ = P[n+1].vX = P[n+1].vY = P[n+1].vZ = 0.0;
		P[n].X = P[n+1].X = rx;
		P[n].Y = ry;
		P[n+1].Y = ry + Re;
		P[n].Z = P[n+1].Z = rz;
		P[n].xp = P[n+1].xp = ((b = int(XF * rx)) >= 0 ? (b < XS ? b : XS - 1) : 0);
		P[n].yp = ((b = int(YF * ry)) >= 0 ? (b < YS ? b : YS - 1) : 0);
		P[n+1].yp = ((b = int(YF * (ry + Re))) >= 0 ? (b < YS ? b : YS - 1) : 0);
		P[n].zp = P[n+1].zp = ((b = int(ZF * rz)) >= 0 ? (b < ZS ? b : ZS - 1) : 0);
		P[n].prev = P[n+1].prev = 0;
		P[n].next = G[P[n].xp][P[n].yp][P[n].zp];
		G[P[n].xp][P[n].yp][P[n].zp] = D[n] = P+n;
		if (P[n].next != 0) P[n].next->prev = D[n];
		P[n+1].next = G[P[n+1].xp][P[n+1].yp][P[n+1].zp];
		G[P[n+1].xp][P[n+1].yp][P[n+1].zp] = D[n+1] = P+n+1;
		if (P[n+1].next != 0) P[n+1].next->prev = D[n+1];
	}
	rxs = MaxX / (double(MXS) + 0.5);
	rys = MaxZ / double(MZS);
	for (n=z=0, rz = 0.0, ry = 0.5 * MaxY; z < MZS; z++, rz += rys) 
		for (x=0, rx = (z % 2 == 0 ? 0.0 : 0.5 * rxs); x < MXS; x++, n++, rx += rxs)
	{
		MvX[n] = MvY[n] = MvZ[n] = 0.0;
		MX[n] = rx;
		MY[n] = ry;
		MZ[n] = rz;
		MG[x][z] = MD[n] = n;
	}
}

void Calculation::run()
{
    // Contains modified version of the rk4 algorithm from Numerical Recipes, Third Edition
	int n, i, x, y, z;
	double hh = 0.5 * h, h6 = h / 6.0, *ax = new double[N], *ay = new double[N];
	double *az = new double[N], *dxm = new double[N], *dym = new double[N], *dzm = new double[N];
	double *dvxm = new double[N], *dvym = new double[N], *dvzm = new double[N];
	double *dxt = new double[N], *dyt = new double[N], *dzt = new double[N];
	double *dvxt = new double[N], *dvyt = new double[N], *dvzt = new double[N];
	double *xt = new double[N], *yt = new double[N], *zt = new double[N], vF;
	double XF = double(XS) / MaxX, YF = double(YS) / MaxY, ZF = double(ZS) / MaxZ;
	double *max = new double[NM], *may = new double[NM], *maz = new double[NM];
	double *mdxm = new double[NM], *mdym = new double[NM], *mdzm = new double[NM];
	double *mdvxm = new double[NM], *mdvym = new double[NM], *mdvzm = new double[NM];
	double *mdxt = new double[NM], *mdyt = new double[NM], *mdzt = new double[NM];
	double *mdvxt = new double[NM], *mdvyt = new double[NM], *mdvzt = new double[NM];
	double *mxt = new double[NM], *myt = new double[NM], *mzt = new double[NM];
	Particle *PB;
	for (i=0, Run = true; Run; i++)
	{
		printf("i=%d, ", i);
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
				if (G[x][y][z] != 0) G[x][y][z]->prev = P+n;
				G[x][y][z] = P+n;
				P[n].xp = x;
				P[n].yp = y;
				P[n].zp = z;
			}
			xt[n] = P[n].X;
			yt[n] = P[n].Y;
			zt[n] = P[n].Z;
		}
		for (n=0, U = T = 0.0; n<N; n++) 
			T += P[n].vX * P[n].vX + P[n].vY * P[n].vY + P[n].vZ * P[n].vZ;
		for (n=0; n < NM; n++) T += MvX[n] * MvX[n] + MvY[n] * MvY[n] + MvZ[n] * MvZ[n];
		geta(xt, yt, zt, ax, ay, az, MX, MY, MZ, max, may, maz);
		if (E == 0.0 || T == 0.0) E = T + U;
		else if (T + U != E && E > U)
		{
			vF = sqrt((E - U) / T);
			for (n=0; n<N; n++)
			{
				P[n].vX *= vF;
				P[n].vY *= vF;
				P[n].vZ *= vF;
			}
			for (n=0; n < NM; n++)
			{
				MvX[n] *= vF;
				MvY[n] *= vF;
				MvZ[n] *= vF;
			}
		}
		for (n=0; n<N; n++)
		{
			xt[n] = P[n].X + hh * P[n].vX;
			yt[n] = P[n].Y + hh * P[n].vY;
			zt[n] = P[n].Z + hh * P[n].vZ;
			dxt[n] = P[n].vX + hh * ax[n];
			dyt[n] = P[n].vY + hh * ay[n];
			dzt[n] = P[n].vZ + hh * az[n];
		}
		for (n=0; n < NM; n++)
		{
			mxt[n] = MX[n] + hh * MvX[n];
			myt[n] = MY[n] + hh * MvY[n];
			mzt[n] = MZ[n] + hh * MvZ[n];
			mdxt[n] = MvX[n] + hh * max[n];
			mdyt[n] = MvY[n] + hh * may[n];
			mdzt[n] = MvZ[n] + hh * maz[n];
		}
		for (n=0, U = T = 0.0; n<N; n++) 
			T += dxt[n] * dxt[n] + dyt[n] * dyt[n] + dzt[n] * dzt[n];
		for (n=0; n < NM; n++) T += mdxt[n] * mdxt[n] + mdyt[n] * mdyt[n] + mdzt[n] * mdzt[n];
		geta(xt, yt, zt, dvxt, dvyt, dvzt, mxt, myt, mzt, mdvxt, mdvyt, mdvzt);
		for (n=0; n<N; n++)
		{
			xt[n] = P[n].X + hh * dxt[n];
			yt[n] = P[n].Y + hh * dyt[n];
			zt[n] = P[n].Z + hh * dzt[n];
			dxm[n] = P[n].vX + hh * dvxt[n];
			dym[n] = P[n].vY + hh * dvyt[n];
			dzm[n] = P[n].vZ + hh * dvzt[n];
		}
		for (n=0; n < NM; n++)
		{
			mxt[n] = MX[n] + hh * mdxt[n];
			myt[n] = MY[n] + hh * mdyt[n];
			mzt[n] = MZ[n] + hh * mdzt[n];
			mdxm[n] = MvX[n] + hh * mdvxt[n];
			mdym[n] = MvY[n] + hh * mdvyt[n];
			mdzm[n] = MvZ[n] + hh * mdvzt[n];
		}
		for (n=0, U = T = 0.0; n<N; n++) 
			T += dxm[n] * dxm[n] + dym[n] * dym[n] + dzm[n] * dzm[n];
		for (n=0; n < NM; n++) T += mdxm[n] * mdxm[n] + mdym[n] * mdym[n] + mdzm[n] * mdzm[n];
		geta(xt, yt, zt, dvxm, dvym, dvzm, mxt, myt, mzt, mdvxm, mdvym, mdvzm);
		for (n=0; n<N; n++)
		{
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
		for (n=0; n < NM; n++)
		{
			mxt[n] = MX[n] + h * mdxm[n];
			myt[n] = MY[n] + h * mdym[n];
			mzt[n] = MZ[n] + h * mdzm[n];
			mdxm[n] += mdxt[n];
			mdym[n] += mdyt[n];
			mdzm[n] += mdzt[n];
			mdxt[n] = MvX[n] + h * mdvxm[n];
			mdyt[n] = MvY[n] + h * mdvym[n];
			mdzt[n] = MvZ[n] + h * mdvzm[n];
			mdvxm[n] += mdvxt[n];
			mdvym[n] += mdvyt[n];
			mdvzm[n] += mdvzt[n];
		}
		for (n=0, U = T = 0.0; n<N; n++) 
			T += dxt[n] * dxt[n] + dyt[n] * dyt[n] + dzt[n] * dzt[n];
		for (n=0; n < NM; n++) T += mdxt[n] * mdxt[n] + mdyt[n] * mdyt[n] + mdzt[n] * mdzt[n];
		geta(xt, yt, zt, dvxt, dvyt, dvzt, mxt, myt, mzt, mdvxt, mdvyt, mdvzt);
		for (n=0; n<N; n++)
		{
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
				printf("Particel %d is nan!\n", n);
				Run = false;
			}*/
		}
		for (n=0; n < NM; n++)
		{
			MX[n] += h6 * (MvX[n] + mdxt[n] + 2.0 * mdxm[n]);
			MY[n] += h6 * (MvY[n] + mdyt[n] + 2.0 * mdym[n]);
			MZ[n] += h6 * (MvZ[n] + mdzt[n] + 2.0 * mdzm[n]);
			MvX[n] += h6 * (max[n] + mdvxt[n] + 2.0 * mdvxm[n]);
			MvY[n] += h6 * (may[n] + mdvyt[n] + 2.0 * mdvym[n]);
			MvZ[n] += h6 * (maz[n] + mdvzt[n] + 2.0 * mdvzm[n]);
			/*if (isnan(MX[n]) || isnan(MY[n]) || isnan(MZ[n]) || isnan(MvX[n]) || isnan(MvY[n]) || isnan(MvZ[n]))
			{
				printf("Mesh particel %d is nan!\n", n);
				Run = false;
			}*/
		}
		for (PB = P; PB != 0; ) for (n=1, PB = 0; n<N; n++) if (D[n]->Z < D[n-1]->Z)
		{
			PB = D[n];
			D[n] = D[n-1];
			D[n-1] = PB;
		}
		for (x=0; x!=-1; ) for (n=1, x=-1; n<N; n++) if (MZ[MD[n]] < MZ[MD[n-1]])
		{
			x = MD[n];
			MD[n] = MD[n-1];
			MD[n-1] = x;
		}
		mutex.lock();
		for (n=x=y=0; n < N + NM; n++)
		{
			if (x < N ? (y < NM ? D[x]->Z < MZ[MD[y]] : true) : false)
			{
				XP[n] = D[x]->X;
				YP[n] = D[x]->Y;
				ZP[n] = D[x++]->Z;
				M[n] = false;
			}
			else
			{
				XP[n] = MX[MD[y]];
				YP[n] = MY[MD[y]];
				ZP[n] = MZ[MD[y++]];
				M[n] = true;
			}
		}
		mutex.unlock();
		emit PictureChanged(XP, YP, ZP, M, N + NM);
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
	delete[] max;
	delete[] may;
	delete[] maz;
	delete[] mdxm;
	delete[] mdym;
	delete[] mdzm;
	delete[] mdvxm;
	delete[] mdvym;
	delete[] mdvzm;
	delete[] mdxt;
	delete[] mdyt;
	delete[] mdzt;
	delete[] mdvxt;
	delete[] mdvyt;
	delete[] mdvzt;
	delete[] mxt;
	delete[] myt;
	delete[] mzt;
}

double Calculation::setEnergy(double newE)
{
	int n;
	double nE = newE, ParE, ParV, ERem = 0.0;
	double VC, A1, A2, RD = M_PI / RAND_MAX, EnDiff = 2.0 * (nE - Energy) / N;
	if (nE > Energy)
	{
		VC = sqrt(EnDiff);
		Energy = nE;
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
		for (n=0; n<N; n++)
		{
			if ((ParE = P[n].vX * P[n].vX + P[n].vY * P[n].vY + P[n].vZ * P[n].vZ) + EnDiff > 0.0)
			{
				ParV = sqrt((ParE + EnDiff) / ParE);
				P[n].vX *= ParV;
				P[n].vY *= ParV;
				P[n].vZ *= ParV;
			}
			else 
			{
				P[n].vX = P[n].vY = P[n].vZ = 0.0;
				ERem -= ParE + EnDiff;
			}
		}
		if (ERem > 0.0) Energy = nE + ERem / double(2);
		else Energy = nE;
	}
	E = Energy;
	return Energy;
}

void Calculation::setStepSize(double nh)
{
	h = nh;
}

void Calculation::stop()
{
	Run = false;
}

Picture::Picture(QWidget* parent): QWidget(parent)
{
	Map = 0;
}

Picture::~Picture()
{
	if (Map != 0) delete Map;
}

QPixmap *Picture::getPixmap()
{
	if (Map == 0) Map = new QPixmap(width(), height());
	return Map;
}

void Picture::paintEvent(QPaintEvent *event)
{
	if (Map == 0) return;
	QPainter P(this);
	QRect R = event->rect();
	P.drawPixmap(R, *Map, R);
}
