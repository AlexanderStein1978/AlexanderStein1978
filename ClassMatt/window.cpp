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
	connect(&Calc, SIGNAL(PictureChanged(double *, double *, double *, int)), 
			this, SLOT(draw(double*,double*,double*,int)));
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

void Window::draw(double* XP, double* YP, double* ZP, int N)
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
		Paint.drawEllipse(QRectF(xc - 5.0, yc - 5.0, 10.0, 10.0));
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
	
	int n, x, y, l;
	double ddR;
	QFile Datei(PotFile);
	Datei.open(QIODevice::ReadOnly);
	QTextStream S(&Datei);
	QStringList SL = S.readLine().split(' ', QString::SkipEmptyParts);
	l = SL[0].toInt();
	Re = SL[1].toDouble();
	Pot = new double[l];
	dPdR = new double[l];
	SL = S.readLine().split(' ', QString::SkipEmptyParts);
	Rm = SL[0].toDouble();
	Pot[0] = SL[1].toDouble();
	SL = S.readLine().split(' ', QString::SkipEmptyParts);
	Pot[1] = SL[1].toDouble();
	dPdR[0] = (Pot[1] - Pot[0]) * (ddR = 1.0 / (SL[0].toDouble() - Rm));
	for (n=2; n<l; n++)
	{
		Pot[n] = (SL = S.readLine().split(' ', QString::SkipEmptyParts))[1].toDouble();
		dPdR[n-1] = (Pot[n] - Pot[n-1]) * ddR;
	}
	dPdR[l-1] = 0.0;
	RM = SL[0].toDouble();
	PS = double(l-1) / (RM - Rm);
	N = nx * ny * nz;
	E = 0.0;
	P = new Particle[N];
	D = new Particle*[N];
	G = new Particle***[XS];
	XP = new double[N];
	YP = new double[N];
	ZP = new double[N];
	for (x=0; x < XS; x++)
	{
		G[x] = new Particle**[YS];
		for (y=0; y < YS; y++) G[x][y] = new Particle*[ZS];
	}
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
	for (x=0; x < XS; x++)
	{
		for (y=0; y < YS; y++) delete[] G[x][y];
		delete[] G[x];
	}
	delete[] G;
}

void Calculation::geta(double* tx, double* ty, double* tz, double* ax, double* ay, double* az)
{
	//printf("geta\n");
	int mx, my, mz, lx, ly, lz, i1, i2, p;
	double r, b, dx, dy, dz, a;
	Particle *P1, *P2;
	for (mx = 0; mx < N; mx++) ax[mx] = ay[mx] = az[mx] = 0.0;
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
			if (r >= RM || r <= Rm) continue;
			p = int((r - Rm) * PS);
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
		for (lz = mz; lz < ZS && lz < mz + GridSizeDiv; lz++)
			for (ly = (lz > mz ? ((ly = my - GridSizeDiv) >= 0 ? ly : 0) : my); 
					ly < YS && ly < my + GridSizeDiv; ly++)
				for (lx = (ly > my || lz > mz ? ((lx = mx - GridSizeDiv) >= 0 ? lx : 0) : mx + 1);
						lx < XS && lx < mx + GridSizeDiv; lx++)
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
			if (r >= RM || r <= Rm) continue;
			p = int((r - Rm) * PS);
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
	T *= 0.5;
	printf("U=%f, T=%f, E=%f\n", U, T, U+T);
	//printf("End geta\n");
}

double Calculation::getEnergy()
{
	int mx, my, mz, lx, ly, lz, i1, i2, p;
	double r, dx, dy, dz, T, U;
	Particle *P1, *P2;
	for (mx = 0, T = U = 0.0; mx < N; mx++) 
		T += P[mx].vX * P[mx].vX + P[mx].vY * P[mx].vY + P[mx].vZ * P[mx].vZ;
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
	return Energy = (T+U)/N;
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
	rxs = 0.5 * (MaxX - Re * double(nx - 1));
	rys = 0.5 * (MaxY - Re * double(ny - 1));
	rz = 0.5 * (MaxZ - Re * double(nz - 1));
	for (x=0; x < XS; x++) for (y=0; y < YS; y++) for (z=0; z < ZS; z++) G[x][y][z] = 0;
	for (n=z=0; z < nz; z++, rz += Re) for (y=0, ry = rys; y < ny; y++, ry += Re)
		for (x=0, rx = rxs; x < nx; x++, n++, rx += Re)
	{
		P[n].vX = P[n].vY = P[n].vZ = 0.0;
		P[n].X = rx;
		P[n].Y = ry;
		P[n].Z = rz;
		P[n].xp = ((b = int(XF * rx)) >= 0 ? (b < XS ? b : XS - 1) : 0);
		P[n].yp = ((b = int(YF * ry)) >= 0 ? (b < YS ? b : YS - 1) : 0);
		P[n].zp = ((b = int(ZF * rz)) >= 0 ? (b < ZS ? b : ZS - 1) : 0);
		P[n].prev = 0;
		P[n].next = G[P[n].xp][P[n].yp][P[n].zp];
		G[P[n].xp][P[n].yp][P[n].zp] = D[n] = P+n;
		if (P[n].next != 0) P[n].next->prev = D[n];
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
		geta(xt, yt, zt, ax, ay, az);
		if (E == 0.0 || T == 0.0) E = T + U;
		else if (T + U != E)
		{
			vF = sqrt((E - U) / T);
			for (n=0; n<N; n++)
			{
				P[n].vX *= vF;
				P[n].vY *= vF;
				P[n].vZ *= vF;
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
		for (n=0, U = T = 0.0; n<N; n++) 
			T += dxt[n] * dxt[n] + dyt[n] * dyt[n] + dzt[n] * dzt[n];
		geta(xt, yt, zt, dvxt, dvyt, dvzt);
		for (n=0; n<N; n++)
		{
			xt[n] = P[n].X + hh * dxt[n];
			yt[n] = P[n].Y + hh * dyt[n];
			zt[n] = P[n].Z + hh * dzt[n];
			dxm[n] = P[n].vX + hh * dvxt[n];
			dym[n] = P[n].vY + hh * dvyt[n];
			dzm[n] = P[n].vZ + hh * dvzt[n];
		}
		for (n=0, U = T = 0.0; n<N; n++) 
			T += dxm[n] * dxm[n] + dym[n] * dym[n] + dzm[n] * dzm[n];
		geta(xt, yt, zt, dvxm, dvym, dvzm);
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
		for (n=0, U = T = 0.0; n<N; n++) 
			T += dxt[n] * dxt[n] + dyt[n] * dyt[n] + dzt[n] * dzt[n];
		geta(xt, yt, zt, dvxt, dvyt, dvzt);
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
		}
		for (PB = P; PB != 0; ) for (n=1, PB = 0; n<N; n++) if (D[n]->Z < D[n-1]->Z)
		{
			PB = D[n];
			D[n] = D[n-1];
			D[n-1] = PB;
		}
		mutex.lock();
		for (n=0; n<N; n++)
		{
			XP[n] = D[n]->X;
			YP[n] = D[n]->Y;
			ZP[n] = D[n]->Z;
		}
		mutex.unlock();
		emit PictureChanged(XP, YP, ZP, N);
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
}

double Calculation::setEnergy(double newE)
{
	int n;
	double nE = newE, ParE, ParV, ERem = 0.0;
	double VC, A1, A2, RD = M_PI / RAND_MAX, EnDiff = 2.0 * (nE - Energy);
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
		if (ERem > 0.0) Energy = nE + ERem / double(2 * N);
		else Energy = nE;
	}
	E = N * Energy;
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
