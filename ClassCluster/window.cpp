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
	L->addWidget(Rotate = new QPushButton("Rotate", this), 0, 2);
	L->addWidget(Move = new QPushButton("Move", this), 0, 3);
	L->addWidget(new QLabel("Speed:", this), 0, 4);
	L->addWidget(Speed = new QLineEdit(QString::number(Calc.getSpeed(), 'f', 3), this), 0, 5);
	L->addWidget(new QLabel("Step size:", this), 0, 6);
	L->addWidget(StepE = new QLineEdit(QString::number(Calc.getStepSize(), 'f', 3), this), 0, 7);
	L->addWidget(new QLabel("Energy:", this), 0, 8);
	L->addWidget(EnE = new QLineEdit(QString::number(Calc.getEnergy(), 'f', 3), this), 0, 9);
	L->addWidget(End = new QPushButton("Quit", this), 0, 10);
	L->addWidget(Pict = new Picture(this), 1, 0, 1, 11);
	L->setRowStretch(1, 1);
	int XSize, YSize;
	Calc.getSize(XSize, YSize);
	setMinimumSize(XSize, YSize);
	setMaximumSize(XSize, YSize);
	connect(Start, SIGNAL(clicked()), this, SLOT(run()));
	connect(Restart, SIGNAL(clicked()), this, SLOT(restart()));
	connect(Rotate, SIGNAL(clicked()), &Calc, SLOT(rotate()));
	connect(Move, SIGNAL(clicked()), this, SLOT(move()));
	connect(Speed, SIGNAL(editingFinished()), this, SLOT(speedChanged()));
	connect(End, SIGNAL(clicked()), this, SLOT(close()));
	connect(&Calc, SIGNAL(PictureChanged(double*, double*, double*, int)), 
			this, SLOT(draw(double*, double*, double*, int)));
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

void Window::move()
{
	if (Calc.getMove()) Move->setText("Move");
	else Move->setText("Hold");
	Calc.move();
}

void Window::speedChanged()
{
	Calc.setSpeed(Speed->text().toDouble());
}


Calculation::Calculation(QObject* parent): QThread(parent)
{
	printf("Calculation::Calculation\n");
	nx = ny = nz = 10;
	double IntDist = 20.0, R, F, P0, st, F0, F1;
	st = 1e-3;
	PS = 1e3; // 1/st
	NPot = 30000;
	Re = 4.0;
	R = pow(Re, 6);
	double De = 1000.0;
	double a = 2.0 * De * R;
	double b = De * R * R;
	double bRep = 8.0 * b;
	Rm = st;
	RM = st * double(NPot);
	h = 0.001; //DefaultStepSize 
	Energy = 0.0; //Default energy per particle
	XS = YS = ZS = 20;
	ScF = 10.0; //Scaling factor
	QString PotFile = "Potential.dat";
	GridSizeDiv = 5;
	rotated = false;
	Move = false;
	Speed = 1e3;
	
	MaxX = IntDist * double(XS / GridSizeDiv);
	MaxY = IntDist * double(YS / GridSizeDiv);
	MaxZ = IntDist * double(ZS / GridSizeDiv);
	
	YMid = 0.5 * MaxY;
	
	int n, x, y;
	/*QFile Datei(PotFile);
	Datei.open(QIODevice::ReadOnly);
	QTextStream S(&Datei);
	QStringList SL = S.readLine().split(' ', QString::SkipEmptyParts);
	NPot = SL[0].toInt();
	Re = SL[1].toDouble();
	SL = S.readLine().split(' ', QString::SkipEmptyParts);
	Rm = SL[0].toDouble();
	P0 = SL[1].toDouble();
	F0 = SL[2].toDouble();
	SL = S.readLine().split(' ', QString::SkipEmptyParts);
	RD = SL[0].toDouble() - Rm;
	P1 = SL[1].toDouble();
	F1 = SL[2].toDouble();
	Pot = new double[NPot += (x = int(Rm / RD)) + 1];
	dPdR = new double[NPot];
	RepF = new double[NPot];
	RepP = new double[NPot];
	for (n=0; n<=x; n++) dPdR[n] = F0;
	for (n=x-1, Pot[x] = P0 + (F = F0 * RD); n>=0; n--) Pot[n] = Pot[n+1] + F;
	for (n=x+3, Pot[x+1] = P0, Pot[x+2] = P1, dPdR[x+1] = F0, dPdR[x+2] = F1; n < NPot; n++)
	{
		Pot[n] = (SL = S.readLine().split(' ', QString::SkipEmptyParts))[1].toDouble();
		dPdR[n] = SL[2].toDouble();
	}
	RM = SL[0].toDouble();
	PS = double(NPot - x - 2) / (RM - Rm);*/
	
	Pot = new double[NPot];
	dPdR = new double[NPot];
	RepF = new double[NPot];
	RepP = new double[NPot];
	for (n=0, R = Rm; n < NPot; n++, R += st)
	{
		F = pow(R, -6.0);
		F0 = 1.0 / R;
		F1 = -a * F;
		P0 = F * F;
		Pot[n] = F1 + b * P0;
		RepP[n] = F1 + bRep * P0;
		F1 *= -6 * F0;
		P0 *= -12 * F0;
		dPdR[n] = F1 + b * P0;
		RepF[n] = F1 + bRep * P0;
	}
	
	PZS = int(round(MaxZ / Re));
	PYS = int(round(MaxY / Re)); 
	PXS = int(round(MaxX / Re));
	N = 2 * PXS * PYS;
	P = new Particle[N];
	D = new Particle*[N];
	G = new Particle***[XS];
	MAR = new double*[N];
	for (x=0; x < XS; x++)
	{
		G[x] = new Particle**[YS];
		for (y=0; y < YS; y++) G[x][y] = new Particle*[ZS];
	}
	for (n=0; n < N; n++) MAR[n] = new double[4];
	XP = new double[N];
	YP = new double[N];
	ZP = new double[N];
	Fixed = new bool[N];
	
	initialize();
}

Calculation::~Calculation()
{
	int x, y, n;
	delete[] Pot;
	delete[] dPdR;
	delete[] RepP;
	delete[] RepF;
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

void Calculation::geta(double *tx, double *ty, double *tz, double *ax, double *ay, 
					   double *az)
{
	//printf("geta\n");
	int mx, my, mz, lx, ly, lz, i1, i2, p, n;
	double r, b, dx, dy, dz, a;
	Particle *PP1, *PP2;
	for (mx = 0; mx < N; mx++) ax[mx] = ay[mx] = az[mx] = 0.0;
	for (n=0; n<N; n++) for (p=0; p<4; p++) MAR[n][p] = RM;
	for (mz = 0; mz < ZS; mz++) for (my = 0; my < YS; my++) for (mx = 0; mx < XS; mx++)
	{
		if (G[mx][my][mz] != 0)
		{
			//printf("l0\n");
			for (PP1 = G[mx][my][mz]; PP1->next != 0; PP1 = PP1->next)
				for (PP2 = PP1->next; PP2 != 0; PP2 = PP2->next)
			{
				i1 = PP1 - P;
				i2 = PP2 - P;
				dx = tx[i2] - tx[i1];
				dy = ty[i2] - ty[i1];
				dz = tz[i2] - tz[i1];
				r = sqrt(dx * dx + dy * dy + dz * dz);
				for (n=0; (n<4 ? r > MAR[i1][n] : false); n++) ;
				for (p=3; p>n; p--) MAR[i1][p] = MAR[i1][p-1];
				if (n<4) MAR[i1][n] = r;
				for (n=0; (n<4 ? r > MAR[i2][n] : false); n++) ;
				for (p=3; p>n; p--) MAR[i2][p] = MAR[i2][p-1];
				if (n<4) MAR[i2][n] = r;
			}
			//printf("l1\n");
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
				//printf("i1=%d, i2=%d\n", i1, i2);
				dx = tx[i2] - tx[i1];
				dy = ty[i2] - ty[i1];
				dz = tz[i2] - tz[i1];
				//printf("tx1=%f, tx2=%f, ty1=%f, ty2=%f, tz1=%f, tz2=%f\n", 
					//   tx[i1], tx[i2], ty[i1], ty[i2], tz[i1], tz[i2]);
				r = sqrt(dx * dx + dy * dy + dz * dz);
				for (n=0; (n<4 ? r > MAR[i1][n] : false); n++) ;
				for (p=3; p>n; p--) MAR[i1][p] = MAR[i1][p-1];
				if (n<4) MAR[i1][n] = r;
				for (n=0; (n<4 ? r > MAR[i2][n] : false); n++) ;
				for (p=3; p>n; p--) MAR[i2][p] = MAR[i2][p-1];
				if (n<4) MAR[i2][n] = r;
			}
		}
	}
	for (mz = 0; mz < ZS; mz++) for (my = 0; my < YS; my++) for (mx = 0; mx < XS; mx++)
	{
		if (G[mx][my][mz] != 0)
		{
			//printf("l0\n");
			for (PP1 = G[mx][my][mz]; PP1->next != 0; PP1 = PP1->next)
				for (PP2 = PP1->next; PP2 != 0; PP2 = PP2->next)
			{
				i1 = PP1 - P;
				i2 = PP2 - P;
				dx = tx[i2] - tx[i1];
				dy = ty[i2] - ty[i1];
				dz = tz[i2] - tz[i1];
				r = sqrt(dx * dx + dy * dy + dz * dz);
				p = int((r - Rm) * PS);
				if (p < 0 || p >= NPot) continue;
				//printf("p=%d\n", p);
				if (r <= MAR[i1][3] && r <= MAR[i2][3])
				{
					a = dPdR[p] / r;
					U += Pot[p];
				}
				else
				{
					a = RepF[p] / r;
					U += RepP[p];
				}
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
						for (PP2 = G[lx][ly][lz]; PP2 != 0; PP2 = PP2->next)
							for (PP1 = G[mx][my][mz]; PP1 != 0; PP1 = PP1->next)
			{
				i1 = PP1 - P;
				i2 = PP2 - P;
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
				if (r <= MAR[i1][3] && r <= MAR[i2][3])
				{
					a = dPdR[p] / r;
					U += Pot[p];
				}
				else
				{
					a = RepF[p] / r;
					U += RepP[p];
				}
				ax[i1] += (b = a * dx);
				ax[i2] -= b;
				ay[i1] += (b = a * dy);
				ay[i2] -= b;
				az[i1] += (b = a * dz);
				az[i2] -= b;
			}
		}
	}
	printf("U=%f\n", U);
	T *= 0.5;
	printf("U=%f, T=%f, U+T=%f, E=%f\n", U, T, U+T, E);
	//printf("End geta\n");
}

double Calculation::getEnergy()
{
	int mx, my, mz, lx, ly, lz, i1, i2, p, n;
	double r, dx, dy, dz, T, U;
	Particle *PP1, *PP2;
	for (mx = 0, T = U = 0.0; mx < N; mx++) 
		T += P[mx].vX * P[mx].vX + P[mx].vY * P[mx].vY + P[mx].vZ * P[mx].vZ;
	for (n=0; n<N; n++) for (p=0; p<4; p++) MAR[n][p] = RM;
	for (mz = 0; mz < ZS; mz++) for (my = 0; my < YS; my++) for (mx = 0; mx < XS; mx++)
	{
		if (G[mx][my][mz] != 0)
		{
			//printf("l0\n");
			for (PP1 = G[mx][my][mz]; PP1->next != 0; PP1 = PP1->next)
				for (PP2 = PP1->next; PP2 != 0; PP2 = PP2->next)
			{
				i1 = PP1 - P;
				i2 = PP2 - P;
				dx = PP2->X - PP1->X;
				dy = PP2->Y - PP1->Y;
				dz = PP2->Z - PP1->Z;
				r = sqrt(dx * dx + dy * dy + dz * dz);
				for (n=0; (n<4 ? r > MAR[i1][n] : false); n++) ;
				for (p=3; p>n; p--) MAR[i1][p] = MAR[i1][p-1];
				if (n<4) MAR[i1][n] = r;
				for (n=0; (n<4 ? r > MAR[i2][n] : false); n++) ;
				for (p=3; p>n; p--) MAR[i2][p] = MAR[i2][p-1];
				if (n<4) MAR[i2][n] = r;
			}
			//printf("l1\n");
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
				//printf("i1=%d, i2=%d\n", i1, i2);
				dx = PP2->X - PP1->X;
				dy = PP2->Y - PP1->Y;
				dz = PP2->Z - PP1->Z;
				//printf("tx1=%f, tx2=%f, ty1=%f, ty2=%f, tz1=%f, tz2=%f\n", 
					//   tx[i1], tx[i2], ty[i1], ty[i2], tz[i1], tz[i2]);
				r = sqrt(dx * dx + dy * dy + dz * dz);
				for (n=0; (n<4 ? r > MAR[i1][n] : false); n++) ;
				for (p=3; p>n; p--) MAR[i1][p] = MAR[i1][p-1];
				if (n<4) MAR[i1][n] = r;
				for (n=0; (n<4 ? r > MAR[i2][n] : false); n++) ;
				for (p=3; p>n; p--) MAR[i2][p] = MAR[i2][p-1];
				if (n<4) MAR[i2][n] = r;
			}
		}
	}
	for (mz = 0; mz < ZS; mz++) for (my = 0; my < YS; my++) for (mx = 0; mx < XS; mx++)
	{
		if (G[mx][my][mz] != 0)
		{
			//printf("l0\n");
			for (PP1 = G[mx][my][mz]; PP1->next != 0; PP1 = PP1->next)
				for (PP2 = PP1->next; PP2 != 0; PP2 = PP2->next)
			{
				i1 = PP1 - P;
				i2 = PP2 - P;
				dx = PP2->X - PP1->X;
				dy = PP2->Y - PP1->Y;
				dz = PP2->Z - PP1->Z;
				r = sqrt(dx * dx + dy * dy + dz * dz);
				p = int((r - Rm) * PS);
				if (p < 0 || p >= NPot) continue;
				//printf("p=%d\n", p);
				if (r <= MAR[i1][3] && r <= MAR[i2][3]) U += Pot[p];
				else U += RepP[p];
			}
			//printf("l1\n");
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
				//printf("i1=%d, i2=%d\n", i1, i2);
				dx = PP2->X - PP1->X;
				dy = PP2->Y - PP1->Y;
				dz = PP2->Z - PP1->Z;
				//printf("tx1=%f, tx2=%f, ty1=%f, ty2=%f, tz1=%f, tz2=%f\n", 
					//   tx[i1], tx[i2], ty[i1], ty[i2], tz[i1], tz[i2]);
				r = sqrt(dx * dx + dy * dy + dz * dz);
				p = int((r - Rm) * PS);
				if (p < 0 || p >= NPot) continue;
				//printf("p=%d, r=%f, Rm=%f, PS=%f\n", p, r, Rm, PS);
				if (r <= MAR[i1][3] && r <= MAR[i2][3]) U += Pot[p];
				else U += RepP[p];
			}
		}
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
	printf("Calculation::initialize\n");
	int n, x, y, z, b;
	double rx, rz, rys = sqrt(0.5) * MaxY / double(PYS), rxs = MaxX / double(PXS), rzs = MaxZ / double(PZS);
	double y1 = 0.5 * MaxY - rys, R, dX, dZ, XMid = 0.5 * MaxX, ZMid = 0.5 * MaxZ;
	double y2 = 0.5 * MaxY + rys;
	double XF = double(XS) / MaxX, YF = double(YS) / MaxY, ZF = double(ZS) / MaxZ;
	for (x=0; x < XS; x++) for (y=0; y < YS; y++) for (z=0; z < ZS; z++) G[x][y][z] = 0;
	for (n=z=0, rz = 0.5 * rzs; z < PZS; z++, rz += rzs) 
	{
		for (x=0, rx = 0.5 * rxs; x < PXS; x++, n++, rx += rxs)
		{
			P[n].vX = P[n].vY = P[n].vZ = 0.0;
			P[n].X = rx;
			P[n].Y = y1;
			P[n].Z = rz;
			P[n].xp = ((b = int(XF * rx)) >= 0 ? (b < XS ? b : XS - 1) : 0);
			P[n].yp = ((b = int(YF * y1)) >= 0 ? (b < YS ? b : YS - 1) : 0);
			P[n].zp = ((b = int(ZF * rz)) >= 0 ? (b < ZS ? b : ZS - 1) : 0);
			P[n].prev = 0;
			P[n].next = G[P[n].xp][P[n].yp][P[n].zp];
			G[P[n].xp][P[n].yp][P[n].zp] = D[n] = P + n;
			if (P[n].next != 0) P[n].next->prev = D[n];
			if (z == 0 || x == 0 || z == PZS - 1 || x == PXS - 1) 
			{
				Fixed[n] = true;
				dX = P[x].X - XMid;
				dZ = P[x].Z - ZMid;
				R = 1.0 / sqrt(dX * dX + dZ * dZ);
				P[x].vX = dX * R;
				P[x].vZ = dZ * R;
			}
			else Fixed[n] = false;
			n++;
			P[n].vX = P[n].vY = P[n].vZ = 0.0;
			P[n].X = rx;
			P[n].Y = y2;
			P[n].Z = rz;
			P[n].xp = ((b = int(XF * rx)) >= 0 ? (b < XS ? b : XS - 1) : 0);
			P[n].yp = ((b = int(YF * y2)) >= 0 ? (b < YS ? b : YS - 1) : 0);
			P[n].zp = ((b = int(ZF * rz)) >= 0 ? (b < ZS ? b : ZS - 1) : 0);
			P[n].prev = 0;
			P[n].next = G[P[n].xp][P[n].yp][P[n].zp];
			G[P[n].xp][P[n].yp][P[n].zp] = D[n] = P + n;
			if (P[n].next != 0) P[n].next->prev = D[n];
			if (z == 0 || x == 0 || z == PZS - 1 || x == PXS - 1) 
			{
				Fixed[n] = true;
				dX = P[x].X - XMid;
				dZ = P[x].Z - ZMid;
				R = 1.0 / sqrt(dX * dX + dZ * dZ);
				P[x].vX = dX * R;
				P[x].vZ = dZ * R;
			}
			else Fixed[n] = false;
		}
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
	int n, i, x, y, z, m;
	double hh = 0.5 * h, h6 = h / 6.0, *ax = new double[N], *ay = new double[N];
	double *az = new double[N], *dxm = new double[N], *dym = new double[N], *dzm = new double[N];
	double *dvxm = new double[N], *dvym = new double[N], *dvzm = new double[N];
	double *dxt = new double[N], *dyt = new double[N], *dzt = new double[N];
	double *dvxt = new double[N], *dvyt = new double[N], *dvzt = new double[N];
	double *xt = new double[N], *yt = new double[N], *zt = new double[N], vF, R, lR, dX, dZ;
	double XF = double(XS) / MaxX, YF = double(YS) / MaxY, ZF = double(ZS) / MaxZ, MA, ZMid = 0.5 * MaxZ, Pi2 = 2.0 * M_PI;
	double XMMiB = Re, XMMaB = MaxX - Re, ZMMiB = Re, ZMMaB = MaxZ - Re, *Angle = new double[N], XMid = 0.5 * MaxX, CA;
	Particle *PB;
	for (n=0; n<N; n++) if (Fixed[n]) Angle[n] = tan((XMid - P[n].X) / (ZMid - P[n].Z));
	for (i=0, Run = true; Run; i++)
	{
		//printf("i=%d, ", i);
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
			/*for (n=0; n<N; n++) if (Fixed[n] && P[n].Y > YMid && P[n].X > XMMiB && P[n].X < XMMaB 
										&& P[n].Z > ZMMiB && P[n].Z < ZMMaB)
			{
				for (i=0, MA = Pi2; i<N; i++) if (Fixed[i] && P[n].Y < YMid)
				{
					if ((CA = fabs(Angle[i] - Angle[n])) < MA)
					{
						MA = CA;
						m = i;
					}
					else if ((CA -= Pi2) < MA)
					{
						MA = CA;
						m = i;
					}
				}
				P[m].Y = P[n].Y;
				P[m].vX = P[n].vX;
				P[m].vZ = P[n].vZ;
				dX = P[n].X - XMid;
				dZ = P[n].Z - ZMid;
				lR = (R = sqrt(dX * dX + dZ * dZ)) + Re;
				P[m].X = XMid - lR * P[m].vX;
				P[m].Z = ZMid - lR * P[m].vZ;
				Fixed[n] = false;
				lR = Speed / R;
				P[n].vX *= lR;
				P[n].vY = 0.0;
				P[n].vZ *= lR;
				for (PB = G[P[m].xp][P[m].yp][P[m].zp], R = MaxX; PB != 0; PB = PB->next) 
					if (PB->Y < YMid && !Fixed[i = P - PB]) 
				{
					dX = PB->X - P[m].X;
					dZ = PB->Z - P[m].Z;
					if ((lR = sqrt(dX * dX + dZ * dZ)) < R)
					{
						R = lR;
						x=i;
					}
				}
				if (R < MaxX)
				{
					Fixed[x] = true;
					dX = P[x].X - XMid;
					dZ = P[x].Z - ZMid;
					R = 1.0 / sqrt(dX * dX + dZ * dZ);
					P[x].vX = dX * R;
					P[x].vZ = dZ * R;
				}
			}*/
		}
		for (n=0, U = T = 0.0; n<N; n++) 
			T += P[n].vX * P[n].vX + P[n].vY * P[n].vY + P[n].vZ * P[n].vZ;
		geta(xt, yt, zt, ax, ay, az);
		if (E == 0.0 || T == 0.0 || Move) E = T + U;
		/*else if (T + U != E && E > U)
		{
			vF = sqrt((E - U) / T);
			for (n=0; n<N; n++)
			{
				P[n].vX *= vF;
				P[n].vY *= vF;
				P[n].vZ *= vF;
			}
		}*/
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
		geta(xt, yt, zt, dvxt, dvyt, dvzt);
		for (n=0; n<N; n++) if (!Fixed[n])
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
		for (n=0; n<N; n++) if (!Fixed[n])
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
		for (n=0; n<N; n++) if (!Fixed[n])
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
