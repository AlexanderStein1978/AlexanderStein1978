//
// C++ Implementation: potentialplot
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2008 - 2020
//
// Copyright: See README file that comes with this source code
//
//


#include "potentialplot.h"
#include "potential.h"
#include "MainWindow.h"
#include "naturalspline.h"
#include "SplinePoint.h"

#include <stdio.h>

#include <QPainter>
#include <QMenu>
#include <QToolButton>


PotentialPlot::PotentialPlot(Potential *potential, MainWindow *MW) : DiagWindow(MDIChild::SimpleDiagWindow, MW)
{
	//printf("PotentialPlot::PotentialPlot\n");
	pot = 0;
	History = true;
	showDiagFuncs = false;
	WFS = WWFS = 0;
	nDFP = 0;
	copyP = 0;
	xStartLabel->setText("R [A]:    from ");
    xStart->setText("0");
    xStopLabel->setText(" to ");
    xStop->setText("100");
    yStartLabel->setText("Energie [cm-1]  :    from ");
    yStart->setText("0");
    yStopLabel->setText(" to ");
    yStop->setText("5000");
	NCopies = copyP = 0;
	potCopies = new Potential*[26];
	plotPotential(potential);
	setWindowTitle(pot != 0 ? "Plot of potential " + pot->getName() : "Potential plot");
	pictureMenu = new QMenu(this);
	addPointAct = new QAction("Add point", this);
	addPointAct->setStatusTip("Add a new point to the spline potential");
	connect(addPointAct, SIGNAL(triggered()), this, SLOT(addPoint()));
	pictureMenu->addAction(addPointAct);
	pointMenu = new QMenu(this);
	removePointAct = new QAction("Remove point", this);
	removePointAct->setStatusTip("Remove the selected point from the spline potential");
	connect(removePointAct, SIGNAL(triggered()), this, SLOT(removePoint()));
	pointMenu->addAction(removePointAct);
	//printf("Ende PotentialPlot::PotentialPlot\n");
}

PotentialPlot::~PotentialPlot()
{
	printf("PotentialPlot::~PotentialPlot\n");
	int n;
	for (n=0; n < NCopies; n++) delete potCopies[n];
	delete[] potCopies;
	delete[] CopyColor;
	if (points != 0) delete[] points;
	if (WFS != 0)
	{
		delete[] WFS;
		delete[] WWFS;
	}
}

void PotentialPlot::addPoint()
{
	if (pot != 0) pot->addPoint((double(cPosx) - XO) / XSF, (YO - double(cPosy)) / YSF);
    cPosx = cPosy = 0;
}

void PotentialPlot::addPotential(Potential *addPot)
{
	printf("PotentialPlot::addPotential\n");
	if (NCopies < 26) potCopies[NCopies++] = new Potential(*addPot);
	else
	{
		if (copyP == 26) copyP = 0;
		delete potCopies[copyP];
		potCopies[copyP++] = new Potential(*addPot);
	}
	SetBoundaries();
	Paint();
}

void PotentialPlot::removePotential(Potential *remPot)
{
    if (nullptr == remPot) return;
    if (pot == remPot) plotPotential();
    else
    {
        int n;
        for (n=0; n < NCopies && potCopies[n]->getFileName() != remPot->getFileName(); ++n) ;
        if (n < NCopies)
        {
            delete potCopies[n];
            for (int m=n+1; m < NCopies; ++m) potCopies[m-1] = potCopies[m];
            --NCopies;
            copyP ^= copyP;
        }
    }
}

void PotentialPlot::clearHistory()
{
	printf("PotentialPlot::clearHistory()\n");
	int n;
	if (NCopies == 26)
	{
		if ((--copyP) == -1) copyP = 25;
		for (n=0; n < 26; n++) if (n != copyP) delete potCopies[n];
		potCopies[0] = potCopies[copyP];
		for (n=1; n < 26; n++) potCopies[n] = 0;
		NCopies = 1;
		copyP = 0;
	}
	else
	{
		for (n=0; n < NCopies - 1; n++) delete potCopies[n];
		potCopies[0] = potCopies[NCopies - 1];
		for (n=1; n < NCopies; n++) potCopies[n] = 0;
		NCopies = 1;
	}
	Paint();
}

void PotentialPlot::closeEvent(QCloseEvent *E)
{
    emit closing();
    E->accept();
}

bool PotentialPlot::getShowDiagFuncs()
{
	return showDiagFuncs;
}

bool PotentialPlot::getShowHistory()
{
	return History;
}

void PotentialPlot::HandleHistoryAfterMoving()
{
    int m, n = (copyP == 0 ? NCopies - 1: copyP - 1);
    delete potCopies[n];
    if (History) potCopies[n] = new Potential(*pot);
    else
    {
        potCopies[n] = 0;
        NCopies--;
        if (copyP != 0)
        {
            Potential **Buff = new Potential*[26];
            for (n = copyP, m=0; m < 25; n++, m++)
            {
                if (n==26) n=0;
                Buff[m] = potCopies[n];
            }
            delete[] potCopies;
            potCopies = Buff;
            copyP = 0;
        }
    }
}

void PotentialPlot::HandleHistoryWhileMoving()
{
    if (!History)
    {
        if (NCopies < 26) potCopies[NCopies++] = new Potential(*pot);
        else
        {
            printf("History, vor delete\n");
            if (copyP == 26) copyP = 0;
            delete potCopies[copyP];
            potCopies[copyP++] = new Potential(*pot);
        }
    }
}

void PotentialPlot::MovePoint()
{
    potCopies[(copyP == 0 ? NCopies - 1 : copyP - 1)]->movePoint(mPoint,
                                                        points[sPoint].x, points[sPoint].y);
}

void PotentialPlot::MovePoint(int i_pointIndex, double i_newX, double i_newY)
{
    pot->movePoint(i_pointIndex, i_newX, i_newY);
}

void PotentialPlot::plotPotential(Potential *potential)
{
	//printf("PotentialPlot::plotPotential, Name=%s\n", potential->getName().toAscii().data());
    if (pot != 0)
    {
        disconnect(pot, SIGNAL(propertiesChanged()), this, SLOT(PotentialChanged()));
        addPotential(pot);
    }
	pot = potential;
	if (pot != 0) 
	{
		connect(pot, SIGNAL(propertiesChanged()), this, SLOT(PotentialChanged()));
		if (showDiagFuncs) potential->setCalcDiagFuncs(true);
	}
	XUnit = "R [A]";
	YUnit = "Energy [cm^{-1}]";
	PotentialChanged();
	ResetZoom();
}

void PotentialPlot::PotentialChanged()
{
	//printf("PotentialPlot::PotentialChanged\n");
	if (pot == 0) return;
	if (History)
	{
		//printf("History, NCopies=%d, copyP=%d\n", NCopies, copyP);
		if (NCopies < 26) potCopies[NCopies++] = new Potential(*pot);
		else
		{
			//printf("History, vor delete\n");
			if (copyP == 26) copyP = 0;
			delete potCopies[copyP];
			//printf("Nach delete\n");
			potCopies[copyP++] = new Potential(*pot);
		}
	}
	else
	{
		if (NCopies > 0) delete potCopies[0];
        else ++NCopies;
		potCopies[0] = new Potential(*pot);
	}
	if (showPoints)
	{
		//printf("Vor delete points\n");
		if (points != 0) delete[] points;
		//printf("Nach delete points\n");
		points = pot->getSplinePoints(numPoints);
	}
	if (showDiagFuncs)
	{
		if (WFS != 0)
		{
			delete[] WFS;
			delete[] WWFS;
		}
        pot->getDiagFuncs(WFS, WWFS, DFRi, DFh, nDFP, NumPoints);
	}
	SetBoundaries();
	Paint();
	//printf("Ende PotentialChanged\n");
}

void PotentialPlot::PotSnapShot()
{
	printf("PotentialPlot::PotSnapShot\n");
	if (History) setShowHistory(false);
	if (NCopies < 26) potCopies[NCopies++] = new Potential(*pot);
	else
	{
		if (copyP == 26) copyP = 0;
		delete potCopies[copyP];
		potCopies[copyP++] = new Potential(*pot);
	}
	Paint();
}

void PotentialPlot::PSpektrum(QPainter &P, const QRect & A, bool /*PrintFN*/ )
{
	//printf("PotentialPlot::PSpektrum, NCopies=%d, copyP=%d, potName=%s\n", 
		//   NCopies, copyP, potCopies[0]->getName().ascii());
	if (pot == 0) 
	{
		printf("Pot = 0!!\n");
		return;
	}
	int i, s, w = A.width() - ScaleYWidth, h = A.height() - ScaleXHeight;
	if (w <= 0 || h <= 0) return;
	int l = A.left() + ScaleYWidth - 1;
	int t = A.top(), r = l + w, lp, n;
	int p=0, b = A.bottom() - ScaleXHeight, NFC, FC;
    double minR = xStart->text().toDouble(), maxR = xStop->text().toDouble();
    double minE = yStart->text().toDouble(), maxE = yStop->text().toDouble();
    double ESc = (double)h / (maxE - minE);
	double *Data;
	//printf("Vor 1. Schleife\n");
	if (showDiagFuncs && WFS != 0 && WWFS != 0)
	{
		double R, S, S2, M, m = 0.0, M2 = 0.0;
		int ax = 0, ay = 0, lx, ly, ay2 = 0, ly2, O2;
		P.setPen(QColor(150, 150, 150));
		for (i=0, M = 0.0; i < nDFP; i++) 
		{
			if (WFS[i] > M) M = WFS[i];
			if (WWFS[i] < m) m = WWFS[i];
			else if (WWFS[i] > M2) M2 = WWFS[i];
		}
		if (DFRi < minR) s = int((minR - DFRi) / DFh) + 1;
		else s=0;
		S2 = double(2*h) / (3.0 * (M2 - m));
		O2 = t + M2 * S2;
		for (i=s, R = DFRi + s * DFh, S = double(h) / (3.0 * M); i < nDFP && R < maxR; i++, R += DFh)
		{
			lx = ax;
			ly = ay;
			ly2 = ay2;
			ax = XO + XSF * R;
			ay = b - S * WFS[i];
			ay2 = O2 - S2 * WWFS[i];
			if (lx != 0 && ly != 0) 
			{
				P.drawLine(lx, ly, ax, ay);
				P.drawLine(lx, ly2, ax, ay2);
			}
		}
		P.drawLine(l, O2, r, O2);
	}
	//printf("Vor 2. Schleife\n");
	for (n = (NCopies < 26 ? 0 : copyP); true; n++)
	{
		if (n == NCopies) if ((n=0) == copyP) break;
		P.setPen(CopyColor[n]);
		for (FC = 0, NFC = potCopies[n]->getNFC(); FC < NFC; FC++)
		{
			Data = potCopies[n]->getPoints(minR, maxR, w+1, FC);
			if (Data == 0) 
			{
				printf("Error, no potential data available for copy %d\n", n);
				continue;
			}
			for (s=0; (s<=w ? Data[s] == 0 : false); s++) ;
			for (i=s+l; i<=l+w; i++)
			{
				lp = p;
				p = b - int(ESc * (Data[i-l] - minE)) + 1;
				if (p > b) p = b + 1;
				if (i>l+s ? (p > lp ? p - lp : lp - p) > 1 : false) P.drawLine(i-1, lp, i, p);
				else P.drawPoint(i, p);
			}
			delete[] Data;
		}
		if (n == copyP - 1) break;
	}
	//printf("Vor 3. Schleife\n");
    DrawPoints(P, A);
	//printf("Ende PSpektrum\n");
}

void PotentialPlot::removePoint()
{
	if (pot != 0 && sPoint >= 0) pot->removePoint(sPoint);
	sPoint = -1;
}

void PotentialPlot::ResetZoom()
{
	xStart->setText(QString::number(XMin, 'g', 5));
    xStop->setText(QString::number(XMax, 'g', 5));
    yStart->setText(QString::number(YMin, 'g', 5));
    yStop->setText(QString::number(YMax, 'g', 5));
	printf("Vor Ende, XMin=%g, XMax=%g, YMin=%g, YMax=%g\n", XMin, XMax, YMin, YMax);
    Paint();
}

void PotentialPlot::SetBoundaries()
{
	printf("PotentialPlot::setBoundaries()\n");
	if (pot == 0) return;
	double R, D, Ym, YM, YMa;
	int n;
	XMin = pot->getMinR();
	XMax = pot->getMaxR();
	printf("XMin=%f, XMax=%f\n", XMin, XMax);
	for (n=0, YMin = 1e99, YMax = -1e99; n < NCopies; n++)
	{
		//printf("XMin=%f, XMax=%f\n", XMin, XMax);
		potCopies[n]->getMinimum(R, Ym);
		//printf("Minimum R=%f, E=%f\n", R, Ym);
		potCopies[n]->getMaximum(R, YM);
		YMa = 1.5 * potCopies[n]->getPoint(XMax) - 0.5 * Ym;
		if (YM > YMa) YM = YMa;
		//printf("Maximum R=%f, E=%f\n", R, YM);
		if (Ym < YMin) YMin = Ym;
		if (YM > YMax) YMax = YM;
	}
	XMin -= 0.1 * (XMax - XMin);
	if (XMin < 0.0) XMin = 0.0;
	YMin -= (D = 0.1 * (YMax - YMin));
	YMax += D;
}

void PotentialPlot::SetPoints()
{
    DiagWindow::SetPoints();
    if (pot != 0) points = pot->getSplinePoints(numPoints);
}

void PotentialPlot::setShowDiagFuncs(bool show)
{
	if (pot != 0) 
	{
		pot->setCalcDiagFuncs(show);
		showDiagFuncs = show;
		if (show) 
		{
			if (WFS != 0)
			{
				delete[] WFS;
				delete[] WWFS;
			}
            pot->getDiagFuncs(WFS, WWFS, DFRi, DFh, nDFP, NumPoints);
		}
		Paint();
	}
}

void PotentialPlot::setShowHistory(bool H)
{
	if (!H && History && NCopies > 1) clearHistory();
	History = H;
}

void PotentialPlot::ShowPopupMenu(const QPoint& i_point)
{
    if (ZoomB->isChecked()) return;
    if (sPoint != -1) pointMenu->popup(i_point);
    else pictureMenu->popup(i_point);
}
