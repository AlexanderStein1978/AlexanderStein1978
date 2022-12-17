//
// C++ Implementation: residualplot
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2009 - 2020
//
// Copyright: See README file that comes with this source code
//
//


#include "residualplot.h"
#include "MainWindow.h"
#include "termtable.h"
#include "linetable.h"
#include "elstate.h"
#include "molecule.h"
#include "utils.h"
#include "fitdata.h"
#include "fit.h"
#include "tableline.h"
#include "PerturbationInfoTable.h"
#include "isotab.h"
#include "Js.h"
#include "point.h"
#include "LocalPerturbation.h"
#include "progression.h"
#include "line.h"

#include <stdio.h>

#include <QComboBox>
#include <QPainter>
#include <QPaintDevice>
#include <QPrintDialog>
#include <QMenu>
#include <QToolButton>
#include <QGridLayout>


ResidualPlot::ResidualPlot(MainWindow *MW) : DiagWindow(ResidPlot, MW),  m_clickedSpline(-1), m_clickedPoint(-1), showResidualFit(false),
    m_ResidualFit(0), m_PerturbationInfoTable(0)
{
	//printf("ResidualPlot::ResidualPlot\n");
	Block = false;
	lineT = 0;
	mol = 0;
	refT = 0;
	state = 0;
	termT = 0;
	FitD = 0;
	v = C = Iso = -1;
	refTIsoT = 0;
	refTcT = 0;
	MrefTcT = 0;
	FSA = false;
	initData();
	setWindowTitle("Residual plot");
	xStartLabel->setText("Min J:");
	xStopLabel->setText("Max J:");
	yStartLabel->setText("Min dev.:");
	yStopLabel->setText("Max dev.:");
	XUnit = "J";
	YUnit = "Residuals [cm^{-1}]";
	MolBox->setEditable(false);
	StateBox->setEditable(false);
	SourceBox->setEditable(false);
	RefBox->setEditable(false);
	IsoBox->setEditable(false);
	vBox->setEditable(false);
	CompBox->setEditable(false);
	connect(MolBox, SIGNAL(currentIndexChanged(int)), this, SLOT(molBoxChanged()));
	connect(StateBox, SIGNAL(currentIndexChanged(int)), this, SLOT(stateBoxChanged(int)));
	connect(SourceBox, SIGNAL(currentIndexChanged(int)), this, SLOT(sourceBoxChanged()));
	connect(RefBox, SIGNAL(currentIndexChanged(int)), this, SLOT(refSourceBoxChanged()));
	connect(IsoBox, SIGNAL(currentIndexChanged(int)), this, SLOT(isoBoxChanged(int)));
	connect(vBox, SIGNAL(currentIndexChanged(int)), this, SLOT(vBoxChanged(int)));
	connect(CompBox, SIGNAL(currentIndexChanged(int)), this, SLOT(compBoxChanged(int)));
	connect(PrintB, SIGNAL(clicked()), this, SLOT(Printall()));
	connect(Bild, SIGNAL(RightClicked(QPoint*)), this, SLOT(showPMenu(QPoint*)));
	connect(Bild, SIGNAL(SelectionChanged(QRect*, bool)), this, SLOT(marked(QRect*, bool)));
	PMenu = new QMenu(this);
	markUncAct = new QAction("Select with &err...", this);
	markUncAct->setStatusTip("Select the levels with the specified uncertainty...");
	PMenu->addAction(markUncAct);
	markCompAct = new QAction("Select of &comp...", this);
	markCompAct->setStatusTip("Select the levels of the specified (fine structure) component...");
	PMenu->addAction(markCompAct);
	markLinesAct = new QAction("&Mark levels", this);
	markLinesAct->setStatusTip("Show the levels in the source table");
	PMenu->addAction(markLinesAct);
	setErrorAct = new QAction("&Set uncertainty...", this);
	setErrorAct->setStatusTip("Enter an uncertainty to set for the selected levels");
	PMenu->addAction(setErrorAct);
	add9toErrorAct = new QAction("&Add 9 to uncertainty", this);
	add9toErrorAct->setStatusTip("Add the value of 9 cm^-1 to the uncertainty of the selected levels");
	PMenu->addAction(add9toErrorAct);
	sub9fErrorAct = new QAction("S&ubtract 9 from uncertainty", this);
	sub9fErrorAct->setStatusTip(
		"Subtract the value of 9 cm^-1 from the uncertainties of the selected levels with uncertainties larger this value");
	PMenu->addAction(sub9fErrorAct);
	deleteLevelsAct = new QAction("&Delete levels", this);
	deleteLevelsAct->setStatusTip("Deletes the selected levels from their data set");
	PMenu->addAction(deleteLevelsAct);
    residualFitSeparatorAct = PMenu->addSeparator();
    residualFitSeparatorAct->setVisible(false);
    showResidualFitAct = new QAction("Show &residual fit", this);
    showResidualFitAct->setStatusTip(
        "Show/Hide the available spline curves of the residual fit which can be used to automatize the update of the level assignment of an added up absorption spectrum.");
    showResidualFitAct->setVisible(false);
    PMenu->addAction(showResidualFitAct);
    showResidualFitSplinePointsAct = new QAction("Show spline points", this);
    showResidualFitSplinePointsAct->setStatusTip("Show the spline points of the residual fit to be able to move it, remove it or add more of them.");
    showResidualFitSplinePointsAct->setVisible(false);
    PMenu->addAction(showResidualFitSplinePointsAct);
    addResidualFitSplinePointAct = new QAction("Add spline point", this);
    addResidualFitSplinePointAct->setStatusTip("Add a spline point to the residual fit at the clicked position.");
    addResidualFitSplinePointAct->setVisible(false);
    PMenu->addAction(addResidualFitSplinePointAct);
    deleteResidualFitSplinePointAct = new QAction("Delete spline point", this);
    deleteResidualFitSplinePointAct->setStatusTip("Remove the marked spline point from the residual fit.");
    deleteResidualFitSplinePointAct->setVisible(false);
    PMenu->addAction(deleteResidualFitSplinePointAct);
    joinResidualFitSplinesAct = new QAction("&Join two splines", this);
    joinResidualFitSplinesAct->setStatusTip("Make one spline out of the first two neighboring splines data points are marked of or the first two splines if no points are marked.");
    joinResidualFitSplinesAct->setVisible(false);
    PMenu->addAction(joinResidualFitSplinesAct);
    createResidualFitAct = new QAction("Create residual &fit", this);
    createResidualFitAct->setStatusTip("Create a new residual fit which can be used to automatize the update of the level assignment of an added up absorption spectrum.");
    createResidualFitAct->setVisible(false);
    PMenu->addAction(createResidualFitAct);
    addLocalPerturbationAct = new QAction("Add local &perturbation", this);
    addLocalPerturbationAct->setStatusTip("Fit a local perturbation to the selected residuals or if not more than one level is selected to the whole band.");
    addLocalPerturbationAct->setVisible(false);
    PMenu->addAction(addLocalPerturbationAct);
    removeLocalPerturbationsAct = new QAction("Remove &local perturbations", this);
    removeLocalPerturbationsAct->setStatusTip("Remove the local perturbations in the selected levels or if no levels are selected all local perturbations in the whole band.");
    removeLocalPerturbationsAct->setVisible(false);
    PMenu->addAction(removeLocalPerturbationsAct);
    showLocalPerturbationInfoTableAct = new QAction("Show local perturbation &info table", this);
    showLocalPerturbationInfoTableAct->setStatusTip("Show information table with most important constants of the local perturbations.");
    showLocalPerturbationInfoTableAct->setVisible(false);
    PMenu->addAction(showLocalPerturbationInfoTableAct);
	connect(markUncAct, SIGNAL(triggered()), this, SLOT(markUnc()));
	connect(markCompAct, SIGNAL(triggered()), this, SLOT(markComp()));
	connect(markLinesAct, SIGNAL(triggered()), this, SLOT(markLines()));
	connect(setErrorAct, SIGNAL(triggered()), this, SLOT(setError()));
	connect(add9toErrorAct, SIGNAL(triggered()), this, SLOT(add9toError()));
	connect(sub9fErrorAct, SIGNAL(triggered()), this, SLOT(sub9fromError()));
	connect(deleteLevelsAct, SIGNAL(triggered()), this, SLOT(deleteLevels()));
    connect(showResidualFitAct, SIGNAL(triggered()), this, SLOT(ShowResidualFit()));
    connect(createResidualFitAct, SIGNAL(triggered()), this, SLOT(CalcResidualFit()));
    connect(addLocalPerturbationAct, SIGNAL(triggered()), this, SLOT(AddLocalPerturbation()));
    connect(removeLocalPerturbationsAct, SIGNAL(triggered()), this, SLOT(RemoveLocalPerturbations()));
    connect(showLocalPerturbationInfoTableAct, SIGNAL(triggered()), this, SLOT(ShowLocalPerturbationInfoTable()));
    connect(addResidualFitSplinePointAct, SIGNAL(triggered()), this, SLOT(addPoint()));
    connect(deleteResidualFitSplinePointAct, SIGNAL(triggered()), this, SLOT(removePoint()));
    connect(joinResidualFitSplinesAct, SIGNAL(triggered()), this, SLOT(JoinResidualFitSplines()));
    connect(showResidualFitSplinePointsAct, SIGNAL(QAction::triggered()), this, SLOT(SwitchShowPoints()));
    

#ifdef TestResidualFit
	SPlot = MW->CreateDiagWindow();
	SPlot->show();
#endif
	
	moleculesChanged();
	//printf("Ende ResidualPlot\n");
}

ResidualPlot::~ResidualPlot()
{
	destroyData();
	if (refTIsoT != 0) delete[] refTIsoT;
	if (refTcT != 0) delete[] refTcT;
	/*if (mol != 0) disconnect(mol, SIGNAL(currentIndexChanged(int)), this, SLOT(molBoxChanged()));
	if (lineT != 0) 
		disconnect(lineT, SIGNAL(currentIndexChanged(int)), this, SLOT(sourceBoxChanged()));
	if (termT != 0)
		disconnect(termT, SIGNAL(currentIndexChanged(int)), this, SLOT(sourceBoxChanged()));
	if (refT != 0) 
		disconnect(refT, SIGNAL(currentIndexChanged(int)), this, SLOT(refSourceBoxChanged()));*/
}

void ResidualPlot::add9toError()
{
	if (lineT == 0 && FitD == 0) return;
	int n, N;
	for (n=N=0; n < nData[C][Iso][v]; n++) if (Marked[n]) N++;
	int *rows = new int[N];
	double *u = new double[N];
	for (n=N=0; n < nData[C][Iso][v]; n++) if (Marked[n])
	{
		rows[N] = Row[C][Iso][v][n];
		u[N++] = Unc[C][Iso][v][n] + 9.0;
	}
	if (lineT != 0) lineT->setUncertainty(rows, u, N);
	else FitD->setUncertainty(rows, u, N);
	delete[] rows;
	delete[] u;
}

void ResidualPlot::AddLocalPerturbation()
{
    if (0 == nData[C][Iso][v] || 0 == state || 0 == mol) return;
    IsoTab *IsoT = mol->getIso();
    if (0 == IsoT) return;
    double Be = state->getBe();
	if (Be <= 0.0)
    {
        QMessageBox::information(this, "MolSpektAnalysis", "No valid rotational constant Be for the current electronic state available. Without this no initial guess neccessary for the fit is possible!");
        return;
    }
    if (!m_ResidualFit->isFitDataAvailable()) SetResidualFitData();
    std::vector<Js> pertFitData = GetPertFitData();
    int n = m_ResidualFit->addLocalPerturbation(IsoT->relRedMass[Iso], state->getOmega(), pertFitData, mol->getJStep(Iso), Be);
    delete IsoT;
    Paint();
    if (n>=0)
    {
        ShowLocalPerturbationInfoTable();
        m_PerturbationInfoTable->SelectPerturbation(n);
    }
}

void ResidualPlot::JoinResidualFitSplines()
{
    if (!m_ResidualFit->isFitDataAvailable()) SetResidualFitData();
    std::vector<Js> pertFitData = GetPertFitData();
    if (!pertFitData.empty()) m_ResidualFit->JoinSplines(-1, pertFitData);
    Paint();
}

std::vector<Js> ResidualPlot::GetPertFitData() const
{
    if (0 == nData[C][Iso][v]) return std::vector<Js>();
    int *SortArray = new int[nData[C][Iso][v]], B=0, n, JStart, NJ = 0;
    for (n=0; n < nData[C][Iso][v]; ++n) SortArray[n] = n;
    while (B>=0) for (n=1, B=-1; n < nData[C][Iso][v]; ++n) if (JData[C][Iso][v][SortArray[n-1]] > JData[C][Iso][v][SortArray[n]])
    {
        B = SortArray[n-1];
        SortArray[n-1] = SortArray[n];
        SortArray[n] = B;
    }
    for (n=0, B=-1; n < nData[C][Iso][v]; ++n)
    {
        if (Marked[SortArray[n]])
        {
            if (B == -1) B = n;
            if (n - B + 1 > NJ)
            {
                NJ = n - B + 1;
                JStart = B;
            }
        }
        else B=-1;
    }
    if (NJ <= 1)
    {
        JStart = 0;
        NJ = nData[C][Iso][v];
        if (NJ <= 1) return std::vector<Js>();
    }
    std::vector<Js> pertFitData;
    for (n = JStart; n < JStart + NJ; ++n)
        pertFitData.push_back(Js(JData[C][Iso][v][SortArray[n]], Obs[C][Iso][v][SortArray[n]] - Data[C][Iso][v][SortArray[n]], Obs[C][Iso][v][SortArray[n]], Unc[C][Iso][v][SortArray[n]]));
    delete[] SortArray;
    return pertFitData;
}

void ResidualPlot::ShowLocalPerturbationInfoTable()
{
    if (0 == m_PerturbationInfoTable)
    {
        m_PerturbationInfoTable = new PerturbationInfoTable;
        MW->showMDIChild(m_PerturbationInfoTable);
    }
    if (!m_PerturbationInfoTable->isVisible()) m_PerturbationInfoTable->show();
    m_PerturbationInfoTable->SetResidualFit(m_ResidualFit);
}

void ResidualPlot::ClearMarked()
{
    if (Marked == 0) return;
	int n;
	for (n=0; n < nData[C][Iso][v]; n++) Marked[n] = false;
	Paint();
}

void ResidualPlot::compBoxChanged(int i)
{
	if (Block) return;
	if (Marked != 0) delete[] Marked;
    if (Data == 0 || Iso < 0)
	{
		Paint();
		Marked = 0;
		return;
	}
	if (FSA) i = CompBox->currentText().toInt();
	else if (i==0 && nData[0][refTIsoT[Iso]][v] == 0 && NC > 1) i=1;
	int J, XM;
	double mR = 1000.0, MR = 0.0, b, mJ, MJ;
	if (nData[C=i][refTIsoT[Iso]][v] > 0) Marked = new bool[nData[C][refTIsoT[Iso]][v]];
	else 
	{
		Paint();
		Marked = 0;
		return;
	}
	for (J=0; J < nData[C][refTIsoT[Iso]][v]; J++)
	{
		if (Data[C][refTIsoT[Iso]][v][J] > MR) MR = Data[C][refTIsoT[Iso]][v][J];
		if (Data[C][refTIsoT[Iso]][v][J] < mR) mR = Data[C][refTIsoT[Iso]][v][J];
		Marked[J] = false;
	}
	CompList.clear();
	if (Comp != 0) for (J=0; J < nData[C][refTIsoT[Iso]][v]; J++) if (!CompList.contains(Comp[C][refTIsoT[Iso]][v][J]))
		CompList.append(Comp[C][refTIsoT[Iso]][v][J]);
	if (CompList.count() > 1) markCompAct->setEnabled(true);
	else markCompAct->setEnabled(false);
	mJ = JData[C][refTIsoT[Iso]][v][0];
	MJ = JData[C][refTIsoT[Iso]][v][nData[C][refTIsoT[Iso]][v] - 1];
	b = 0.1 * (MJ - mJ);
	XMin = mJ - b;
	XMax = MJ + 5.0 * b;
	XM = int(MJ + b);
	xStart->setText(QString::number(int(XMin) - 1));
	xStop->setText(QString::number(XM + 1));
	b = 0.1 * (MR - mR);
	YMin = mR - b;
	YMax = MR + 5.0 * b;
	yStart->setText(QString::number(YMin, 'g', 5));
	yStop->setText(QString::number(MR + 2 * b, 'g', 5));
    if (0 != m_ResidualFit) disconnect(m_ResidualFit, SIGNAL(Changed()), this, SLOT(Paint()));
    m_ResidualFit = (FitD != 0 ? FitD->getResidualFit(state, Iso, v, C) : 0);
    if (0 != m_ResidualFit) connect(m_ResidualFit, SIGNAL(Changed()), this, SLOT(Paint()));
    if (0 != m_PerturbationInfoTable && 0 != m_ResidualFit && m_ResidualFit->getNumberOfLocalPerturbations() > 0) m_PerturbationInfoTable->SetResidualFit(m_ResidualFit);
	Paint();
}

void ResidualPlot::CalcResidualFit()
{
    if (FitD != 0)
    {
        if (m_ResidualFit == 0)
        {
            m_ResidualFit = new ResidualFit(state);
            QString stateName(state != 0 ? state->getName() : "");
            m_ResidualFit->setAssignment(&stateName, Iso, v, C);
            FitD->setResidualFit(m_ResidualFit);
        }
        Point *points = getResidualFitData();
        m_ResidualFit->FitTo(nData[C][refTIsoT[Iso]][v], points);
        delete[] points;
        showResidualFit = true;
        Paint();

#ifdef TestResidualFit
        Spline *Sp = TestResFit.getSpline(0);
        int n, N = Sp->getNSFuncs();
        double **Data = Create(1000, 2), r, R0 = Sp->getx0(), RN = Sp->getxN();
        double h = (RN - R0) / 999;
        SPlot->clear();
        for (n=0, r = R0; n < 999; n++, r+=h) Data[n][0] = r;
        Data[n][0] = RN;
        for (n=0; n<N; n++)
        {
            Sp->GetSFuncs(n, Data, 1000);
            SPlot->addData(Data, 1000);
        }
        Destroy(Data, 1000);
#endif
    }
}

void ResidualPlot::deleteLevels()
{
	if (FitD == 0 && lineT == 0) return;
	int n, N;
	for (n=N=0; n < nData[C][Iso][v]; n++) if (Marked[n]) N++;
	int *rows = new int[N];
	for (n=N=0; n < nData[C][Iso][v]; n++) if (Marked[n]) rows[N++] = Row[C][Iso][v][n];
	if (lineT != 0) lineT->deleteRows(rows, N);
	else FitD->deleteRows(rows, N);
	delete[] rows;
}

void ResidualPlot::destroyData()
{
	int c, n, i;
	if (nData != 0) Destroy(nData, NC, NIso);
	if (Data != 0) 
	{
		for (c=0; c < NC; c++) 
		{
			for (i = 0; i < NIso; i++) 
			{
				for (n=0; n < Nv; n++) if (Data[c][i][n] != 0) 
				{
					delete[] Data[c][i][n];
					delete[] JData[c][i][n];
					delete[] Unc[c][i][n];
                    delete[] Obs[c][i][n];
					if (Row != 0) delete[] Row[c][i][n];
					if (Comp != 0) delete[] Comp[c][i][n];
				}
				delete[] Data[c][i];
				delete[] JData[c][i];
				delete[] Unc[c][i];
                delete[] Obs[c][i];
				if (Row != 0) delete[] Row[c][i];
				if (Comp != 0) delete[] Comp[c][i];
			}
			delete[] Data[c];
			delete[] JData[c];
			delete[] Unc[c];
            delete[] Obs[c];
			if (Row != 0) delete[] Row[c];
			if (Comp != 0) delete[] Comp[c];
		}
		delete[] Data;
		delete[] JData;
		delete[] Unc;
        delete[] Obs;
		if (Row != 0) delete[] Row;
		if (Comp != 0) delete[] Comp;
	}
	if (nv != 0) delete[] nv;
	if (Marked != 0) delete[] Marked;
	nData = 0;
	Data = 0;
	nv = 0;
	Marked = 0;
}

void ResidualPlot::drawDataPoint(QPainter *P, int x, int y, int u1, int u2, int e)
{
	if (u1==u2)
	{
		P->drawLine(x, y-2, x, y+2);
		P->drawLine(x-2, y, x+2, y);
	}
	else
	{
		P->drawRect(x-2, y-2, 4, 4);
		P->drawLine(x, u1, x, u2);
		if (e < 2) P->drawLine(x-2, u1, x+2, u1);
		else e-=2;
		if (e==0) P->drawLine(x-2, u2, x+2, u2);
	}
}

Point* ResidualPlot::getResidualFitData() const
{
    Point *points = new Point[nData[C][refTIsoT[Iso]][v]];
    for (int J=0; J < nData[C][refTIsoT[Iso]][v]; J++)
    {
        points[J].x = JData[C][refTIsoT[Iso]][v][J];
        points[J].y = Data[C][refTIsoT[Iso]][v][J];
        points[J].unc = Unc[C][refTIsoT[Iso]][v][J];
        points[J].obs = Obs[C][refTIsoT[Iso]][v][J];
        if (points[J].unc > 0) points[J].sig = 1.0 / points[J].unc;
        else points[J].sig = 1e-10;
    }
    return points;
}

void ResidualPlot::SetResidualFitData()
{
    IsoTab* IsoT = mol->getIso();
    Point* points = getResidualFitData();
    Block = true;
    m_ResidualFit->setFitData(nData[C][Iso][v], points, IsoT->relRedMass[Iso], state->getOmega(), mol->getJStep(Iso), state->getBe());
    Block = false;
    delete[] points;
    delete IsoT;
}

void ResidualPlot::ShowResidualFit()
{
    showResidualFit = !showResidualFit;
    if (showResidualFit && !m_ResidualFit->isFitDataAvailable() && m_ResidualFit->getNumberOfLocalPerturbations() > 0) SetResidualFitData();
    Paint();
}

void ResidualPlot::initData()
{
	Data = 0;
	NC = 0;
	nData = 0;
	JData = 0;
	NIso = 0;
	NJ = 0;
	Nv = 0;
	nv = 0;
	Marked = 0;
	Row = 0;
	Comp = 0;
}

void ResidualPlot::isoBoxChanged(int i)
{
	if (Block) return;
	//printf("ResidualPlot::isoBoxChanged\n");
	if (nv == 0 || i==-1)
	{
		Iso = i;
		vBox->clear();
		return;
	}
	int n, m, c, I;
	bool V;
	for (n=I=0; n<=i && I < NIso; I++) if (nv[I] > 0) n++;
	Iso = I-1;
	//printf("Iso=%d, NIso=%d, n=%d, i=%d, I=%d, nv[6]=%d\n", Iso, NIso, n, i, I, nv[6]);
	vBox->clear();
    if (refTIsoT == 0 || refTIsoT[Iso] == -1)
	{
		Iso = -1;
		return;
	}
	Block = true;
	for (n=m=0; n < Nv; n++)
	{
		for (V = false, c=0; c < NC; c++) if (nData[c][refTIsoT[Iso]][n] > 0) V = true;
		if (V) 
		{
			vBox->addItem(QString::number(n));
			if (n<v) m++;
		}
	}
	if (m >= vBox->count()) m = vBox->count() - 1;
	if (m>=0) vBox->setCurrentIndex(m);
	Block = false;
	//printf("Nv=%d, m=%d\n", vBox->count(), m);
	vBoxChanged(m);
	//printf("Ende isoBoxChanged\n");
}

void ResidualPlot::mark(int x1, int x2, int y1, int y2, bool CP)
{
	int n, x, y;
	for (n=0; n < nData[C][Iso][v]; n++)
	{
		x = int(XO + XSF * double(JData[C][Iso][v][n]));
		y = int(YO - YSF * Data[C][Iso][v][n]);
		if (x >= x1 && x <= x2 && y >= y1 && y <= y2) Marked[n] = true;
		else if (!CP) Marked[n] = false;
	}
	Paint();
}

void ResidualPlot::marked(QRect* R, bool CP)
{
	if (!ZoomB->isChecked()) mark(R->left(), R->right(), R->top(), R->bottom(), CP);
}

void ResidualPlot::markComp()
{
	int n, c;
	QDialog *D = new QDialog(this);
	QGridLayout *L;
	QComboBox *B;
	QPushButton *O, *Ca;
	D->setWindowTitle("Mark lines of component");
	L = new QGridLayout(D);
	L->addWidget(new QLabel("Component:", D), 0, 0);
	L->addWidget(B = new QComboBox(D), 0, 1);
	L->setRowMinimumHeight(1, 20);
	L->addWidget(O = new QPushButton("OK", D), 2, 0);
	L->addWidget(Ca = new QPushButton("Cancel", D), 2, 1);
	for (n=0; n < CompList.count(); n++) B->addItem(QString::number(CompList[n]));
	B->setEditable(false);
	connect(O, SIGNAL(clicked()), D, SLOT(accept()));
	connect(Ca, SIGNAL(clicked()), D, SLOT(reject()));
	if (D->exec() == QDialog::Accepted) 
	{
		for (n=0, c = B->currentText().toInt(); n < nData[C][refTIsoT[Iso]][v]; n++)
			Marked[n] = (Comp[C][refTIsoT[Iso]][v][n] == c ? true : false);
		Paint();
	}
	delete D;
}

void ResidualPlot::markLines()
{
	if (FitD == 0 && lineT == 0) return;
	int n, N;
	for (n=N=0; n < nData[C][Iso][v]; n++) if (Marked[n]) N++;
	int *rows = new int[N];
	for (n=N=0; n < nData[C][Iso][v]; n++) if (Marked[n]) rows[N++] = Row[C][Iso][v][n];
	if (lineT != 0) lineT->MarkLines(rows, N);
	else FitD->MarkLines(rows, N);
	delete[] rows;
}

void ResidualPlot::markUnc()
{
	QDialog *D = new QDialog(this);
	QGridLayout *L = new QGridLayout(D);
	QLineEdit *UE = new QLineEdit("0.02", D);
	QPushButton *OK = new QPushButton("OK", D), *Cancel = new QPushButton("Cancel", D);
	QComboBox *SelE = new QComboBox(D);
	double U;
	int n, s;
	D->setWindowTitle("Mark lines with uncertainty");
	L->addWidget(new QLabel("Uncertainty [cm^-1]:"), 0, 0);
	SelE->addItems(QStringList() << "<" << "<=" << "=" << ">=" << ">");
	SelE->setEditable(false);
	L->addWidget(SelE, 0, 1);
	L->addWidget(UE, 0, 2);
	L->setRowMinimumHeight(1, 20);
	L->addWidget(OK, 2, 0);
	L->addWidget(Cancel, 2, 2);
	connect(OK, SIGNAL(clicked()), D, SLOT(accept()));
	connect(Cancel, SIGNAL(clicked()), D, SLOT(reject()));
	if (D->exec() == QDialog::Accepted)
	{
		U = UE->text().toDouble();
		s = SelE->currentIndex();
		switch (s)
		{
			case 0:
				for (n=0; n < nData[C][Iso][v]; n++)
				{
					if (Unc[C][Iso][v][n] < U) Marked[n] = true;
					else Marked[n] = false;
				}
				break;
			case 1:
				for (n=0; n < nData[C][Iso][v]; n++)
				{
					if (Unc[C][Iso][v][n] <= U) Marked[n] = true;
					else Marked[n] = false;
				}
				break;
			case 2:
				for (n=0; n < nData[C][Iso][v]; n++)
				{
					if (Unc[C][Iso][v][n] == U) Marked[n] = true;
					else Marked[n] = false;
				}
				break;
			case 3:
				for (n=0; n < nData[C][Iso][v]; n++)
				{
					if (Unc[C][Iso][v][n] >= U) Marked[n] = true;
					else Marked[n] = false;
				}
				break;
			case 4:
				for (n=0; n < nData[C][Iso][v]; n++)
				{
					if (Unc[C][Iso][v][n] > U) Marked[n] = true;
					else Marked[n] = false;
				}
				break;
		}
		Paint();
	}
	delete D;
}

void ResidualPlot::molBoxChanged()
{
	if (Block) return;
	Molecule *nmol = MW->getMolecule(MolBox->currentText());
	if (mol != nmol)
	{
		if (mol != 0) disconnect(mol, SIGNAL(propertiesChanged()), this, SLOT(molBoxChanged()));
		if ((mol = nmol) != 0) connect(mol, SIGNAL(propertiesChanged()), this, SLOT(molBoxChanged()));
	}
	if (mol == 0)
	{
		StateBox->clear();
		return;
	}
	Block = true;
	StateBox->clear();
	int n, m=0, N = mol->getNumStates();
	ElState *SP;
	for (n=0; n<N; n++)
	{
		StateBox->addItem((SP = mol->getStateP(n))->getName());
		if (SP == state) StateBox->setCurrentIndex(m=n);
	}
	Block = false;
	stateBoxChanged(m);
}

void ResidualPlot::moleculesChanged()
{
	Block = true;
	Molecule *nmol;
	int n, m=-1, N = MW->getNumMolecules();
	MolBox->clear();
	for (n=0; n<N; n++) 
	{
		MolBox->addItem((nmol = MW->getMolecule(n))->getName());
		if (nmol == mol) m=n;
	}
	Block = false;
	if (m==-1) molBoxChanged();
}

void ResidualPlot::PictureClicked(QPoint* P, bool CP)
{
	mark(P->x() - 5, P->x() + 5, P->y() - 5, P->y() + 5, CP);
	DiagWindow::PictureClicked(P);
}

void ResidualPlot::PSpektrum(QPainter &P, const QRect &A, bool /*PrintFN*/)
{
	if (Block || C == -1 || C >= NC || Iso >= NIso || Iso == -1 
		|| v == -1 || v >= Nv) return;
	//printf("Residualplot::PSpektrum\n");
	int j=0, x, y, s = int(xStart->text().toDouble()), e = int(xStop->text().toDouble());
	int xm = ScaleYWidth + 2, xM = A.width() - 2, unc, unc1, unc2, end, NView = 0;
	int ym = ScaleTopOffset + 2, yM = A.height() - ScaleXHeight - 2, *viewR = new int[nData[C][Iso][v]];
	while (j < nData[C][Iso][v] ? JData[C][Iso][v][j] < s : false) j++;
	//printf("j=%d, nData=%d, JData=%d\n", j, nData[C][Iso][v], JData[C][Iso][v][j]);
	for ( ; j < nData[C][Iso][v] ? JData[C][Iso][v][j] <= e : false; j++)
	{
		x = int(XO + XSF * double(JData[C][Iso][v][j]));
		y = int(YO - YSF * Data[C][Iso][v][j]);
		unc = int(YSF * Unc[C][Iso][v][j]);
		end = 0;
		//printf("j=%d, x=%f, y=%f\n", j, x, y);
		if ((unc1 = y - unc) < ScaleTopOffset)
		{
			unc1 = ScaleTopOffset;
			end = 2;
		}
		if ((unc2 = y + unc) > A.height() - ScaleXHeight)
		{
			unc2 = A.height() - ScaleXHeight;
			end++;
		}
		if (Marked[j]) P.setPen(QColor(255, 0, 0));
		else P.setPen(QColor(0, 0, 0));
		if (x > xm && x < xM && y > ym && y < yM) 
		{
			drawDataPoint(&P, x, y, unc1, unc2, end);
			viewR[NView++] = (termT == 0 ? Row[C][Iso][v][j] : JData[C][Iso][v][j]);
		}
	}
    if (showResidualFit && NULL != m_ResidualFit)
    {
        int m, ly = 0, NPert = m_ResidualFit->getNumberOfLocalPerturbations();
        double Jval;
        for (x = xm; x < xM; x++)
        {
            y = (YO - YSF * m_ResidualFit->GetPointFromSplines((double(x) - XO) / XSF));
            if (y < ym) y = ym;
            else if (y > yM) y = yM;
            P.drawLine(x-1, ly, x, y);
            ly = y;
        }
        if (0 < NPert)
        {
            if (!m_ResidualFit->isFitDataAvailable()) SetResidualFitData();
            P.setPen(QColor(0, 0, 255));
            for (m=0; m < NPert; ++m) for (x = xm; x < xM; x++)
            {
                Jval = (double(x) - XO) / XSF;
                y = (YO - YSF * (m_ResidualFit->GetPointPerturber(Jval, m)));
                if (y < ym) y = ym;
                else if (y > yM) y = yM;
                if ((y > ym && y < yM) || (ly > ym && ly < yM)) P.drawLine(x-1, ly, x, y);
                ly = y;
            }
            P.setPen(QColor(0, 255, 0));
            for (m=0; m <= NPert; ++m) for (x = xm; x < xM; x++)
            {
                Jval = (double(x) - XO) / XSF;
                y = (YO - YSF * (m_ResidualFit->GetPoint(Jval, m)));
                if (y < ym) y = ym;
                else if (y > yM) y = yM;
                if ((y > ym && y < yM) || (ly > ym && ly < yM)) P.drawLine(x-1, ly, x, y);
                ly = y;
            }
        }
        DrawPoints(P, A);
    }
	if (termT != 0) termT->setViewnLevels(this, C, Iso, v, viewR, NView);
	else if (lineT != 0) lineT->setViewnRows(this, NView, viewR);
	else FitD->setViewnRows(this, NView, viewR);
	//printf("Ende PSpektrum\n");
}

void ResidualPlot::Printall()
{
	QPrinter Pr;
	QPrintDialog *PD = new QPrintDialog(&Pr);
	IsoTab *IsoT = mol->getIso();
	QString Caption;
	int v, Nv, I, NI = IsoBox->count(), c, Nc;
	if (PD->exec() == QDialog::Accepted)
	{
		QPainter P(&Pr);
        QRect R = Pr.pageLayout().paintRectPixels(Pr.resolution());
		for (I = 0; I < NI; I++)
		{
			IsoBox->setCurrentIndex(I);
			for (v=0, Nv = vBox->count(); v < Nv; v++)
			{
				vBox->setCurrentIndex(v);
				for (c=0, Nc = CompBox->count(); c < Nc; c++)
				{
					CompBox->setCurrentIndex(c);
					if (I>0 || v>0 || c>0) Pr.newPage();
					PaintScale(P, R);
					Caption = "Isotopologue " + IsoT->texName[Iso] + ", v'=" + vBox->currentText();
					if (NC > 1) Caption += ", " + CompBox->currentText();
					WriteText(P, R.center().x() - TextWidth(ScaleFont, Caption), 
							  TextHeight(ScaleFont, Caption), Caption, ScaleFont, 0);
					PSpektrum(P, R, false);
				}
			}
		}
	}
	delete PD;
	delete IsoT;
}

void ResidualPlot::refSourceBoxChanged()
{
	if (Block) return;
	int i = RefBox->currentIndex(); 
    TermTable *nRef = (state != 0 && RefBox->currentText() != "Potential" && i < state->getNumTermTables() && i>=0 ?
                state->getTermTable(i >= addedToFront ? i - addedToFront : i) : 0);
	if (refT != nRef)
	{
		if (refT != 0) 
			disconnect(refT, SIGNAL(propertiesChanged()), this, SLOT(refSourceBoxChanged()));
		if ((refT = nRef) != 0) 
			connect(refT, SIGNAL(propertiesChanged()), this, SLOT(refSourceBoxChanged()));
	}
	sourceBoxChanged();
}

void ResidualPlot::RemoveLocalPerturbations()
{
    int N = m_ResidualFit->getNumberOfLocalPerturbations(), n, m, M=0;
    double Center[N];
    bool Selected[N];
    for (n=0; n<N; ++n)
    {
        Center[n] = m_ResidualFit->getLocalPerturbation(n)->GetCenter();
        Selected[N] = false;
    }
    for (m=0; m < nData[C][Iso][v]; ++m) if (Marked[m])
    {
        ++M;
        if (m>0 && Marked[m-1]) for (n=0; n<N; ++n) if (JData[C][Iso][v][m-1] < Center[n] && JData[C][Iso][v][m] > Center[n]) Selected[n] = true;
    }
    for (n=N-1; n>=0; n--) if (M==0 || Selected[n]) m_ResidualFit->RemoveLocalPerturbation(n);
    Paint();
}

void ResidualPlot::sub9fromError()
{
	if (lineT == 0 && FitD == 0) return;
	int n, N;
	for (n=N=0; n < nData[C][Iso][v]; n++) if (Marked[n]) N++;
	int *rows = new int[N];
	double *u = new double[N];
	for (n=N=0; n < nData[C][Iso][v]; n++) if (Marked[n])
	{
		u[N] = (Unc[C][Iso][v][n] > 9.0 ? Unc[C][Iso][v][n] - 9.0 : Unc[C][Iso][v][n]);
		rows[N++] = Row[C][Iso][v][n];
	}
	if (lineT != 0) lineT->setUncertainty(rows, u, N);
	else FitD->setUncertainty(rows, u, N);
	delete[] rows;
	delete[] u;
}

void ResidualPlot::setError()
{
	if (FitD == 0 && lineT == 0) return;
	QDialog *D = new QDialog(this);
	QGridLayout *L = new QGridLayout(D);
	QLineEdit *E = new QLineEdit("0.01", D);
	QPushButton *OK = new QPushButton("OK", this), *Cancel = new QPushButton("Cancel", this);
	D->setWindowTitle("Uncertainty to set");
	L->addWidget(new QLabel("Uncertainty [cm^-1]"), 0, 0);
	L->addWidget(E, 0, 1);
	L->setRowMinimumHeight(1, 20);
	L->addWidget(OK, 2, 0);
	L->addWidget(Cancel, 2, 1);
	connect(OK, SIGNAL(clicked()), D, SLOT(accept()));
	connect(Cancel, SIGNAL(clicked()), D, SLOT(reject()));
	if (D->exec() == QDialog::Accepted)
	{
		int n, N;
		for (n=N=0; n < nData[C][Iso][v]; n++) if (Marked[n]) N++;
		int *rows = new int[N];
		double *unc = new double[N], nu = E->text().toDouble();
		for (n=N=0; n < nData[C][Iso][v]; n++) if (Marked[n])
		{
			rows[N] = Row[C][Iso][v][n];
			unc[N++] = nu;
		}
		if (lineT != 0) lineT->setUncertainty(rows, unc, N);
		else FitD->setUncertainty(rows, unc, N);
		delete[] rows;
		delete[] unc;
	}
	delete D;
}

void ResidualPlot::showPMenu(QPoint* P)
{
    if (!ZoomB->isChecked())
    {
        residualFitSeparatorAct->setVisible(false);
        showResidualFitAct->setVisible(false);
        showResidualFitSplinePointsAct->setVisible(false);
        createResidualFitAct->setVisible(false);
        addResidualFitSplinePointAct->setVisible(false);
        deleteResidualFitSplinePointAct->setVisible(false);
        joinResidualFitSplinesAct->setVisible(false);
        addLocalPerturbationAct->setVisible(false);
        removeLocalPerturbationsAct->setVisible(false);
        showLocalPerturbationInfoTableAct->setVisible(false);
        if (FitD != 0)
        {
            residualFitSeparatorAct->setVisible(true);
            if (m_ResidualFit != 0)
            {
                showResidualFitAct->setText(showResidualFit ? "Hide &residual fit" : "Show &residual fit");
                showResidualFitAct->setVisible(true);
                createResidualFitAct->setText("Re&create residual fit");
                if (showResidualFit)
                {
                    showResidualFitSplinePointsAct->setText(showPoints ? "Hide spline points" : "Show spline points");
                    showResidualFitSplinePointsAct->setVisible(true);
                    joinResidualFitSplinesAct->setVisible(m_ResidualFit->getNumSplines() > 1);
                    addLocalPerturbationAct->setVisible(true);
                    if (m_ResidualFit->getNumberOfLocalPerturbations() > 0)
                    {
                        removeLocalPerturbationsAct->setVisible(true);
                        showLocalPerturbationInfoTableAct->setVisible(true);
                    }
                    if (showPoints)
                    {
                        if (sPoint == -1) addResidualFitSplinePointAct->setVisible(true);
                        else deleteResidualFitSplinePointAct->setVisible(true);
                    }
                }
            }
            createResidualFitAct->setVisible(true);
        }
        PMenu->popup(*P);
    }
}

void ResidualPlot::sourceBoxChanged()
{
	if (Block) return;
	//printf("ResidualPlot::sourceBoxChanged\n");
    int n, t, m, j, k, r, N, I = SourceBox->currentIndex() - addedToFront;
	destroyData();
	TermTable *nTerm = 0;
	LineTable *nLine = 0;
    FitData *nFitD = (FitD != 0 && SourceBox->currentText() == FitD->getName() ? FitD : 0);
	Transition *T=0;
	Marked = 0;
    if (state != 0 && I >= 0 && nFitD == 0)
	{
		n = state->getNumTermTables();
		if (I<n) nTerm = state->getTermTable(I);
		else if (I < (m = n + state->getNumFitDataSets())) nFitD = state->getFitData(I-n); 
		else
		{
			for (t=0; m<=I; t++) 
			{
				n=m;
				T = mol->getTransitionP(t);
				if (T->getLowerState() == state || T->getUpperState() == state)
					m = n + T->getNumLineTables();
			}
			nLine = T->getLineTable(I-n);
		}
		//printf("I=%d, n=%d, nTerm=%d, nLine=%d, t=%d\n", I, n, nTerm, nLine, t);
	}
	//printf("Nach state\n");
	if (nTerm != termT || nLine != lineT || nFitD != FitD)
	{
		//printf("Vor disconnect, termT=%d\n", termT);
		if (termT != 0) 
		{
			disconnect(termT, SIGNAL(propertiesChanged()), this, SLOT(sourceBoxChanged()));
			termT->setViewnRows(this, 0, 0);
		}
		if (FitD != 0) 
		{
			disconnect(FitD, SIGNAL(propertiesChanged()), this, SLOT(sourceBoxChanged()));
			disconnect(FitD, SIGNAL(AssignmentsAccepted(FitData*)), this, SLOT(sourceBoxChanged()));
			FitD->setViewnRows(this, 0, 0);
		}
		//printf("Vor lineT\n");
		if (lineT != 0) 
		{
			disconnect(lineT, SIGNAL(propertiesChanged()), this, SLOT(sourceBoxChanged()));
			disconnect(lineT, SIGNAL(DataChanged()), this, SLOT(sourceBoxChanged()));
			lineT->setViewnRows(this, 0, 0);
		}
		//printf("vor connect\n");
		if ((termT = nTerm) != 0) 
			connect(termT, SIGNAL(propertiesChanged()), this, SLOT(sourceBoxChanged()));
		if ((FitD = nFitD) != 0) 
		{
			connect(FitD, SIGNAL(propertiesChanged()), this, SLOT(sourceBoxChanged()));
			connect(FitD, SIGNAL(AssignmentsAccepted(FitData*)), this, SLOT(sourceBoxChanged()));
		}
		if ((lineT = nLine) != 0)
		{
			connect(lineT, SIGNAL(propertiesChanged()), this, SLOT(sourceBoxChanged()));
			connect(lineT, SIGNAL(DataChanged()), this, SLOT(sourceBoxChanged()));
		}
		//printf("Nach connect\n");
	}
	initData();
	if ((refT == 0 && FitD == 0) || (lineT == 0 && termT == 0 && FitD == 0))
	{
		//printf("Vor initdata\n");
		IsoBox->clear();
		//printf("ResidualPlot::SourceBoxChanged: No data!\n");
		return;
	}
	//printf("Vor getData\n");
	double ****RData;
	int rNC, rNIso, rNv, rNJ;
	if (refT != 0)
	{
		RData = refT->getData();
		rNC = refT->getNumComp();
		rNIso = refT->getNumIso();
		rNv = refT->getMaxv() + 1;
		rNJ = refT->getMaxJ() + 1;
		if (refTIsoT != 0) delete[] refTIsoT;
		if (refTcT != 0) delete[] refTcT;
		refTIsoT = refT->getIsoT();
		refT->getCompT(MrefTcT, refTcT);
	}
	else
	{
		RData = 0;
		rNC = rNIso = rNv = rNJ = 0;
	}
	if (termT != 0)
	{
		//printf("Term energies\n");
		double ****tData = termT->getData();
		NC = termT->getNumComp();
		if ((NIso = termT->getNumIso()) > rNIso) NIso = rNIso;
		if ((Nv = termT->getMaxv() + 1) > rNv) Nv = rNv;
		if ((NJ = termT->getMaxJ() + 1) > rNJ) NJ = rNJ;
		nData = CreateInt(NC, NIso, Nv);
		Data = new double***[NC];
		JData = new int***[NC];
		Unc = new double***[NC];
        Obs = new double***[NC];
		for (n=0; n < NC; n++) 
		{
			r = (n < rNC ? n : 0);
			Data[n] = new double**[NIso];
			JData[n] = new int**[NIso];
			Unc[n] = new double**[NIso];
            Obs[n] = new double**[NIso];
			for (m=0; m < NIso; m++) 
			{
				Data[n][m] = new double*[Nv];
				JData[n][m] = new int*[Nv];
				Unc[n][m] = new double*[Nv];
                Obs[n][m] = new double*[Nv];
				for (t=0; t < Nv; t++) 
				{
					for (j=0, nData[n][m][t] = 0; j < NJ; j++)
						if (tData[n][m][t][j] != 0.0 && RData[r][m][t][j] != 0.0) nData[n][m][t]++;
					if (nData[n][m][t] > 0)
					{
						Data[n][m][t] = new double[nData[n][m][t]];
						JData[n][m][t] = new int[nData[n][m][t]];
						Unc[n][m][t] = new double[nData[n][m][t]];
                        Obs[n][m][t] = new double[nData[n][m][t]];
						for (j=k=0; j < NJ; j++) 
							if (tData[n][m][t][j] != 0.0 && RData[r][m][t][j] != 0.0)
						{
							Data[n][m][t][k] = tData[n][m][t][j] - RData[r][m][t][j];
							Unc[n][m][t][k] = 0.0;
                            JData[n][m][t][k] = j;
                            Obs[n][m][t][k++] = tData[n][m][t][j];
						}
					}
					else
					{
						Data[n][m][t] = 0;
						JData[n][m][t] = 0;
						Unc[n][m][t] = 0;
                        Obs[n][m][t] = 0;
					}
				}
			}
		}
	}
	else if (FitD != 0)
	{
		TableLine *Lines;
		int PN, Js, vs, IsoN;
		double ES;
		bool FSA = false, FSP = (state->getS() > 0);
        FitD->getData(Lines, N, -1, -2, -2, 0, -1, state);
        for (n = NIso = Nv = NJ = 0, NC = 1; n<N; n++) if (Lines[n].State == 0 || Lines[n].State == state)
		{
			if (Lines[n].Js == Lines[n].Jss && Lines[n].isTE) NC = 2;
			if (Lines[n].Iso >= NIso) NIso = Lines[n].Iso + 1;
			if (Lines[n].vss >= Nv) Nv = Lines[n].vss + 1;
			if (Lines[n].Jss >= NJ) NJ = Lines[n].Jss + 1;
			if (Lines[n].FC >= NC) NC = Lines[n].FC + 1;
			if (FSP && Lines[n].FC > 0) FSA = true;
		}
		if (NC > 10) NC = 10;
		if (RData != 0)
		{
			if (NIso > rNIso) NIso = rNIso;
			if (Nv > rNv) Nv = rNv;
			if (NJ > rNJ) NJ = rNJ;
		}
		nData = CreateInt(NC, NIso, Nv);
		Data = new double***[NC];
		JData = new int***[NC];
		Unc = new double***[NC];
        Obs = new double***[NC];
		Row = new int***[NC];
		Comp = new int***[NC];
		for (n=0; n < NC; n++)
		{
			Data[n] = new double**[NIso];
			JData[n] = new int**[NIso];
			Unc[n] = new double**[NIso];
            Obs[n] = new double**[NIso];
			Row[n] = new int**[NIso];
			Comp[n] = new int**[NIso];
			for (m=0; m < NIso; m++)
			{
				Data[n][m] = new double*[Nv];
				JData[n][m] = new int*[Nv];
				Unc[n][m] = new double*[Nv];
                Obs[n][m] = new double*[Nv];
				Row[n][m] = new int*[Nv];
				Comp[n][m] = new int*[Nv];
				for (t=0; t < Nv; t++) nData[n][m][t] = 0;
			}
		}
        if (refTIsoT != 0) for (n=0; n<N; n++) if (Lines[n].State == 0 || Lines[n].State == state)
		{
			t = (FSA ? (Lines[n].FC >= 0 && Lines[n].FC < 10 ? Lines[n].FC : 0) 
					 : (Lines[n].Jss == Lines[n].Js && NC > 1 && Lines[n].isTE ? 1 : 0));
			r = (t < rNC ? t : 0);
            if (refTIsoT[Lines[n].Iso] != -1 && Lines[n].vss < Nv && Lines[n].vss >= 0 && Lines[n].Jss < NJ
				   && (RData == 0 || RData[r][refTIsoT[Lines[n].Iso]][Lines[n].vss][Lines[n].Jss] != 0.0))
				nData[t][refTIsoT[Lines[n].Iso]][Lines[n].vss]++;
		}
		for (n=0; n < NC; n++) for (m=0; m < NIso; m++) for (t=0; t < Nv; t++)
		{
			if (nData[n][m][t] > 0)
			{
				Data[n][m][t] = new double[nData[n][m][t]];
				JData[n][m][t] = new int[nData[n][m][t]];
				Unc[n][m][t] = new double[nData[n][m][t]];
                Obs[n][m][t] = new double[nData[n][m][t]];
				Row[n][m][t] = new int[nData[n][m][t]];
				Comp[n][m][t] = new int[nData[n][m][t]];
				nData[n][m][t] = 0;
			}
			else
			{
				Data[n][m][t] = 0;
				JData[n][m][t] = 0;
				Unc[n][m][t] = 0;
                Obs[n][m][t] = 0;
				Row[n][m][t] = 0;
				Comp[n][m][t] = 0;
			}
		}
        for (n=0; n<N; n++) if (Lines[n].State == 0 || Lines[n].State == state)
		{
			t = (FSA ? (Lines[n].FC >= 0 && Lines[n].FC < 10 ? Lines[n].FC : 0) 
			         : (Lines[n].Jss == Lines[n].Js && NC > 1 && Lines[n].isTE ? 1 : 0));
			r = (t < rNC ? t : 0);
			if (!Lines[n].isTE && RData != 0)
			{
				for (m=n, Js = Lines[n].Js, vs = Lines[n].vs, IsoN = Lines[n].Iso, 
						PN = Lines[n].PN; 
					 (n-1<N ? Lines[n+1].PN == PN && Lines[n+1].Iso == IsoN 
						&& Lines[n+1].vs == vs && Lines[n+1].Js == Js : false); n++) ;
				for (j=m, k=0, ES = 0.0; j<=n; j++)
					if (refTIsoT[Lines[n].Iso] != -1 && Lines[n].vss < Nv && Lines[n].vss >= 0 && Lines[n].Jss < NJ ?
						RData[r][refTIsoT[Lines[n].Iso]][Lines[n].vss][Lines[n].Jss] != 0.0 : false)
				{
					ES += RData[r][refTIsoT[Lines[n].Iso]][Lines[n].vss][Lines[n].Jss] + Lines[n].WN;
					k++;
				}
				ES /= k;
				for (j=m, k=0; j<=n; j++)
					if (refTIsoT[Lines[j].Iso] != -1 && Lines[j].vss < Nv && Lines[n].vss >= 0 && Lines[j].Jss < NJ ?
						RData[r][refTIsoT[Lines[j].Iso]][Lines[j].vss][Lines[j].Jss] != 0.0 : false)
				{
					Data[t][refTIsoT[Lines[j].Iso]][Lines[j].vss][nData[t][refTIsoT[Lines[j].Iso]][Lines[j].vss]]
						= ES - Lines[j].WN - RData[r][refTIsoT[Lines[j].Iso]][Lines[j].vss][Lines[j].Jss];
					JData[t][refTIsoT[Lines[j].Iso]][Lines[j].vss][nData[t][refTIsoT[Lines[j].Iso]][Lines[j].vss]]
						= Lines[j].Jss;
					Unc[t][refTIsoT[Lines[j].Iso]][Lines[j].vss][nData[t][refTIsoT[Lines[j].Iso]][Lines[j].vss]] = Lines[j].err;
                    Obs[t][refTIsoT[Lines[j].Iso]][Lines[j].vss][nData[t][refTIsoT[Lines[j].Iso]][Lines[j].vss]] = ES - Lines[j].WN;
					Comp[t][refTIsoT[Lines[j].Iso]][Lines[j].vss][nData[t][refTIsoT[Lines[j].Iso]][Lines[j].vss]] = 0;
					Row[t][refTIsoT[Lines[j].Iso]][Lines[j].vss][nData[t][refTIsoT[Lines[j].Iso]][Lines[j].vss]++] = Lines[j].Row;
				}
			}
            else if (refTIsoT != 0 && (RData != 0 ? (refTIsoT[Lines[n].Iso] != -1 && Lines[n].vss < Nv && Lines[n].vss >= 0 && Lines[n].Jss < NJ ?
                RData[r][refTIsoT[Lines[n].Iso]][Lines[n].vss][Lines[n].Jss] != 0.0 : false) : Lines[n].vss >= 0))
			{
				Data[t][refTIsoT[Lines[n].Iso]][Lines[n].vss][nData[t][refTIsoT[Lines[n].Iso]][Lines[n].vss]]
					= (RData != 0 ? Lines[n].WN - RData[r][refTIsoT[Lines[n].Iso]][Lines[n].vss][Lines[n].Jss] : Lines[n].dev);
				Unc[t][refTIsoT[Lines[n].Iso]][Lines[n].vss][nData[t][refTIsoT[Lines[n].Iso]][Lines[n].vss]] = Lines[n].err;
                Obs[t][refTIsoT[Lines[n].Iso]][Lines[n].vss][nData[t][refTIsoT[Lines[n].Iso]][Lines[n].vss]] = Lines[n].WN;
				JData[t][refTIsoT[Lines[n].Iso]][Lines[n].vss][nData[t][refTIsoT[Lines[n].Iso]][Lines[n].vss]] 
					= Lines[n].Jss;
				Comp[t][refTIsoT[Lines[n].Iso]][Lines[n].vss][nData[t][refTIsoT[Lines[n].Iso]][Lines[n].vss]] = -1 - Lines[n].vs;
				Row[t][refTIsoT[Lines[n].Iso]][Lines[n].vss][nData[t][refTIsoT[Lines[n].Iso]][Lines[n].vss]++] = Lines[n].Row;
			}
		}
		delete[] Lines;
	}
	else 
	{
		//printf("Lines\n");
		if (state == T->getUpperState())
		{
			TermTable *lT = (T->getLowerState() != 0 ? T->getLowerState()->getTermTable() : 0);
			if (lT == 0)
			{
				initData();
				IsoBox->clear();
				printf("ResidualPlot::SourceBoxChanged: No data for lower state!\n");
				return;
			}
			N = lineT->getAnzahlLinien();
			int **Z = CreateInt(N, 6), lNIso = lT->getNumIso();
			int lNv = lT->getMaxv() + 1, lNJ = lT->getMaxJ() + 1;
			double *E = new double[N], ****lData = lT->getData(), *Un = new double[N];
			int MlTcT, *lTcT, *lTIsoT = lT->getIsoT(), lt;
			lT->getCompT(MlTcT, lTcT);
			lineT->getLines(Z, E, Un);
			for (n = NIso = Nv = NJ = 0, NC = 1; n<N; n++)
			{
				if (Z[n][2] == Z[n][4] && NC == 1) NC = 2;
				if (Z[n][0] >= NIso) NIso = Z[n][0] + 1;
				if (Z[n][1] >= Nv) Nv = Z[n][1] + 1;
				if (Z[n][2] >= NJ) NJ = Z[n][2] + 1;
				if (Z[n][5] >= NC) NC = Z[n][5] + 1;
			}
			if (NIso > rNIso) NIso = rNIso;
			if (NIso > lNIso) NIso = lNIso;
			if (Nv > rNv) Nv = rNv;
			if (NJ > rNJ) NJ = rNJ;
			if (NJ > lNJ + 1) NJ = lNJ + 1;
			//printf("NJ=%d, NIso=%d, NC=%d, Nv=%d\n", NJ, NIso, NC, Nv);
			nData = CreateInt(NC, NIso, Nv);
			Data = new double***[NC];
			JData = new int***[NC];
			Unc = new double***[NC];
            Obs = new double***[NC];
			Row = new int***[NC];
			for (n=0; n < NC; n++)
			{
				Data[n] = new double**[NIso];
				JData[n] = new int**[NIso];
				Unc[n] = new double**[NIso];
                Obs[n] = new double**[NIso];
				Row[n] = new int**[NIso];
				for (m=0; m < NIso; m++)
				{
					Data[n][m] = new double*[Nv];
					JData[n][m] = new int*[Nv];
					Unc[n][m] = new double*[Nv];
                    Obs[n][m] = new double*[Nv];
					Row[n][m] = new int*[Nv];
					for (t=0; t < Nv; t++) nData[n][m][t] = 0;
				}
			}
			FSA = false;
			for (n=0; n<N; n++) 
			{
				if (Z[n][5] >= 0)
				{
					t = Z[n][5];
					r = ((t <= MrefTcT ? refTcT[t] >= 0 : false) ? refTcT[t] : 0); 
					lt = ((t <= MlTcT ? lTcT[t] >= 0 : false) ? lTcT[t] : 0);
					FSA = true;
				}
				else
				{
					t = (Z[n][2] == Z[n][4] && NC > 1 ? 1 : 0); 
					r = (t < rNC ? t : 0);
					lt = 0;
				}
				if ((lTIsoT[Z[n][0]] >= 0 && Z[n][3] < lNv && Z[n][4] < lNJ ?
							lData[lt][lTIsoT[Z[n][0]]][Z[n][3]][Z[n][4]] != 0.0 : false) 
						&& (refTIsoT[Z[n][0]] >= 0 && Z[n][1] < Nv && Z[n][1] >= 0 && Z[n][2] < NJ ? 
							RData[r][refTIsoT[Z[n][0]]][Z[n][1]][Z[n][2]] != 0.0 : false))
					nData[t][refTIsoT[Z[n][0]]][Z[n][1]]++;
			}
			for (n=0; n < NC; n++) for (m=0; m < NIso; m++) for (t=0; t < Nv; t++)
			{
				if (nData[n][m][t] > 0)
				{
					Data[n][m][t] = new double[nData[n][m][t]];
					JData[n][m][t] = new int[nData[n][m][t]];
					Unc[n][m][t] = new double[nData[n][m][t]];
                    Obs[n][m][t] = new double[nData[n][m][t]];
					Row[n][m][t] = new int[nData[n][m][t]];
					nData[n][m][t] = 0;
				}
				else
				{
					Data[n][m][t] = 0;
					JData[n][m][t] = 0;
					Unc[n][m][t] = 0;
                    Obs[n][m][t] = 0;
					Row[n][m][t] = 0;
				}
			}
			for (n=0; n<N; n++)
			{
				//printf("N=%d, ", n, N);
				//printf("Z[%d][0]=%d, Z[%d][1]=%d, Z[%d][3]=%d, Z[%d][4]=%d\n", 
					//	   n, Z[n][0], n, Z[n][1], n, Z[n][3], n, Z[n][4]);
				if (Z[n][5] >= 0)
				{
					t = Z[n][5];
					r = ((t <= MrefTcT ? refTcT[t] >= 0 : false) ? refTcT[t] : 0); 
					lt = ((t <= MlTcT ? lTcT[t] >= 0 : false) ? lTcT[t] : 0);
				}
				else
				{
					t = (Z[n][2] == Z[n][4] && NC > 1 ? 1 : 0); 
					r = (t < rNC ? t : 0);
					lt = 0;
				}
				if ((lTIsoT[Z[n][0]] >= 0 && Z[n][3] < lNv && Z[n][4] < lNJ ?
							lData[lt][lTIsoT[Z[n][0]]][Z[n][3]][Z[n][4]] != 0.0 : false) 
						&& (refTIsoT[Z[n][0]] >= 0 && Z[n][1] < Nv && Z[n][1] >= 0 && Z[n][2] < NJ ? 
							RData[r][refTIsoT[Z[n][0]]][Z[n][1]][Z[n][2]] != 0.0 : false))
				{
					Data[t][refTIsoT[Z[n][0]]][Z[n][1]][nData[t][refTIsoT[Z[n][0]]][Z[n][1]]] = E[n]
							- RData[r][refTIsoT[Z[n][0]]][Z[n][1]][Z[n][2]] 
							+ lData[lt][lTIsoT[Z[n][0]]][Z[n][3]][Z[n][4]];
					Unc[t][refTIsoT[Z[n][0]]][Z[n][1]][nData[t][refTIsoT[Z[n][0]]][Z[n][1]]] = Un[n];
                    Obs[t][refTIsoT[Z[n][0]]][Z[n][1]][nData[t][refTIsoT[Z[n][0]]][Z[n][1]]] = E[n] + lData[lt][lTIsoT[Z[n][0]]][Z[n][3]][Z[n][4]];
					JData[t][refTIsoT[Z[n][0]]][Z[n][1]][nData[t][refTIsoT[Z[n][0]]][Z[n][1]]] = Z[n][2];
					Row[t][refTIsoT[Z[n][0]]][Z[n][1]][nData[t][refTIsoT[Z[n][0]]][Z[n][1]]++] = n;
				}
			}
			Destroy(Z, N);
			delete[] E;
			delete[] Un;
			delete[] lTIsoT;
			delete[] lTcT;
		}
		else
		{
			Progression *P;
			lineT->getProgressions(N, P);
			for (n = NIso = Nv = NJ = 0, NC = 1; n<N; n++)
			{
				if (P[n].Js != P[n].L[0].Jss) NC = 2;
				if (P[n].Iso >= NIso) NIso = P[n].Iso + 1;
				for (m=0; m < P[n].N; m++) 
				{
					if (P[n].L[m].vss >= Nv) Nv = P[n].L[m].vss + 1;
					if (P[n].L[m].Jss >= NJ) NJ = P[n].L[m].Jss + 1;
				}
			}
			//printf("NIso=%d, rNIso=%d, Nv=%d, rNv=%d, NJ=%d, rNJ=%d, N=%d\n", 
				//   NIso, rNIso, Nv, rNv, NJ, rNJ, N);
			if (NIso > rNIso) NIso = rNIso;
			if (Nv > rNv) Nv = rNv;
			if (NJ > rNJ) NJ = rNJ;
			nData = CreateInt(NC, NIso, Nv);
			Data = new double***[NC];
			JData = new int***[NC];
			Unc = new double***[NC];
            Obs = new double***[NC];
			Row = new int***[NC];
			for (n=0; n < NC; n++)
			{
				Data[n] = new double**[NIso];
				JData[n] = new int**[NIso];
				Unc[n] = new double**[NIso];
                Obs[n] = new double**[NIso];
				Row[n] = new int**[NIso];
				for (m=0; m < NIso; m++)
				{
					Data[n][m] = new double*[Nv];
					JData[n][m] = new int*[Nv];
					Unc[n][m] = new double*[Nv];
                    Obs[n][m] = new double*[Nv];
					Row[n][m] = new int*[Nv];
					for (t=0; t < Nv; t++) nData[n][m][t] = 0;
				}
			}
			//printf("1.\n");
			for (n=0; n<N; n++) 
			{
				if (P[n].Js != P[n].L[0].Jss && NC > 1) t=1;
				else t=0;
				r = (t < rNC ? t : 0);
				if (P[n].Iso < NIso) for (m=0; m < P[n].N; m++) 
						if (P[n].L[m].vss < Nv && P[n].L[m].Jss < NJ ? 
								RData[r][P[n].Iso][P[n].L[m].vss][P[n].L[m].Jss] != 0.0 : false)
							nData[t][P[n].Iso][P[n].L[m].vss]++;
			}
			//printf("2.\n");
			for (n=0; n < NC; n++) for (m=0; m < NIso; m++) for (t=0; t < Nv; t++)
			{
				if (nData[n][m][t] > 0)
				{
					//printf("nData[%d][%d][%d]=%d\n", n, m, t, nData[n][m][t]);
					Data[n][m][t] = new double[nData[n][m][t]];
					JData[n][m][t] = new int[nData[n][m][t]];
					Unc[n][m][t] = new double[nData[n][m][t]];
                    Obs[n][m][t] = new double[nData[n][m][t]];
					Row[n][m][t] = new int[nData[n][m][t]];
					nData[n][m][t] = 0;
				}
				else
				{
					Data[n][m][t] = 0;
					JData[n][m][t] = 0;
					Unc[n][m][t] = 0;
                    Obs[n][m][t] = 0;
					Row[n][m][t] = 0;
				}
			}
			//printf("3.\n");
			double ES;
			for (n=0; n<N; n++)
			{
				//printf("n=%d, N=%d\n", n, N);
				if (P[n].Js != P[n].L[0].Jss && NC > 1) t=1;
				else t=0;
				r = (t < rNC ? t : 0);
				if (P[n].Iso < NIso)
				{
					for (m=j=0, ES = 0.0; m < P[n].N; m++)
						if (P[n].L[m].vss < Nv && P[n].L[m].Jss < NJ ?
								RData[r][P[n].Iso][P[n].L[m].vss][P[n].L[m].Jss] != 0.0 : false)
					{
						ES += P[n].L[m].E + RData[r][P[n].Iso][P[n].L[m].vss][P[n].L[m].Jss];
						j++;
					}
					ES /= j;
					//printf("ES=%f, j=%d\n", ES, j);
					//printf("t=%d, vss=%d, Iso=%d, Jss=%d, P[%d].N=%d\n", t, P[n].L[0].vss, P[n].Iso,
						//   P[n].L[0].Jss, n, P[n].N);
					for (m=0; m < P[n].N; m++)
						if (P[n].L[m].vss < Nv && P[n].L[m].Jss < NJ ?
								RData[r][P[n].Iso][P[n].L[m].vss][P[n].L[m].Jss] != 0.0 : false)
					{
						Data[t][P[n].Iso][P[n].L[m].vss][nData[t][P[n].Iso][P[n].L[m].vss]] =
							ES - P[n].L[m].E - RData[r][P[n].Iso][P[n].L[m].vss][P[n].L[m].Jss];
						Unc[t][P[n].Iso][P[n].L[m].vss][nData[t][P[n].Iso][P[n].L[m].vss]] = P[n].L[m].err;
                        Obs[t][P[n].Iso][P[n].L[m].vss][nData[t][P[n].Iso][P[n].L[m].vss]] = ES - P[n].L[m].E;
						JData[t][P[n].Iso][P[n].L[m].vss][nData[t][P[n].Iso][P[n].L[m].vss]] =
								P[n].L[m].Jss;
						Row[t][P[n].Iso][P[n].L[m].vss][nData[t][P[n].Iso][P[n].L[m].vss]++] = P[n].L[m].row;
					}
				}
			}
			delete[] P;
		}
	}
	if (termT == 0)
	{
		//printf("4.\n");
		double B;
		for (n=0; n < NC; n++) for (m=0; m < NIso; m++) for (t=0; t < Nv; t++) for (j=0; j!=-1; )
			for (j=-1, k=1; k < nData[n][m][t]; k++) if (JData[n][m][t][k-1] > JData[n][m][t][k])
		{
			j = JData[n][m][t][k-1];
			JData[n][m][t][k-1] = JData[n][m][t][k];
			JData[n][m][t][k] = j;
			B = Data[n][m][t][k-1];
			Data[n][m][t][k-1] = Data[n][m][t][k];
			Data[n][m][t][k] = B;
			B = Unc[n][m][t][k-1];
			Unc[n][m][t][k-1] = Unc[n][m][t][k];
			Unc[n][m][t][k] = B;
            B = Obs[n][m][t][k-1];
            Obs[n][m][t][k-1] = Obs[n][m][t][k];
            Obs[n][m][t][k] = B;
			j = Row[n][m][t][k-1];
			Row[n][m][t][k-1] = Row[n][m][t][k];
			Row[n][m][t][k] = j;
		}
	}
	//printf("Vor Iso\n");
	nv = new int[NIso];
	for (n=0; n < NIso; n++) nv[n] = 0;
	for (n=0; n < NC; n++) for (m=0; m < NIso; m++) for (t=0; t < Nv; t++) if (nData[n][m][t] > 0) 
					nv[m]++;
	Block = true;
	IsoBox->clear();
	IsoTab *IT = mol->getIso();
	int *IsoZ;
	if (refT != 0 ? (IsoZ = refT->getIsoZ()) != 0 : false)
	{
		for (n=m=0; n < NIso; n++) if (nv[n] > 0) 
		{
			//printf("addIso %d, m=%d\n", n, m);
			IsoBox->addItem(IT->getIsoName(IsoZ[n]));
			if (n < Iso) m++;
		}
	}
	else for (n=m=0; n < NIso; n++) if (nv[n] > 0)
	{
		IsoBox->addItem(IT->getIsoName(n));
		if (n < Iso) m++;
	}
	if (m >= IsoBox->count()) m=0;
	//printf("Iso=%d, m=%d\n", Iso, m);
	delete IT;
	IsoBox->setCurrentIndex(m);
	Block = false;
	//printf("Vor IsoBox\n");
	isoBoxChanged(m);
	//printf("Ende sourceBoxChanged, Nv=%d, NJ=%d, NC=%d, NIso=%d\n", Nv, NJ, NC, NIso);
}

void ResidualPlot::stateBoxChanged(int i)
{
	if (Block) return;
	Block = true;
    state = (mol != 0 ? mol->getStateP(i) : 0);
    bool SV = true, RV = true, containsState = (FitD != 0 && state != 0 && FitD->containsState(state));
	int n, N, m, c, M;
    QString Name, SourceText = SourceBox->currentText(), RefText = RefBox->currentText();
	Transition *T;
	SourceBox->clear();
	RefBox->clear();
	if (state != 0)
	{
        if (containsState)
        {
            SourceBox->addItem(SourceText);
            RefBox->addItem(RefText);
            addedToFront = 1;
        }
        else addedToFront = 0;
        for (n=0, c = N = state->getNumTermTables(); n<N; n++)
		{
			Name = state->getTermTableName(n);
			if (0 == addedToFront || Name != SourceText)
			{
				SourceBox->addItem(Name);
				if (refT != 0 ? Name == refT->getName() : false)
				{
					RefBox->setCurrentIndex(n);
					RV = false;
				}
			}
			if (0 == addedToFront || Name != RefText)
			{
				RefBox->addItem(Name);
				if (termT != 0 ? Name == termT->getName() : false)
				{
					SourceBox->setCurrentIndex(n);
					SV = false;
				}
			}
		}
		for (n=0, N = state->getNumFitDataSets(); n<N; n++)
		{
			Name = state->getFitDataName(n);
			if (0 == addedToFront || Name != SourceText)
			{
				SourceBox->addItem(Name);
				if (FitD != 0 ? Name == FitD->getName() : false)
				{
					SourceBox->setCurrentIndex(c+n);
					SV = false;
				}
			}
		}
		if (N>0) RefBox->addItem("Potential");
		for (m=0, c+=N, M = mol->getNumTransitions(); m<M; m++, c+=N) 
		{
			T = mol->getTransitionP(m);
			if (T->getLowerState() == state || T->getUpperState() == state)
				for (n=0, N = T->getNumLineTables(); n<N; n++)
			{
				Name = T->getLineTableName(n);
				if (0 == addedToFront || Name != SourceText)
				{
					SourceBox->addItem(Name);
					if (lineT != 0 ? Name == lineT->getName() : false)
					{
						SourceBox->setCurrentIndex(c+n);
						SV = false;
					}
				}
			}
		}
	}
	Block = false;
	if (RV) refSourceBoxChanged();
	else if (SV) sourceBoxChanged();
}

void ResidualPlot::vBoxChanged(int i)
{
	if (Block || i == -1) return;
	//printf("ResidualPlot::vBoxChanged\n");
	int oC = 0;
	Block = true;
	v = vBox->currentText().toInt();
	if (FSA)
	{
		QString CT = CompBox->currentText(), B;
		int c;
		CompBox->clear();
		for (c=0; c < NC; c++) if (nData[c][refTIsoT[Iso]][v] > 0)
		{
			CompBox->addItem(B = QString::number(c));
			if (B == CT) oC = CompBox->count() - 1;
		}
	}
	else
	{
		if (CompBox->currentText() == "f levels") oC = 1;
		CompBox->clear();
		if (nData[0][refTIsoT[Iso]][v] > 0) CompBox->addItem("e levels");
		else oC = 0;
		if (NC > 1 && nData[1][refTIsoT[Iso]][v] > 0) CompBox->addItem("f levels");
		else oC = 0;
		CompBox->setCurrentIndex(oC);
	}
	Block = false;
	compBoxChanged(oC);
	//printf("Ende vBoxChanged\n");
}

void ResidualPlot::SetPoints()
{
    DiagWindow::SetPoints();
    if (0 != m_ResidualFit)
    {
        int NSplines = m_ResidualFit->getNumSplines();
        for (int si = 0; si < NSplines; ++si)
        {
            Spline* currentSpline = m_ResidualFit->getSpline(si);
            numPoints += currentSpline->getNPoints();
        }
        points = new SplinePoint[numPoints];
        for (int si = 0, pi = 0; si < NSplines; ++si)
        {
            Spline* currentSpline = m_ResidualFit->getSpline(si);
            for (int spi = 0, N = currentSpline->getNPoints(); spi < N; ++spi) points[pi++] = currentSpline->getSplinePoint(spi);
        }
    }
}

void ResidualPlot::DetermineClickedSplineAndClickedPoint(const int i_overalIndex)
{
    if (0 != m_ResidualFit)
    {
        int N, si, NumSplines = m_ResidualFit->getNumSplines(), lastN = 0;
        for (si = N = 0; N < i_overalIndex && si < NumSplines; ++si)
        {
            lastN = N;
            N += m_ResidualFit->getSpline(si)->getNPoints();
        }
        if (si < NumSplines)
        {
            m_clickedSpline = si;
            m_clickedPoint = i_overalIndex - lastN;
        }
        else m_clickedSpline = m_clickedPoint = -1;
    }
}

void ResidualPlot::MovePoint(int i_pointIndex, double i_newX, double i_newY)
{
    if (0 <= m_clickedSpline && 0 <= m_clickedPoint && mPoint == i_pointIndex && 0 != m_ResidualFit
            && m_clickedSpline < m_ResidualFit->getNumSplines() && m_clickedPoint < m_ResidualFit->getSpline(m_clickedSpline)->getNPoints())
    {
        m_ResidualFit->getSpline(m_clickedSpline)->MovePoint(m_clickedPoint, i_newX, i_newY);
        m_ResidualFit->getSpline(m_clickedSpline)->Refit();
    }
}

void ResidualPlot::removePoint()
{
    DetermineClickedSplineAndClickedPoint(sPoint);
    if (0 != m_ResidualFit && m_clickedSpline >= 0 && m_clickedSpline < m_ResidualFit->getNumSplines())
    {
        Spline* spline = m_ResidualFit->getSpline(m_clickedSpline);
        if (spline != 0 && m_clickedPoint >= 0 && m_clickedPoint < spline->getNPoints()) spline->RemovePoint(m_clickedPoint);
    }
    sPoint = -1;
}

void ResidualPlot::addPoint()
{
    if (0 == m_ResidualFit) return;
    int i, N = m_ResidualFit->getNumSplines();
    double minDist = 1e99, cDiff, x0, xN;
    Spline* currentSpline, *bestSpline;
    for (i=0; i<N; ++i)
    {
        currentSpline = m_ResidualFit->getSpline(i);
        if ((x0 = currentSpline->getx0()) <= cPosx && (xN = currentSpline->getxN()) >= cPosx)
        {
            currentSpline->AddPoint(cPosx, cPosy);
            break;
        }
        if ((cDiff = x0 - cPosx) > 0.0 && cDiff < minDist)
        {
            minDist = cDiff;
            bestSpline = currentSpline;
        }
        else if ((cDiff = cPosx - xN) > 0.0 && cDiff < minDist)
        {
            minDist = cDiff;
            bestSpline = currentSpline;
        }
    }
    if (i == N) bestSpline->AddPoint(cPosx, cPosy);
    cPosx = cPosy = -1;
}
