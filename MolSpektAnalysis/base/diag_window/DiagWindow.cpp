//
// C++ Implementation: DiagWindow
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2006 - 2020
//
// Copyright: See README file that comes with this source code
//
//

#include "DiagWindow.h"
#include "Picture.h"
#include "MainWindow.h"
#include "tools.h"
#include "utils.h"
#include "naturalspline.h"
#include "datensatz.h"
#include "marker.h"
#include "SplinePoint.h"

#include <qvariant.h>
#include <qscrollbar.h>
#include <qlabel.h>
#include <qlineedit.h>
#include <qlayout.h>
#include <qtooltip.h>
#include <qaction.h>
#include <qmenubar.h>
#include <qimage.h>
#include <qpixmap.h>
#include <QFileDialog>
#include <qstring.h>
#include <qfile.h>
#include <qmessagebox.h>
#include <QTextStream>
#include <qstringlist.h>
#include <qpainter.h>
#include <qvalidator.h>
#include <QPrinter>
#include <QKeyEvent>
#include <QGridLayout>
#include <QMouseEvent>
#include <QPrintDialog>
#include <QImage>
#include <QToolButton>
#include <QImageWriter>

#include <cmath>
#include <limits.h>
#include <stdio.h>

using std::isinf;
using std::isnan;

DiagWindow::DiagWindow(Type type, MainWindow *mw, QString Filter, QString FileExt, int o)
    : MDIChild(type, mw, Filter, FileExt), HScroll(0), xStartLabel(0), xStart(0), xStopLabel(0), xStop(0), VScroll(0), yStartLabel(0), yStart(0)
    , yStopLabel(0), yStop(0), Bild(0), SpektrumLayout(0), spacer1(0), IsoBox(0), StateBox(0), JBox(0), vBox(0), SimBox(0), TUnit(0), MolBox(0)
    , SourceBox(0), RefBox(0), CompBox(0), Evs(0), Temp(0), LaserF(0), WindowS(0), WindowE(0), Calc(0), PrintB(0), m_contrastSlider(0)
    , m_intensitySlider(0), points(0), XMax(0.0), XMin(0.0), YMax(0.0), YMin(0.0), XO(0.0), YO(0.0), XSF(0.0), YSF(0.0), aPosx(0.0), aPosy(0.0), m_minSelectedFrequency(-1.0)
    , m_maxSelectedFrequency(-1.0), ScaleMaxLargeSteps(5), ScaleSmallStepDiv(2), ScaleNAllowedSteps(3), ScaleAllowedSteps(new int[3]), ScaleXHeight(0), ScaleYWidth(0)
    , ScaleMinYWidth(0), ScaleTopOffset(0), numPoints(0), XUnit("Frequency [cm^{-1}]"), YUnit("Intensity [arbitrary units]"), UnitFont()
    , ScaleFont(), Daten(new Datensatz[26]), DMarkers(false)
    , ShowAssignmentsOnTop(false), showPoints(false), marker(0), AnzahlMarker(0), nDatenS(0), cPosx(0), cPosy(0), sPoint(-1), mPoint(-1)
    , AssFont(), SpektColor(0, 0, 0), CopyColor(new QColor[26]), ZoomB(0), Image(0)
{
    //printf("DiagWindow::DiagWindow, MW=%d\n", mw);
	ShowMarkerLabels = true;
	setMouseTracking(true);
	ScaleAllowedSteps[0] = 1;
	ScaleAllowedSteps[1] = 2;
	ScaleAllowedSteps[2] = 5;
	//ScaleFont = QFont();
	CopyColor[0].setRgb(0, 0, 0);
	CopyColor[1].setRgb(255, 0, 0);
	CopyColor[2].setRgb(0, 255, 0);
	CopyColor[3].setRgb(0, 0, 255);
	CopyColor[4].setRgb(255, 255, 0);
	CopyColor[5].setRgb(255, 0, 255);
	CopyColor[6].setRgb(0, 255, 255);
	CopyColor[7].setRgb(127, 0, 0);
	CopyColor[8].setRgb(0, 127, 0);
	CopyColor[9].setRgb(0, 0, 127);
	CopyColor[10].setRgb(127, 127, 0);
	CopyColor[11].setRgb(127, 0, 127);
	CopyColor[12].setRgb(0, 127, 127);
	CopyColor[13].setRgb(127, 127, 127);
	CopyColor[14].setRgb(255, 127, 0);
	CopyColor[15].setRgb(255, 0, 127);
	CopyColor[16].setRgb(127, 255, 0);
	CopyColor[17].setRgb(127, 0, 255);
	CopyColor[18].setRgb(0, 127, 255);
	CopyColor[19].setRgb(0, 255, 127);
	CopyColor[20].setRgb(255, 127, 127);
	CopyColor[21].setRgb(127, 255, 127);
	CopyColor[22].setRgb(127, 127, 255);
	CopyColor[23].setRgb(255, 255, 127);
	CopyColor[24].setRgb(255, 127, 255);
	CopyColor[25].setRgb(127, 255, 255);
	int width = 500, height = 500;
	if (MW != 0)
	{
		if (MW->size().width() < 625) width = (4 * MW->size().width()) / 5;
		if (MW->size().height() < 625) height = (4 * MW->size().height()) / 5;
	}
	else printf("Error: MW == 0! Can't resize correctly!\n");
	resize(width, height);
	SpektrumLayout = new QGridLayout(this); 
	QLabel *Labelt2,  *Labelt1, *Labelw4,  *Labelw3,  *Labelw2, *Labelw1, *L1, *L2, *L3, *L, *SimL, *MolL, *TempL, *LaserL, *WindS, *WindE;
	QSpacerItem *spacert1;
	QGridLayout *Layout, *Layout2;
	switch (type)
	{
		case TermEnergyPlot:
			Layout = new QGridLayout;
			Labelt1 = new QLabel("Zuordnung v':", this);
			Layout->addWidget( Labelt1, 0, 0 );
			//Label1->setGeometry(10, 10, 100, 20);
			Evs = new QLineEdit("", this);
			Layout->addWidget(Evs, 0, 1);
			//vs->setGeometry(120, 10, 100, 20);
			spacert1 = new QSpacerItem( 16, 21, QSizePolicy::Expanding, QSizePolicy::Minimum );
			Layout->addItem( spacert1, 0, 2 );
			Labelt2 = new QLabel("Show isotopomer: ", this);
			Layout->addWidget(Labelt2, 0, 4);
			IsoBox = new QComboBox(this);
			Layout->addWidget(IsoBox, 0, 5);
			SpektrumLayout->addLayout(Layout, 0, 0, 1, 9);
			o++;
			break;
		case WaveFunctionPlot:
			Layout = new QGridLayout;
			Labelw1 = new QLabel("Potential:", this);
			Layout->addWidget(Labelw1, 0, 0);
			StateBox = new QComboBox(this);
			StateBox->setEditable(false);
			Layout->addWidget(StateBox, 0, 1);
			Labelw2 = new QLabel("Isotopomer:", this);
			Layout->addWidget(Labelw2, 0, 2);
			IsoBox = new QComboBox(this);
			IsoBox->setEditable(false);
			Layout->addWidget(IsoBox, 0, 3);
			Labelw3 = new QLabel("J:", this);
			Layout->addWidget(Labelw3, 0, 4);
			JBox = new QComboBox(this);
			Layout->addWidget(JBox, 0, 5);
			Labelw4 = new QLabel("v:", this);
			Layout->addWidget(Labelw4, 0, 6);
			vBox = new QComboBox(this);
			vBox->setEditable(false);
			Layout->addWidget(vBox, 0, 7);
			SpektrumLayout->addLayout(Layout, 0, 0, 1, 9);
			o++;
			break;
		case DataSetPlot:
			Layout = new QGridLayout;
			L1 = new QLabel("Molecule:", this);
			Layout->addWidget(L1, 0, 0);
			MolBox = new QComboBox(this);
			Layout->addWidget(MolBox, 0, 1);
			L2 = new QLabel("Electronic state:", this);
			Layout->addWidget(L2, 0, 2);
			StateBox = new QComboBox(this);
			Layout->addWidget(StateBox, 0, 3);
			L3 = new QLabel("Data source:", this);
			Layout->addWidget(L3, 0, 4);
			SourceBox = new QComboBox(this);
			Layout->addWidget(SourceBox, 0, 5);
			SpektrumLayout->addLayout(Layout, 0, 0, 1, 9);
			o++;
			break;
		case ResidPlot:
			setMinimumSize(800, 600);
			Layout = new QGridLayout;
			L = new QLabel("Molecule:", this);
			Layout->addWidget(L, 0, 0);
			MolBox = new QComboBox(this);
			Layout->addWidget(MolBox, 0, 1);
			Layout->setColumnStretch(1, 1);
			L = new QLabel("Electronic state:", this);
			Layout->addWidget(L, 0, 2);
			StateBox = new QComboBox(this);
			Layout->addWidget(StateBox, 0, 3);
			Layout->setColumnStretch(3, 1);
			L = new QLabel("Data source:", this);
			Layout->addWidget(L, 0, 4);
			SourceBox = new QComboBox(this);
			Layout->addWidget(SourceBox, 0, 5);
			Layout->setColumnStretch(5, 1);
			L = new QLabel("Reference termTable:", this);
			Layout->addWidget(L, 0, 6);
			RefBox = new QComboBox(this);
			Layout->addWidget(RefBox, 0, 7);
			Layout->setColumnStretch(7, 1);
			SpektrumLayout->addLayout(Layout, 0, 0, 1, 9);
			Layout = new QGridLayout;
			L = new QLabel("Isotopologue:", this);
			Layout->addWidget(L, 0, 0);
			IsoBox = new QComboBox(this);
			Layout->addWidget(IsoBox, 0, 1);
			Layout->setColumnStretch(1, 1);
			L = new QLabel("v:", this);
			Layout->addWidget(L, 0, 2);
			vBox = new QComboBox(this);
			Layout->addWidget(vBox, 0, 3);
			Layout->setColumnStretch(3, 1);
			L = new QLabel("Component:", this);
			Layout->addWidget(L, 0, 4);
			CompBox = new QComboBox(this);
			Layout->addWidget(CompBox, 0, 5);
			Layout->setColumnStretch(5, 1);
			PrintB = new QPushButton("Print all...");
			Layout->addWidget(PrintB, 0, 6);
			SpektrumLayout->addLayout(Layout, 1, 1, 1, 9);
			o+=2;
			break;
		case AddSpect:
			setMinimumSize(800, 600);
			Layout = new QGridLayout;
			Layout->addWidget(new QLabel("Source for energy offset:", this), 0, 0);
			Layout->addWidget(SourceBox = new QComboBox(this), 0, 1);
			Layout->addWidget(new QLabel("use:"), 0, 2);
			Layout->addWidget(CompBox = new QComboBox(this), 0, 3);
			SpektrumLayout->addLayout(Layout, 0, 0, 1, 9);
            Layout = new QGridLayout;
            Layout->addWidget(new QLabel("ln2 contrast:",  this), 0, 0);
            Layout->addWidget(m_contrastSlider = new QSlider(Qt::Horizontal, this), 0, 1);
            m_contrastSlider->setMinimum(-10);
            m_contrastSlider->setMaximum(10);
            Layout->setColumnStretch(1, 1);
            Layout->addWidget(new QLabel("intens offset:", this), 0, 2);
            Layout->addWidget(m_intensitySlider = new QSlider(Qt::Horizontal, this), 0, 3);
            m_intensitySlider->setMinimum(-256);
            m_intensitySlider->setMaximum(512);
            Layout->setColumnStretch(3, 1);
            SpektrumLayout->addLayout(Layout, 1, 0, 1, 9);
            o+=2;
			SourceBox->setEditable(false);
			CompBox->setEditable(false);
			break;
		case SpectrumSimulation:
			Layout = new QGridLayout;
			SimL = new QLabel("Simulate:", this);
			Layout->addWidget(SimL, 0, 0);
			SimBox = new QComboBox(this);
			SimBox->setEditable(false);
			Layout->addWidget(SimBox, 0, 1);
			MolL = new QLabel("Molecule:", this);
			Layout->addWidget(MolL, 0, 2);
			MolBox = new QComboBox(this);
			MolBox->setEditable(false);
			Layout->addWidget(MolBox, 0, 3);
			TempL = new QLabel("T:", this);
			Layout->addWidget(TempL, 0, 4);
			Temp = new QLineEdit("20", this);
			Layout->addWidget(Temp, 0, 5);
			TUnit = new QComboBox(this);
			TUnit->setEditable(false);
			Layout->addWidget(TUnit, 0, 6);
			Calc = new QPushButton("Calculate", this);
			Layout->addWidget(Calc, 0, 7);
			Layout2 = new QGridLayout;
			LaserL = new QLabel("Laser frequency [cm^-1]:", this);
			Layout2->addWidget(LaserL, 0, 0);
			LaserF = new QLineEdit("17000", this);
			Layout2->addWidget(LaserF, 0, 1);
			WindS = new QLabel("Window:", this);
			Layout2->addWidget(WindS, 0, 2);
			WindowS = new QLineEdit("16600", this);
			Layout2->addWidget(WindowS, 0, 3);
			WindE = new QLabel("to", this);
			Layout2->addWidget(WindE, 0, 4);
			WindowE = new QLineEdit("16601", this);
			Layout2->addWidget(WindowE, 0, 5);
			SpektrumLayout->addLayout(Layout, 0, 0, 1, 9);
			SpektrumLayout->addLayout(Layout2, 1, 1, 1, 9);
			o+=2;
			break;
		default:
			break;
	}
	HScroll = new QScrollBar(this);
    HScroll->setMaximum( 100000 );
    HScroll->setPageStep( 9 );
    HScroll->setOrientation( Qt::Horizontal );

    SpektrumLayout->addWidget( HScroll, o+2, 0, 1, 9 );

    xStartLabel = new QLabel("xStart", this);

    SpektrumLayout->addWidget( xStartLabel, o, 0 );

    xStart = new QLineEdit("0", this);

    SpektrumLayout->addWidget( xStart, o, 1 );

    xStopLabel = new QLabel("xStop", this);

    SpektrumLayout->addWidget( xStopLabel, o, 2 );

    xStop = new QLineEdit("1000", this);

    SpektrumLayout->addWidget( xStop, o, 3 );
    spacer1 = new QSpacerItem( 16, 21, QSizePolicy::Expanding, QSizePolicy::Minimum );
    SpektrumLayout->addItem( spacer1, o, 4 );

    VScroll = new QScrollBar(this);
    VScroll->setMaximum(100000);
    VScroll->setPageStep( 9 );
    VScroll->setOrientation( Qt::Vertical );

    SpektrumLayout->addWidget( VScroll, o+1, 9 );

    yStartLabel = new QLabel("yStart", this);

    SpektrumLayout->addWidget( yStartLabel, o, 5 );

    yStart = new QLineEdit("0", this);

    SpektrumLayout->addWidget( yStart, o, 6 );

    yStopLabel = new QLabel("yStop", this);

    SpektrumLayout->addWidget( yStopLabel, o, 7 );

    yStop = new QLineEdit("1000", this);

    SpektrumLayout->addWidget( yStop, o, 8 );
	
	ZoomB = new QToolButton(this);
	ZoomB->setToolButtonStyle(Qt::ToolButtonTextOnly);
	ZoomB->setCheckable(true);
	ZoomB->setText("Z");
	ZoomB->setChecked(true);
	SpektrumLayout->addWidget(ZoomB, o, 9);
	
    Bild = new Picture(this, "Bild" );
    Bild->setFocusPolicy( Qt::StrongFocus );

    SpektrumLayout->addWidget(Bild, o+1, 0, 1, 9);
	//SpektrumLauout->setRowStretch(
	
    //languageChange();
    //resize( QSize(1033, 882).expandedTo(minimumSizeHint()) );
    //clearWState( WState_Polished );
    
    connect( Bild, SIGNAL( Resized(QSize*) ), this, SLOT( Paint() ) );
    connect( Bild, SIGNAL( SelectionChanged(QRect*) ), this, SLOT( Zoom(QRect*) ) );
    connect( yStart, SIGNAL( returnPressed() ), this, SLOT( Paint() ) );
    connect( yStop, SIGNAL( returnPressed() ), this, SLOT( Paint() ) );
    connect( xStart, SIGNAL( returnPressed() ), this, SLOT( Paint() ) );
    connect( xStop, SIGNAL( returnPressed() ), this, SLOT( Paint() ) );
    connect( VScroll, SIGNAL( valueChanged(int) ), this, SLOT( VScrolled(int) ) );
    connect( HScroll, SIGNAL( valueChanged(int) ), this, SLOT( HScrolled(int) ) );
    connect( Bild, SIGNAL( KeyPressed(QKeyEvent*) ), this, SLOT( KeyPressed(QKeyEvent*) ) );
    connect( Bild, SIGNAL( RightClicked() ), this, SLOT( ZoomOut() ) );
    connect( Bild, SIGNAL( LeftClicked(QPoint*) ), this, SLOT( PictureClicked(QPoint*) ) );
	connect(Bild, SIGNAL(MouseMoved(int, int)), this, SLOT(MouseMovedOverPicture(int, int)));
    
    // tab order
    setTabOrder( Bild, xStart );
    setTabOrder( xStart, xStop );
    setTabOrder( xStop, yStart );
    setTabOrder( yStart, yStop );
    init();
	//printf("Ende DiagWindow, MW=%d\n", MW);
}

DiagWindow::~DiagWindow()
{
	//printf("DiagWindow::~DiagWindow\n");
	delete[] ScaleAllowedSteps;
	delete[] Daten;
	//printf("Ende ~DiagWindow\n");
}

void DiagWindow::ClearMarked()
{
	emit clearMarked();
}

void DiagWindow::languageChange()
{
    
}

void DiagWindow::init()
{
     VWNVon = new QDoubleValidator(xStart);
     xStart->setValidator(VWNVon);
     VWNBis = new QDoubleValidator(xStop);
     xStop->setValidator(VWNBis);
     VEVon = new QDoubleValidator(yStart);
     yStart->setValidator(VEVon);
     VEBis = new QDoubleValidator(yStop);
     yStop->setValidator(VEBis);
     //DebugPrint = false;
  }

void DiagWindow::destroy()
{

}

void DiagWindow::Print(QPrinter &P)
{
	QPrintDialog *PD = new QPrintDialog(&P, this);
    int i;
    //P.setMinMax(1, 1);
    if (!PD->exec()) return;
    QRect R(0, 0, P.width(), P.height());
	QPainter Pt(&P);
    for (i=0; i < P.copyCount(); i++)
    {
		if (i>0) while (!P.newPage()) 
		{
	    	if (!P.outputFileName().isEmpty()) break;
	    	if (QMessageBox::warning(this, "Fehler beim Drucken", "Neuer Versuch?",
						      QMessageBox::Retry, QMessageBox::Abort) ==
					 QMessageBox::Abort) return;
		}
		PaintScale(Pt, R);
        setClippingRect(Pt, P.width(), P.height());
		PSpektrum(Pt, R, false);
    }
    Pt.end();
	delete PD;
}

void DiagWindow::setClippingRect(QPainter& Pt, int width, int height)
{
    Pt.setClipping(true);
    Pt.setClipRect(ScaleYWidth, 0, width - ScaleYWidth, height - ScaleXHeight);
}

void DiagWindow::Paint()
{
     double xmax, xmin, ymax, ymin;
     Test(xmin, xmax, xStart, xStop);
     Test(ymin, ymax, yStart, yStop);
     ScrollsEnabled = false;
     double HRange = 10 * (XMax - XMin - xmax + xmin) / (xmax - xmin);
     double VRange = 10 * (YMax - YMin - ymax + ymin) / (ymax - ymin);
     //printf("YMax - YMin = %e | ymax - ymin = %e | VRange = %d\n", 
		//   YMax - YMin, ymax - ymin, Runden(VRange));
     HScroll->setRange(0, Runden(HRange));
     VScroll->setRange(0, Runden(VRange));
     int HValue = Runden(10 * (xmin - XMin) / (xmax - xmin));
     int VValue = Runden(VRange - 10 * (ymin - YMin) / (ymax - ymin));
     //printf("YMin=%e, ymin=%e, VValue=%d\n", YMin, ymin, VValue);
     HScroll->setValue(HValue);
     VScroll->setValue(VValue);
     ScrollsEnabled = true;
	 QPainter P(Bild->getPixmap());
     PaintScale(P, Bild->contentsRect());
     setClippingRect(P, Bild->width(), Bild->height());
     PSpektrum(P, Bild->contentsRect(), false);
     Bild->update();
}

void DiagWindow::PaintScale(QPainter &P, const QRect &R)
{
    //printf("DiagWindow::PaintScale\n");
	QFontMetrics SFM(ScaleFont);
	QString *XST = 0, *YST = 0;
	double *xls = 0, *yls = 0;
	double *XS = 0, *YS = 0;
	bool Success = false, XL = true, YL = true;
    double XM, YM, XStep, YStep;
    mXStart = xStart->text().toDouble();
    mYStart = yStart->text().toDouble();
    mXStop = xStop->text().toDouble();
    mYStop = yStop->text().toDouble();
    if (isinf(mXStart) || isnan(mXStart) || isinf(mXStop) || isnan(mXStop) || isinf(mYStart) || isnan(mYStop)
           || isnan(mYStart) || isinf(mYStop) || mXStop == mXStart || mYStart == mYStop)
	{
		printf("Diagwindow::PaintScale: error invalid boundaries: \n");
        printf("XStart=%f, XStop=%f, YStart=%f, YStop=%f\n", mXStart, mXStop, mYStart, mYStop);
		return;
	}
	//printf("XStart=%f, XStop=%f, YStart=%f, YStop=%f\n", XStart, XStop, YStart, YStop);
    double XWidth = mXStop - mXStart, YWidth = mYStop - mYStart;
	double XSS, YSS, xso, xlo, yso, ylo, W;
	int YMLS = ScaleMaxLargeSteps, XMLS = ScaleMaxLargeSteps, j, se, le, fe, w, rw, fh;
	int XSSD = ScaleSmallStepDiv, YSSD = ScaleSmallStepDiv, i, nsx, nlx, nsy, nly, nd, l, h, x, y;
	//printf("Vor loop\n");
	while (!Success)
	{
		XM = XWidth / XMLS;
		YM = YWidth / YMLS;
		//printf("XM=%f, XWidth=%f, XMLS=%d\n", XM, XWidth, XMLS);
		//printf("YM=%f, YWidth=%f, YMLS=%d\n", YM, YWidth, YMLS);
		for (XStep = 1.0; XStep < XM; XStep *= 10.0) ;
		for (YStep = 1.0; YStep < YM; YStep *= 10.0) ;
		//printf("1. XStep=%f\n", XStep);
		while (XStep > XM) XStep *= 0.1;
		while (YStep > YM) YStep *= 0.1;
		//printf("2. XStep=%f\n", XStep);
		for (i=0; (i < ScaleNAllowedSteps ? XStep * ScaleAllowedSteps[i] < XM : false); i++) ;
		if (i == ScaleNAllowedSteps) XStep *= 10;
		else XStep *= ScaleAllowedSteps[i];
		for (i=0; (i < ScaleNAllowedSteps ? YStep * ScaleAllowedSteps[i] < YM : false); i++) ;
		if (i == ScaleNAllowedSteps) YStep *= 10;
		else YStep *= ScaleAllowedSteps[i];
		//printf("1. XStep=%f\n", XStep);
		XSS = XStep / XSSD;
		YSS = YStep / YSSD;
		//printf("YSS=%f, YSSD=%d, YStep=%f\n", YSS, YSSD, YStep);
        xso = XSS * ceil(mXStart / XSS);
        xlo = XStep * ceil(mXStart / XStep);
        yso = YSS * ceil(mYStart / YSS);
        ylo = YStep * ceil(mYStart / YStep);
        nsx = floor((mXStop - xso) / XSS) + 1;
        nlx = floor((mXStop - xlo) / XStep) + 1;
        nsy = floor((mYStop - yso) / YSS) + 1;
        nly = floor((mYStop - ylo) / YStep) + 1;
		//printf("YSS=%f, YSSD=%d\n", YSS, YSSD);
		//printf("YStart=%f, YStep=%f, ylo=%f, YStop=%f\n", YStart, YStep, ylo, YStop);
		//printf("Nach Beginn, nsx=%d, nsy=%d, nlx=%d, nly=%d\n", nsx, nsy, nlx, nly);
		if (nlx == 0)
		{
			XMLS++;
			XL = false;
			continue;
		}
		if (nly == 0)
		{
			YMLS++;
			YL = false;
			continue;
		}
		//printf("nsy=%d\n", nsy);
		if (XS != 0) delete[] XS;
		if (YS != 0) delete[] YS;
		XS = new double[nsx];
		YS = new double[nsy];
		for (i=1, XS[0] = xso; i < nsx; i++) XS[i] = XS[i-1] + XSS;
		for (i=1, YS[0] = yso; i < nsy; i++) YS[i] = YS[i-1] + YSS;
		if (xls != 0) delete[] xls;
		if (yls != 0) delete[] yls;
		xls = new double[nlx];
		yls = new double[nly];
		for (i=1, xls[0] = xlo; i < nlx; i++) xls[i] = xls[i-1] + XStep;
		for (i=1, yls[0] = ylo; i < nly; i++) yls[i] = yls[i-1] + YStep;
		//printf("Nach XS,xls\n");
		if (XST != 0) delete[] XST;
		if (YST != 0) delete[] YST;
		XST = new QString[nlx];
		YST = new QString[nly];
		//printf("Nach Definition\n");
		nd = int(floor(log10(XStep)));
		if (nd > 0) nd = 0;
		else nd *= -1;
		for (i=0; i < nlx; i++)
			XST[i] = QString::number(xls[i], 'f', nd);
		ScaleXHeight = 15 + SFM.height() + (3 * TextHeight(UnitFont, XUnit)) / 2;
		nd = int(floor(log10(YStep)));
		if (nd > 0) nd = 0;
		else nd *= -1;
		for (i=0; i < nly; i++)
			YST[i] = QString::number(yls[i], 'f', nd);
        ScaleYWidth = SFM.horizontalAdvance(YST[0]);
        if (ScaleYWidth < SFM.horizontalAdvance(YST[nly - 1])) ScaleYWidth = SFM.horizontalAdvance(YST[nly - 1]);
		ScaleYWidth += 15 + 2 * TextHeight(UnitFont, YUnit);
		if (ScaleYWidth < ScaleMinYWidth) ScaleYWidth = ScaleMinYWidth;
		XSF = (R.width() - ScaleYWidth) / XWidth;
		YSF = (R.height() - ScaleXHeight - ScaleTopOffset) / YWidth;
		//printf("XSF=%f, R.width=%d, ScaleYWidth=%d, XWidth=%f\n", XSF, R.width(), ScaleYWidth, XWidth);
		//printf("YSF=%f, R.height=%d, ScaleXHeight=%d, ScaleTopOffset=%d, YWidth=%d\n",
			//   YSF, R.height(), ScaleXHeight, ScaleTopOffset, YWidth);
		if (nlx > 1 && XL)
		{
			W = XSF * (xls[1] - xls[0]);
			//printf("W=%f, xls[0]=%f, xls[1]=%f\n", W, xls[0], xls[1]);
			//printf("XST[0]=%s, XST[1]=%s, XST[%d]=%s, XST[%d]=%s\n", XST[0].ascii(), XST[1].ascii(),
				//   nlx-2, XST[nlx-2].ascii(), nlx-1, XST[nlx-1].ascii());
            if ((W < 0.5 * (SFM.horizontalAdvance(XST[0]) + SFM.horizontalAdvance(XST[1]))
                      || W < 0.5 * (SFM.horizontalAdvance(XST[nlx - 2]) + SFM.horizontalAdvance(XST[nlx - 1]))) && XMLS > 1)
				XMLS--;
			else
			{
				if ((XS[1] - XS[0]) * XSF < 2 && XSSD > 1)
				{
					if (XSSD == 5) XSSD = 2;
					else XSSD /= 2;
				}
				else Success = true;
			}
		}
		else Success = true;
		if (nly > 1 && YL)
		{
			if (YSF * (yls[1] - yls[0]) < 1.5 * SFM.height() && YMLS > 1)
			{
				//printf("YSF=%f, yls[1]=%f, yls[0]=%f, SFM.height()=%d\n", 
					//   YSF, yls[1], yls[0], SFM.height());
				YMLS--;
				Success = false;
			}
			else
			{
				if ((YS[1] - YS[0]) * YSF < 2 && YSSD > 1)
				{
					if (YSSD == 5) YSSD = 2;
					else YSSD /= 2;
					Success = false;
				}
			}
		}
		//printf("Nach Test\n");
		if (!Success) continue;
		l = ScaleYWidth - 1;
		h = R.height() - ScaleXHeight + 1;
        XO = double(ScaleYWidth) - XSF * mXStart;
        YO = double(h) - 1.0 + YSF * mYStart;
		P.fillRect(R.left(), R.top(), R.width(), R.height(), QColor(255, 255, 255));
		P.drawLine(l, 0.0, l, h + 19);
		P.drawLine(ScaleYWidth - 20, h, R.width(), h);
		rw = R.width();
		//printf("Vor X loop\n");
		for (i=0, se = h + 4, le = h + 9, fe = h + 14 + SFM.ascent(); i < nsx; i++)
		{
			x = int(XO + XSF * XS[i]);
			//printf("i=%d, nsx=%d, j=%d, XO=%d, XSF=%f, XS=%f\n", i, nsx, j, XO, XSF, XS[i]);
			//printf("x=%d, h=%d, se=%d, fe=%d, nlx=%d\n", x, h, se, fe, nlx);
			P.drawLine(x, h, x, se);
		}
		for (j=0; j < nlx; j++)
		{
			x = int(XO + XSF * xls[j]);
			P.drawLine(x, h, x, le);
            x -= (w = SFM.horizontalAdvance(XST[j])) / 2;
            if (x<=l ? (nlx > 1 ? l + w + 20 < XO + XSF * xls[1] - SFM.horizontalAdvance(XST[1]) / 2 : true)
					: false) 
				x = l + 10;
			else if (x + w >= rw ? (nlx > 1 ?
                            rw - w - 20 > XO + XSF * xls[nlx - 2] + SFM.horizontalAdvance(XST[nlx - 2]) / 2
							: true) : false)
								 x = rw - w - 10;
			if (x > l && x + w < rw) P.drawText(x, fe, XST[j]);
		}
		//printf("Vor y loop\n");
		fh = SFM.ascent();
		for (i=0, se = l - 4, le = l - 9, fe = l - 14; i < nsy; i++)
		{
			y = int(YO - YSF * YS[i]);
			P.drawLine(se, y, l, y);
		}
		for (j=0; j < nly; j++)
		{
			y = int(YO - YSF * yls[j]);
			P.drawLine(le, y, l, y);
			y += fh/2;
			if (y <= fh ? (nly > 1 ? fh < yls[nly - 2] * YSF - fh / 2 : true) : false) y = fh;
			else if (y>=h ? (nly > 1 ? h - fh > yls[1] + fh / 2 : true) : false) y=h;
            if (y<h && y > fh) P.drawText(fe - SFM.horizontalAdvance(YST[j]), y, YST[j]);
		}
		WriteText(P, l + (rw - l - TextWidth(UnitFont, XUnit)) / 2, 
				  R.height() - TextHeight(UnitFont, XUnit) / 2, XUnit, UnitFont, 0);
		WriteText(P, (5 * TextHeight(UnitFont, YUnit)) / 4, (h + TextWidth(UnitFont, YUnit)) / 2,
				  YUnit, UnitFont, 1);
		//printf("Vor Ende\n");
	}
	Bild->setMouseCrossBorders(ScaleYWidth, ScaleTopOffset + 1, 
							   Bild->width() - 2, Bild->height() - ScaleXHeight - 2);
	if (XS != 0) delete[] XS;
	if (YS != 0) delete[] YS;
	if (xls != 0) delete[] xls;
	if (yls != 0) delete[] yls;
	XST = new QString[nlx];
	YST = new QString[nly];
	//printf("Paintscale Ende\n");
}

void DiagWindow::drawPoint(QPainter &P, double X, double Y)
{
    if (X >= mXStart && X <= mXStop && Y >= mYStart && Y <= mYStop) P.drawPoint(static_cast<int>(XO + XSF * X), static_cast<int>(YO - YSF * Y));
}

void DiagWindow::startLine(double X, double Y)
{
    mLastX = X;
    mLastY = Y;
}

void DiagWindow::continueLine(QPainter &P, double X, double Y)
{
    int x1 = static_cast<int>(XO + XSF * mLastX);
    int y1 = static_cast<int>(YO - YSF * mLastY);
    int x2 = static_cast<int>(XO + XSF * X);
    int y2 = static_cast<int>(YO - YSF * Y);
    if (x1 != x2 || y1 != y2) P.drawLine(x1, y1, x2, y2);
    else P.drawPoint(X, Y);
    startLine(X, Y);
}

void DiagWindow::mouseMoveEvent(QMouseEvent *)
{
	if (MW != 0) MW->showStatusText();
	Bild->mouseOutside();
}

void DiagWindow::MouseMovedOverPicture(int x, int y)
{
	//printf("MouseMoved\n");
	if (MW == 0) return;
	if (x < ScaleYWidth || y > Bild->contentsRect().height() - ScaleXHeight)
	{
		MW->showStatusText();
		return;
	}
	double X = double(x - XO) / XSF, Y = -double(y - YO) / YSF;
    MW->showStatusText(getName() + ": " + XUnit + " = " + QString::number(X, 'g', 9)
			+ ", " + YUnit + " = " + QString::number(Y, 'g', 6));
}

void DiagWindow::PictureClicked(QPoint *P)
{
	//printf("DiagWindow::PictureClicked\n");
	emit clicked(double(P->x() - XO) / XSF, -double(P->y() - YO) / YSF);
}

void DiagWindow::Zoom( QRect *M )
{
    if (!ZoomB->isChecked()) return;
	QRect CR = Bild->contentsRect();
    double axmin = xStart->text().toDouble();
    double axmax = xStop->text().toDouble();
    double aymin = yStart->text().toDouble();
    double aymax = yStop->text().toDouble();
    //printf("aymin=%e, aymax=%e\n", aymin, aymax);
    double xm = CR.left() + ScaleYWidth;
    double ym = CR.top() + ScaleTopOffset;
    double xw = CR.width() - ScaleYWidth;
    double yw = CR.height() - ScaleXHeight - ScaleTopOffset;
    double nxi = M->left();
    double nyi = M->top();
    double nxa = M->right();
    double nya = M->bottom();
    double nxmin = axmin + (nxi - xm)  / xw  * (axmax - axmin);
    double nxmax = axmin + (nxa - xm) / xw * (axmax - axmin);
    double nymin = aymax + (ym - nya) / yw * (aymax - aymin);
    double nymax = aymax + (ym - nyi) / yw * (aymax - aymin);
	printf("nymin=%f\n", nymin);
    xStart->setText(QString::number(nxmin, 'g', 11));
    xStop->setText(QString::number(nxmax, 'g', 11));
    yStart->setText(QString::number(nymin, 'g', 5));
    yStop->setText(QString::number(nymax, 'g', 5));
    Paint();
}


void DiagWindow::HScrolled(int V)
{
    if (!ScrollsEnabled) return;
    //printf("HScrolled, V=%d\n", V);
    double axmin = xStart->text().toDouble();
    double axmax = xStop->text().toDouble();
    double d = axmax - axmin;
    double nxmin = XMin + V * d * 0.1;
    xStart->setText(QString::number(nxmin, 'g', 11));
    xStop->setText(QString::number(nxmin + d, 'g', 11));
    Paint();
}

void DiagWindow::VScrolled( int V )
{
    if (!ScrollsEnabled) return;
    //printf("VScrolled, V=%d\n", V);
    double aymin = yStart->text().toDouble();
    double aymax = yStop->text().toDouble();
    //printf("V=%d; aymin=%e; aymax=%e\n", V, aymin, aymax);
    double d = aymax - aymin;
    double nymin = YMin + (VScroll->maximum() - V) * d * 0.1;
    //printf("d=%e; nymin=%e\n", d, nymin);
    yStart->setText(QString::number(nymin, 'g', 7));
    yStop->setText(QString::number(nymin + d, 'g', 7));
    Paint();
}

void DiagWindow::setData(double **Data, int numRows)
{
	int r;
	double Min = 0.0, Max = 0.0, d;
	Daten->reinit();
	for (r=0; r < numRows; r++) 
	{
		Daten->AddValue(Data[r][0], Data[r][1], false);
		if (Data[r][1] > Max) Max = Data[r][1];
		if (Data[r][1] < Min) Min = Data[r][1];
	}
	d = (Max - Min) * 0.05;
	XMin = Data[0][0];
	XMax = Data[numRows - 1][0];
	YMin = Min - d;
	YMax = Max + d;
	if (mRescaleOnSetAndAdd)
	{
		xStart->setText(QString::number(XMin, 'g', 11));
		xStop->setText(QString::number(XMax, 'g', 11));
		yStart->setText(QString::number(YMin));
		yStop->setText(QString::number(YMax));
	}
	nDatenS = 1;
	Paint();
	Changed();
	setNewCreated();
}

void DiagWindow::setUnits(QString X, QString Y)
{
	XUnit = X;
	YUnit = Y;
	xStartLabel->setText(X + ": from");
	xStopLabel->setText("to");
	yStartLabel->setText(Y + ": from");
	yStopLabel->setText("to");
}

void DiagWindow::addData(double** Data, int numRows)
{
	if (nDatenS == 26)
	{
		QMessageBox::information(this, "MolSpektAnalysis", 
								 "The maximum number of functions has been reached!");
		return;
	}
	int r;
	if (nDatenS == 0) YMax = YMin = 0.0;
	double d = (YMax - YMin) / 22.0; 
	double Min = YMin + d, Max = YMax - d;
	Daten[nDatenS].reinit();
	for (r=0; r < numRows; r++)
	{
		Daten[nDatenS].AddValue(Data[r][0], Data[r][1], false);
		if (Data[r][1] > Max) Max = Data[r][1];
		if (Data[r][1] < Min) Min = Data[r][1];
	}
	nDatenS++;
	d = (Max - Min) * 0.05;
	if (Data[0][0] < XMin) XMin = Data[0][0];
	if (Data[numRows - 1][0] > XMax) XMax = Data[numRows - 1][0];
	YMin = Min - d;
	YMax = Max + d;
	if (mRescaleOnSetAndAdd)
	{
		xStart->setText(QString::number(XMin, 'g', 11));
		xStop->setText(QString::number(XMax, 'g', 11));
		yStart->setText(QString::number(YMin));
		yStop->setText(QString::number(YMax));
	}
	Paint();
	Changed();
}

void DiagWindow::clearData()
{
	nDatenS = 1;
	Daten->reinit();
}

void DiagWindow::clear()
{
	int n;
	for (n=0; n < nDatenS; n++) Daten[n].reinit();
	nDatenS = 0;
}

bool DiagWindow::writeData(QString FileName)
{
	int r, N = Daten->GetDSL();
	QFile Datei(FileName);
	if (!write(&Datei)) return false;
	QTextStream S(&Datei);
	for (r=0; r<N; r++) 
		S << Daten->GetValue(r, 0) << "\t" << Daten->GetValue(r, 1) 
				<< (Daten->GetMarked(r) ? "\tM\n" : "\n");
	Saved();
	return true;
}

void DiagWindow::getData(double **&Data, int &N)
{
	if ((N = Daten->GetDSL()) == 0)
	{
		Data = 0;
		return;
	}
	Data = Create(N, 2);
	int n;
	for (n=0; n<N; n++)
	{
		Data[n][0] = Daten->GetValue(n, 0);
		Data[n][1] = Daten->GetValue(n, 1);
	}
}

void DiagWindow::ZoomOut()
{
    if (!ZoomB->isChecked()) return;
	double xmin = xStart->text().toDouble();
    double xmax = xStop->text().toDouble();
    double ymin = yStart->text().toDouble();
    double ymax = yStop->text().toDouble();
    double dx = 0.5 * (xmax - xmin);
    double dy = 0.5 * (ymax - ymin);
    xmin -= dx;
    if (xmin < XMin) xmin = XMin;
    xmax += dx;
    if (xmax > XMax) xmax = XMax;
    ymin -= dy;
    if (ymin < YMin) ymin = YMin;
    ymax += dy;
    if (ymax > YMax) ymax = YMax;
    xStart->setText(QString::number(xmin, 'g', 11));
    xStop->setText(QString::number(xmax, 'g', 11));
    yStart->setText(QString::number(ymin, 'g', 7));
    yStop->setText(QString::number(ymax, 'g', 7));
    Paint();
}

void DiagWindow::setShowMouseCross(bool show)
{
	Bild->setMouseCross(show);
	if (show) 
		Bild->setMouseCrossBorders(ScaleYWidth, ScaleTopOffset + 1, 
								   Bild->width() - 2, Bild->height() - ScaleXHeight - 2);
}

bool DiagWindow::getShowMouseCross() const
{
	return Bild->getMouseCross();
}

void DiagWindow::KeyPressed( QKeyEvent *K )
{
    int V;
    switch (K->key())
    {
    case Qt::Key_Left:
	V = HScroll->value();
	if (V > 0) HScroll->setValue(V - 1);
	break;
    case Qt::Key_Right:
	V = HScroll->value();
	if (V < HScroll->maximum()) HScroll->setValue(V + 1);
	break;
    case Qt::Key_Up:
	V = VScroll->value();
	if (V > 0) VScroll->setValue(V - 1);
	break;
    case Qt::Key_Down:
	V = VScroll->value();
	if (V < VScroll->maximum()) VScroll->setValue(V + 1);
	break;
    case Qt::Key_Minus:
	ZoomOut();
	break;
    case Qt::Key_Plus:
	ZoomIn();
	break;
    default:
	KeyPressEvent(K);
	break;
    }
}

void DiagWindow::KeyPressEvent(QKeyEvent *K)
{
    K->ignore();
}

void DiagWindow::ZoomIn()
{
    if (!ZoomB->isChecked()) return;
	double xmin = xStart->text().toDouble();
    double xmax = xStop->text().toDouble();
    double ymin = yStart->text().toDouble();
    double ymax = yStop->text().toDouble();
    double dx = 0.25 * (xmax - xmin);
    double dy = 0.25 * (ymax - ymin);
    xmin += dx;
    xmax -= dx;
    ymin += dy;
    ymax -= dy;
    xStart->setText(QString::number(xmin, 'g', 11));
    xStop->setText(QString::number(xmax, 'g', 11));
    yStart->setText(QString::number(ymin, 'g', 7));
    yStop->setText(QString::number(ymax, 'g', 7));
    Paint();
}

void DiagWindow::PSpektrum(QPainter &P, const QRect &A, bool PrintFN )
{
    //printf("DiagWindow::PSpektrum\n");
	double xmin = xStart->text().toDouble(), xmax = xStop->text().toDouble();
	double ymin = yStart->text().toDouble(), ymax = yStop->text().toDouble();
	if (Image != 0)
	{
		double xsc = Image->width() / (XMax - XMin);
		double ysc = Image->height() / (YMax - YMin);
		P.drawImage(QRectF(A.left() + ScaleYWidth - 1, A.top() + ScaleTopOffset, A.width() - ScaleYWidth + 1, 
						   A.height() - ScaleXHeight), Image->mirrored(),
					QRectF((xmin - XMin) * xsc, Image->height() - (ymax - YMin) * ysc,
						   (xmax - xmin) * xsc, (ymax - ymin) * ysc));
	}
    if (nDatenS == 0) return;
	//printf("Nach Painter\n");
    int i, n, l;
    Marker *M;
	bool Draw = true;
	double lmin, lmax, lx, ly, ax, ay;
    int BIndex, DIndex, lBIndex;
	QString Text;
    for (i=0; i < AnzahlMarker; i++) marker[i].Visible = false;
    int Bxmin = A.left() + ScaleYWidth - 1;
    int Bymin = A.top() + ScaleTopOffset;
	int BHeight = A.height() - ScaleTopOffset - ScaleXHeight;
	int BWidth = A.width() - ScaleYWidth;
	int BTop = Bymin + BHeight, BRight = Bxmin + BWidth;
	int x1, y1, x2, y2, DStart;
    double LMLX = 0.0, LMLY = 0.0;
    int LMLB = 0;
    Marker *NM = NULL;
	//printf("Vor Daten, l=%d\n", l);
    for (n=1; true; n++)
	{
		if (n == nDatenS) n=0;
		if ((l = Daten->GetDSL()) == 0) 
		{
			if (n==0) return;
			continue;
		}
		printf("PSpektrum: n=%d, l=%d\n", n, l);
		P.setPen(CopyColor[n]);
		DStart = int(XO + XSF * Daten[n].GetValue(0, 0));
		for (DIndex = 0, lBIndex = -1; Daten[n].GetValue(DIndex,0) < xmin && DIndex < l - 1; DIndex++) ;
		if (DIndex < 1) DIndex = 1;
		lx = Daten[n].GetValue(DIndex - 1, 0);
		ly = Daten[n].GetValue(DIndex - 1, 1);
		//printf("Vor der Hauptschleife\n");
		for (BIndex = (Bxmin > DStart ? Bxmin : DStart); BIndex <= BRight && DIndex < l; BIndex++)
		{
			//printf("BIndex=%d, Bxmin=%d, BRight=%d\n", BIndex, Bxmin, BRight);
			//if (BIndex==685) printf("BIndex=685, l=%d, DIndex=%d\n", l, DIndex);
			ax = Daten[n].GetValue(DIndex, 0);
			if (XO + ax * XSF < BIndex + 1 || BIndex == BRight)
			{
				//if (BIndex==334) printf("Vor dem ersten Zeichnen.\n");
				ay = Daten[n].GetValue(DIndex, 1);
				lmin = lmax = ay;
				if (DMarkers && !ShowAssignmentsOnTop && n==0 && ShowMarkerLabels)
				{
					if (Daten->GetMarked(DIndex)) 
					{
						M = Daten->GetMarker(DIndex);
						if (lmin > LMLY)
						{
							LMLX = ax;
							LMLY = lmin;
							LMLB = BIndex;
							NM = M;
						}
						if (M->Marked)
						{
							M->Visible = true;
							M->x1 = BIndex - 4;
							M->y1 = (int)(YO - lmin * YSF) - 10;
							M->x2 = M->x1 + 8;
							M->y2 = M->y2 - 70;
						}
					}
					//if (BIndex==334) printf("Mitte DMarker\n");
					if (BIndex - LMLB > 8 && LMLB > 6)
					{
						NM->Visible = true;
						//if (BIndex==334) printf("NM=%d, AnzahlMarker=%d\n", NM - marker, AnzahlMarker);
						NM->x1 = LMLB - 4;
						NM->y1 = (int)(YO - YSF * LMLY) - 10;
						NM->x2 = NM->x1 + 8;
						NM->y2 = NM->y1 - 70;
						Text = QString::number(LMLX, 'g', 11);
						//if (BIndex==334) printf("Vor if\n");
						if (NM->DisplayData && !NM->Marked)
						{
							//if (BIndex==334) printf("Beginn if\n");
							Text += " " + NM->IsoName;
							if (NM->vss != -1) Text +=  " v''=" + QString::number(NM->vss);
							if (NM->Jss != -1) Text +=  " J''=" + QString::number(NM->Jss);
							if (NM->vs != -1) Text += " v'=" + QString::number(NM->vs);
							if (NM->Js != -1) Text += " J'=" + QString::number(NM->Js);
							if (NM->FC != -1) Text += " F_" + QString::number(NM->FC + 1);
							Text += " Diff=" + QString::number(NM->DD, 'g', 6);
							/*if (BIndex==334) 
							{
								printf("Vor NM->Line\n");
								printf("NM->Line[0]=%f\n", NM->Line[0]);
							}*/
							//if (ELU != NULL) Text += " OT=" 
							//	+ QString::number(NM->Line[0] + ELU[NM->Iso][NM->vss][NM->Jss], 'g', 11);
							P.setPen(QColor(0, 0, 255));
						}
						//if (BIndex==334) printf("Nach if\n");
						if (NM->y1 < BTop) WriteText(P, NM->x2, NM->y1, Text, ScaleFont, 1);
						//printf("NM->x2=%d, NM->y1=%d\n", NM->x2, NM->y1);
						LMLX = LMLY = 0.0;
						LMLB = 0;
					}
				}
				//if (BIndex==334) printf("Mitte 2. Zeichnen\n");
				if (DIndex > 0) 
				{
					x1 = lBIndex;
					y1 = (int)(YO - YSF * ly);
					x2 = BIndex; 
					y2 = (int)(YO - YSF * ay);
					if (lBIndex < Bxmin)
					{
						x1 = Bxmin;
						lBIndex = int(XO + XSF * lx);
						if (BIndex != lBIndex) 
							y1 += (y2 - y1) * (Bxmin - lBIndex) / (BIndex - lBIndex);
					}
					else if (BIndex == BRight)
					{
						x2 = int(XO + XSF * ax);
						if (x2 > A.right())
						{
							y2 -= (y2 - y1) * (x2 - A.right()) / (x2 - x1);
							x2 = A.right();
						}
					} 
					if (y1 < Bymin)
					{
						if (y2 < Bymin) Draw = false;
						else
						{
							x1 += (x2 - x1) * (Bymin - y1) / (y2 - y1);
							y1 = Bymin;
						}
					}
					else if (y1 > BTop)
					{
						if (y2 > BTop) Draw = false;
						else
						{
							x1 += (x2 - x1) * (y1 - BTop) / (y1 - y2);
							y1 = BTop;
						}
					}
					if (Draw)
					{
						if (y2 < Bymin)
						{
							x2 -= (x2 - x1) * (Bymin - y2) / (y1 - y2);
							y2 = Bymin;
						}
						else if (y2 > BTop)
						{
							x2 -= (x2 - x1) * (y2 - BTop) / (y2 - y1);
							y2 = BTop;
						}
                        if (ax >= m_minSelectedFrequency && ax <= m_maxSelectedFrequency) P.setPen(QColor(255, 0, 0));
                        else P.setPen(CopyColor[n]);
						P.drawLine(x1, y1 , x2, y2);
					}
					else Draw = true;
				}
				lBIndex = BIndex;
				lmin = lmax = Daten[n].GetValue(DIndex, 1);
				for (DIndex++; (DIndex < l ? XO + XSF * Daten[n].GetValue(DIndex,0) < BIndex + 1 : false);
					   DIndex++) 
				{
					ax = Daten[n].GetValue(DIndex, 0);
					ay = Daten[n].GetValue(DIndex, 1);
					if (ay < lmin) lmin = ay;
					if (ay > lmax) lmax = ay;
					if (DMarkers && Daten[n].GetMarked(DIndex) && !ShowAssignmentsOnTop && n==0) 
					{
						M = Daten->GetMarker(DIndex);
						if (ay > LMLY)
						{
							LMLX = M->Line[0];
							LMLY = ay;
							LMLB = BIndex;
							NM = M;
						}
						if (M->Marked)
						{
							M->Visible = true;
							M->x1 = BIndex - 4;
							M->y1 = (int)(YO - YSF * M->Line[1]) - 10;
							M->x2 = M->x1 + 8;
							M->y2 = M->y2 - 70;
						}
					}
				}
				//if (BIndex==334) printf("Vor dem zweiten Zeichnen.\n");
				//DIndex--;
				if (lmin != lmax)
				{
					y1 = (int)(YO - YSF * lmin);
					y2 = (int)(YO - YSF * lmax);
					if (y1 > Bymin && y2 < BTop)
					{
						if (y2 < Bymin) y2 = Bymin;
						if (y1 > BTop) y1 = BTop;
						P.drawLine(BIndex, y1, BIndex, y2);
						//printf("BIndex=%d, Bxmin=%d, y1=%d, y2=%d, width=%d, height=%d\n", 
							//   BIndex, Bxmin, y1, y2, A.width(), A.height());
					}
				}
				lx = ax;
				ly = ay;
			}
		}
		if (n==0) break;
	}
    if (DMarkers && !ShowAssignmentsOnTop)
    {
		P.setPen(QColor(255, 0, 0));
		//P.rotate(270);
		for (i=0; i<AnzahlMarker; i++) if (marker[i].Marked && marker[i].Visible)
		{
		    Text = QString::number(marker[i].Line[0], 'g', 11);
	    	if (marker[i].DisplayData)
		    {
				Text += marker[i].IsoName + " v''=";
				Text += QString::number(marker[i].vss) + " J''=";
				Text += QString::number(marker[i].Jss);
				if (marker[i].vs != -1)
				{
		    		Text += " v'=";
		    		Text += QString::number(marker[i].vs);
				}
				Text += " J'=";
				Text += QString::number(marker[i].Js);
				if (marker[i].FC != -1) Text += " F_" + QString::number(marker[i].FC + 1);
				Text += " Diff=";
				Text += QString::number(marker[i].DD, 'g', 6);
				//if (ELU != NULL) Text += " OT=" + QString::number(marker[i].Line[0] 
					//	+ ELU[marker[i].Iso][marker[i].vss][marker[i].Jss], 'g', 11);
				if (marker[i].oc != 0.0) Text += " o-c=" + QString::number(marker[i].oc, 'g', 6);
	    	}
			else if (marker[i].DD != -1.0) Text += " dist=" + QString::number(marker[i].DD, 'g', 9);
	    	if (marker[i].y1 < BTop) WriteText(P, marker[i].x1 + 8, marker[i].y1, Text, QFont(), 1);
		}
		//P.rotate(90);
		P.setPen(QColor(0, 0, 0));
    }
    if (PrintFN)
    {
		QString FN = windowTitle(); 
		FN = FN.right(FN.length() - FN.indexOf(' '));
		P.drawText(10, 10, FN);
    }
	//printf("Ende von DiagWindow::PSpektrum\n");
}

void DiagWindow::setImage(QImage* nImage)
{
	Image = nImage;
	Paint();
}

void DiagWindow::setRanges(double xmin, double xmax, double ymin, double ymax)
{
	XMin = xmin;
	XMax = xmax;
	YMin = ymin;
	YMax = ymax;
    setCurrentZoomRange(xmin, xmax, ymin, ymax);
}

void DiagWindow::setCurrentZoomRange(const double xmin, const double xmax, const double ymin, const double ymax)
{
    xStart->setText(QString::number(xmin, 'g', 6));
    xStop->setText(QString::number(xmax, 'g', 6));
    yStart->setText(QString::number(ymin, 'g', 6));
    yStop->setText(QString::number(ymax, 'g', 6));
    Paint();
}

bool DiagWindow::getShowMarkerLabels()
{
	return ShowMarkerLabels;
}

void DiagWindow::setShowMarkerLabels(bool Show)
{
	ShowMarkerLabels = Show;
	Paint();
}

void DiagWindow::exportPicture()
{
	QImageWriter IW("test");
	int n;
	QList<QByteArray> Formats = IW.supportedImageFormats();
	QString Filter, rFilter, FileName;
	//printf("exportPicture, n=%d\n", Formats.count());
	for (n=0; n < Formats.count(); n++) 
	{
		//printf(Formats[n].data());
		//printf("\n");
		IW.setFormat(Formats[n]);
		if (IW.canWrite())
		{
			if (!Filter.isEmpty()) Filter += ";;";
			Filter += Formats[n];
			Filter += " (*.";
			Filter += Formats[n];
			Filter += ")";
		}
	}
	FileName = QFileDialog::getSaveFileName(this, "Select file name and format",
											MW->getDir(Pict), Filter, &rFilter);
	if (FileName.isEmpty()) return;
	if (!Bild->savePicture(FileName, rFilter.left(rFilter.lastIndexOf(' '))))
		QMessageBox::information(this, "MolSpektAnalysis", 
								 "Error: the picture could not be saved!");
	MW->setDir(FileName.left(FileName.lastIndexOf(QRegExp("[\\/]"))), Pict);
}

void DiagWindow::SetPoints()
{
    if (points != 0) delete[] points;
    points = 0;
    numPoints = 0;
}

void DiagWindow::setShowPoints(bool show)
{
    if (show)
    {
        SetPoints();
        connect(Bild, SIGNAL(mouseMoved(QMouseEvent*)), this, SLOT(mouseMoved(QMouseEvent*)));
        connect(Bild, SIGNAL(mousePressed(QMouseEvent*)), this, SLOT(mousePressed(QMouseEvent*)));
        connect(Bild, SIGNAL(mouseReleased(QMouseEvent*)), this, SLOT(mouseReleased(QMouseEvent*)));
    }
    else if (showPoints)
    {
        if (points != 0)
        {
            delete[] points;
            points = 0;
            numPoints = 0;
        }
        disconnect(Bild, SIGNAL(mouseMoved(QMouseEvent*)), this, SLOT(mouseMoved(QMouseEvent*)));
        disconnect(Bild, SIGNAL(mousePressed(QMouseEvent*)), this, SLOT(mousePressed(QMouseEvent*)));
        disconnect(Bild, SIGNAL(mouseReleased(QMouseEvent*)), this, SLOT(mouseReleased(QMouseEvent*)));
    }
    showPoints = show;
    Paint();
}

void DiagWindow::mouseMoved(QMouseEvent *e)
{
    if (points != 0 && sPoint >= 0 && cPosx != 0 && cPosy != 0)
    {
        e->setAccepted(true);
        points[sPoint].x = aPosx + double(e->x() - cPosx) / XSF;
        points[sPoint].y = aPosy + double(cPosy - e->y()) / YSF;
        MovePoint();
        Paint();
    }
}

void DiagWindow::mousePressed(QMouseEvent *e)
{
    if (!showPoints || points == 0) return;
    int mx = e->x(), my = e->y(), px, py, n;
    //printf("mousePressed, mx=%d, my=%d\n", mx, my);
    sPoint = -1;
    for (n=0; n < numPoints; n++)
    {
        px = XO + points[n].x * XSF;
        py = YO - points[n].y * YSF;
        //printf("n=%d, px=%d, py=%d\n", n, px , py);
        if (mx >= px - 2 && mx <= px + 2 && my >= py - 2 && my <= py + 2) sPoint = n;
    }
    if (e->button() == Qt::RightButton)
    {
        if (sPoint == -1)
        {
            cPosx = mx;
            cPosy = my;
        }
        ShowPopupMenu(e->globalPos());
        e->setAccepted(true);
    }
    else if (e->button() == Qt::LeftButton && sPoint >= 0)
    {
        cPosx = mx;
        cPosy = my;
        mPoint = sPoint;
        aPosx = points[sPoint].x;
        aPosy = points[sPoint].y;
        e->setAccepted(true);
        HandleHistoryWhileMoving();
    }
}

void DiagWindow::mouseReleased(QMouseEvent *e)
{
    if (points != 0 && sPoint >= 0 && cPosx != 0 && cPosy != 0)
    {
        //printf("PotentialPlot::mouseReleased\n");
        double x = aPosx + double(e->x() - cPosx) / XSF;
        double y = aPosy + double(cPosy - e->y()) / YSF;
        HandleHistoryAfterMoving();
        MovePoint(sPoint, x, y);
        e->setAccepted(true);
        sPoint = -1;
        cPosx = cPosy = 0;
    }
}

void DiagWindow::DrawPoints(QPainter &P, const QRect &A)
{
    if (showPoints && points != 0)
    {
        int n, x, y;
        int l = A.left() + ScaleYWidth - 1, w = A.width() - ScaleYWidth;
        int t = A.top(), r = l + w, b = A.bottom() - ScaleXHeight;
        P.setPen(QColor(0, 0, 0));
        for (n=0; n < numPoints; n++)
        {
            x = XO + XSF * points[n].x;
            y = YO - YSF * points[n].y;
            //printf("n=%d, x=%f, y=%f\n", n, points[n].x, points[n].y);
            if (x > l && x < r && y > t && y < b) P.fillRect(x-2, y-2, 5, 5, QColor(0, 0, 0));
        }
    }
}
