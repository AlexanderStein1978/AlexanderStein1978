//
// C++ Interface: DiagWindow
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2006 - 2016
//
// Copyright: See README file that comes with this source code
//
//

#ifndef DIAGWINDOW_H
#define DIAGWINDOW_H

#include <QVariant>
#include <QPixmap>
#include <QLabel>
#include <QLayout>
#include <QValidator>
#include <QScrollBar>
#include <QGridLayout>
#include <QKeyEvent>
#include <QComboBox>
#include <QPushButton>
#include "tools.h"
#include "Picture.h"
#include "mdichild.h"

class MainWindow;
class Datensatz;
class Marker;

class QPainter;
class QPrinter;
class QFont;
class QMouseEvent;
class QImage;
class QToolButton;

struct SplinePoint;

class DiagWindow : public MDIChild
{
    Q_OBJECT

public:
    DiagWindow(Type type = SimpleDiagWindow, MainWindow *MW = 0, 
               QString Filter = "Data files (*.dat)", QString FileExt = ".dat", int o=0);
    virtual ~DiagWindow();
	virtual void setData(double **Data, int numRows);
	void addData(double **Data, int numRows);
	void clearData();
	void clear();
	void setShowMouseCross(bool show);
	bool getShowMouseCross() const;
	void setUnits(QString XUnit, QString YUnit);
	void setImage(QImage *Image);
	void setRanges(double xmin, double xmax, double ymin, double ymax);
    void setCurrentZoomRange(const double xmin, const double xmax, const double ymin, const double ymax);
	void setShowMarkerLabels(bool Show);
	bool getShowMarkerLabels();
	void exportPicture();
	virtual void ClearMarked();
    void setShowPoints(bool show);

    inline bool getShowPoints()
    {
        return showPoints;
    }
	
public slots:
	virtual void Print(QPrinter &P);
    virtual void Paint();
    virtual void Zoom( QRect * M );
    virtual void HScrolled( int V );
    virtual void VScrolled( int V );
    virtual void ZoomOut();
    virtual void KeyPressed( QKeyEvent * K );
    virtual void KeyPressEvent(QKeyEvent *K);
    virtual void ZoomIn();
    virtual void PictureClicked( QPoint * P );
	bool writeData(QString FileName = "");
	void getData(double **&Data, int &NPoints);

signals:
    void quit();
	void clicked(double x, double y);
	void clearMarked();

protected:
	void mouseMoveEvent(QMouseEvent *event);

    virtual void SetPoints();
    virtual void MovePoint(){}
    virtual void MovePoint(int /*i_pointIndex*/, double /*i_newX*/, double /*i_newY*/){}
    virtual void ShowPopupMenu(const QPoint& /*i_point*/){}
    virtual void HandleHistoryWhileMoving(){}
    virtual void HandleHistoryAfterMoving(){}
    virtual void PaintScale(QPainter &P, const QRect &A);
    virtual void PSpektrum(QPainter &P, const QRect &A, bool PrintFN );
    void DrawPoints(QPainter &P, const QRect &A);
		
	QScrollBar* HScroll;
    QLabel* xStartLabel;
    QLineEdit* xStart;
    QLabel* xStopLabel;
    QLineEdit* xStop;
    QScrollBar* VScroll;
    QLabel* yStartLabel;
    QLineEdit* yStart;
    QLabel* yStopLabel;
    QLineEdit* yStop;
    Picture* Bild;
    QGridLayout* SpektrumLayout;
    QSpacerItem* spacer1;
	QComboBox *IsoBox, *StateBox, *JBox, *vBox, *SimBox, *TUnit, *MolBox, *SourceBox, *RefBox, *CompBox;
    QLineEdit *Evs, *Temp, *LaserF, *WindowS, *WindowE;
	QPushButton *Calc, *PrintB;
    QSlider *m_contrastSlider, *m_intensitySlider;
    SplinePoint *points;
    double XMax, XMin, YMax, YMin, XO, YO, XSF, YSF, aPosx, aPosy, m_minSelectedFrequency, m_maxSelectedFrequency;
	int ScaleMaxLargeSteps, ScaleSmallStepDiv, ScaleNAllowedSteps, *ScaleAllowedSteps;
    int ScaleXHeight, ScaleYWidth, ScaleMinYWidth, ScaleTopOffset, numPoints;
	QString XUnit, YUnit;
	QFont UnitFont, ScaleFont;
	Datensatz *Daten;
    bool DMarkers, ShowAssignmentsOnTop, showPoints;
	Marker *marker;
    int AnzahlMarker, nDatenS, cPosx, cPosy, sPoint, mPoint;
	QFont AssFont;
	QColor SpektColor;
	QColor *CopyColor;
	QToolButton *ZoomB;
	QImage *Image;
	
protected slots:
    virtual void languageChange();
    void mouseReleased(QMouseEvent *e);
    void mouseMoved(QMouseEvent *e);
    void mousePressed(QMouseEvent *e);

    virtual inline void addPoint()
    {
        cPosx = cPosy = 0;
    }

    virtual inline void removePoint()
    {
        sPoint = -1;
    }

private slots:
	void MouseMovedOverPicture(int x, int y);
	
private:
    QPixmap image0;
    bool ScrollsEnabled, ShowMarkerLabels;
    QDoubleValidator *VEVon, *VEBis, *VWNVon, *VWNBis;

    void init();
    void destroy();

};

#endif // SPEKTRUM_H
