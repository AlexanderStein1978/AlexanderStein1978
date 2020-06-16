//
// C++ Interface: residualplot
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2009 - 2020
//
// Copyright: See README file that comes with this source code
//
//

#ifndef RESIDUALPLOT_H
#define RESIDUALPLOT_H

//#define TestResidualFit


#include <QTableWidget>

#include "DiagWindow.h"
#include "ResidualFit.h"


class MainWindow;
class Molecule;
class TermTable;
class LineTable;
class FitData;
class PerturbationInfoTable;

struct Point;
struct Js;

class QPaintDevice;
class QPainter;


class ResidualPlot : public DiagWindow
{
	Q_OBJECT
	
	public:
    	ResidualPlot(MainWindow *MW);
		~ResidualPlot();
		void PSpektrum(QPainter &P, const QRect & A, bool PrintFN);
        void ClearMarked();
		
	public slots:
		void molBoxChanged();
		void sourceBoxChanged();
		void refSourceBoxChanged();
		void moleculesChanged();
		
	private slots:
		void stateBoxChanged(int Index);
		void isoBoxChanged(int Index);
		void vBoxChanged(int Index);
		void compBoxChanged(int Index);
		void Printall();
		void markLines();
		void setError();
		void add9toError();
		void sub9fromError();
		void deleteLevels();
		void showPMenu(QPoint *P);
		void PictureClicked(QPoint *P, bool ControlPressed);
		void marked(QRect *R, bool ControlPressed);
		void markUnc();
		void markComp();
        void CalcResidualFit();
        void AddLocalPerturbation();
        void RemoveLocalPerturbations();
        void ShowLocalPerturbationInfoTable();
        void ShowResidualFit();
        void JoinResidualFitSplines();
        
        inline void SwitchShowPoints()
        {
            setShowPoints(!showPoints);
        };
		
	private:
		void drawDataPoint(QPainter *P, int x, int y, int unc1, int unc2, int end);
		void destroyData();
		void initData();
		void mark(int xi, int xo, int yi, int yo, bool ControlPressed);
        Point* getResidualFitData() const;
        void SetResidualFitData();
        void SetPoints();
        void MovePoint(int i_pointIndex, double i_newX, double i_newY);
        void DetermineClickedSplineAndClickedPoint(const int i_overalIndex);
        void removePoint();
        void addPoint();
        std::vector<Js> GetPertFitData() const;

        inline void HandleHistoryWhileMoving()
        {
            DetermineClickedSplineAndClickedPoint(mPoint);
        }

        inline void MovePoint()
        {
            if (mPoint >= 0 && m_clickedSpline >= 0 && m_clickedPoint >= 0 && m_ResidualFit != 0)
                m_ResidualFit->getSpline(m_clickedSpline)->MovePoint(m_clickedPoint, points[mPoint].x, points[mPoint].y);
        }

		Molecule *mol;
		ElState *state;
		TermTable *refT, *termT;
		LineTable *lineT;
		FitData *FitD;
		QMenu *PMenu;
        QAction *setErrorAct, *add9toErrorAct, *sub9fErrorAct, *deleteLevelsAct, *markLinesAct, *markUncAct, *markCompAct;
        QAction *residualFitSeparatorAct, *showResidualFitAct, *createResidualFitAct, *addLocalPerturbationAct, *removeLocalPerturbationsAct, *showLocalPerturbationInfoTableAct;
        QAction *showResidualFitSplinePointsAct, *addResidualFitSplinePointAct, *deleteResidualFitSplinePointAct, *joinResidualFitSplinesAct;
		QList<int> CompList;
        double ****Data, ****Unc, ****Obs;
        int ***nData, *nv, ****JData, Nv, NC, NJ, NIso, v, C, Iso, ****Row, ****Comp, m_clickedSpline, m_clickedPoint;
        int MrefTcT, *refTcT, *refTIsoT, addedToFront;
        bool Block, *Marked, FSA, showResidualFit;
        ResidualFit *m_ResidualFit;
        PerturbationInfoTable* m_PerturbationInfoTable;

#ifdef TestResidualFit
		DiagWindow *SPlot;
#endif
};

#endif
