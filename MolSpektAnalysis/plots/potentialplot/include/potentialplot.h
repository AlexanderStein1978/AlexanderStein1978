//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//

#ifndef POTENTIALPLOT_H
#define POTENTIALPLOT_H

#include "DiagWindow.h"

class Potential;
class MainWindow;
struct SplinePoint;


class PotentialPlot : public DiagWindow
{
	Q_OBJECT

public:
    PotentialPlot(Potential *potential = 0, MainWindow *MW = 0);
    ~PotentialPlot();
	void plotPotential(Potential *potential = 0);
	void SetBoundaries();
    void PSpektrum(QPainter &P, const QRect & A, bool PrintFN ) override;
	void ResetZoom();
	void clearHistory();
	void addPotential(Potential *pot);
    void removePotential(Potential *pot);
	void setShowHistory(bool H);
	bool getShowHistory();
	void PotSnapShot();
	bool getShowDiagFuncs();
	void setShowDiagFuncs(bool show);
    void setRRange(const double RMin, const double RMax);
	
public slots:
	void PotentialChanged();

protected:
    void closeEvent(QCloseEvent *E) override;

signals:
    void closing();
	
private slots:
	void addPoint() override;
	void removePoint() override;
	
private:
    void SetPoints() override;
    void MovePoint() override;
    void MovePoint(int i_pointIndex, double i_newX, double i_newY) override;
    void ShowPopupMenu(const QPoint &i_point) override;
    void HandleHistoryWhileMoving() override;
    void HandleHistoryAfterMoving() override;

	Potential *pot, **potCopies;
    int NCopies, copyP;
    bool History, showDiagFuncs;
    int nDFP;
	QMenu *pointMenu, *pictureMenu;
	QAction *removePointAct, *addPointAct;
    double *WFS, *WWFS, DFRi, DFh;
};

#endif
