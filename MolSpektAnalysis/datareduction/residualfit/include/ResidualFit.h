//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//

#ifndef RESIDUALFIT_H
#define RESIDUALFIT_H


#include <QObject>

#include "Spline.h"

class ConnectSpline;
class ElState;

struct Point;
struct Js;


class ResidualFit : public QObject
{
    Q_OBJECT

public:
    ResidualFit(const ElState* i_state);
	virtual ~ResidualFit();
	void FitTo(int N, Point *points);
    void writeData(QTextStream* i_stream) const;
    bool readData(QTextStream* i_stream);
    void setAssignment(const QString* i_stateName, const int i_iso, const int i_v, const int i_comp);
    void RemoveLocalPerturbation(const int i_index);
    int addLocalPerturbation(const double i_IsoF, const double i_Omega, const std::vector<Js>& i_data, const int i_JStep,
                                     const double i_BeCurState);
    void setFitData(const int N, const Point *const points, const double i_IsoF, const double i_Omega, const int i_JStep, const double i_BeCurState);
    void DetermineUncertaintiesByBeVariation(LocalPerturbation& io_perturbation);
    double GetPoint(const double JValue, const int pertIndex);
    double GetPoint(const double JValue);
    double GetPointPerturber(const double JValue, const int pertIndex);
    double GetPointFromSplines(const double JValue);
    void JoinSplines(int i_indexFirstSpline, const std::vector<Js> &i_data);

    inline int getNumSplines()
	{
		return NumSplines;
	}
	
	inline Spline* getSpline(int n)
	{
		return splines[n];
	}

    inline int getIso() const
    {
        return iso;
    }

    inline int getv() const
    {
        return v;
    }

    inline int getComp() const
    {
        return comp;
    }

    inline const QString* getStateName() const
    {
        return stateName;
    }

    inline void SetState(ElState * const i_state)
    {
        m_state = i_state;
    }

    inline const LocalPerturbation* getLocalPerturbation(const int i_index) const
    {
        return m_localPerturbations[i_index];
    }

    inline LocalPerturbation* getLocalPerturbation(const int i_index)
    {
        return m_localPerturbations[i_index];
    }

    inline int getNumberOfLocalPerturbations() const
    {
        return m_localPerturbations.size();
    }

    inline bool isFitDataAvailable() const
    {
        return (NumSplines > 0 && splines[0]->isFitDataAvailable());
    }

signals:
    void Changed();

private:
    ResidualFit(const ResidualFit&);
    ResidualFit& operator=(const ResidualFit&);

	void clear();
	bool inCycle(int N, SplinePoint* spoints);
	void clearCycleDetectList();
	void setAverage(int N, SplinePoint* spoints);
    void DetermineUncertaintiesByBeVariation(LocalPerturbation& io_perturbation, const Spline& i_spline,
         const int i_nData, Point *i_backupData, Point *i_splineFitData, const std::vector<Js> &i_data,
         std::vector<Js> i_pertFitData) const;
    bool DetermineExtremalBeOneStep(double &o_H12Val, double *const o_BandCValues, bool& o_maxReached, double& o_startSigma,
         LocalPerturbation &i_localPert, double i_initStep, const double i_startBe) const;
    bool DetermineExtremalBeOneDirection(double& o_H12Val, double* o_BandCVal,
         LocalPerturbation i_perturbation, Spline i_spline, const double i_dir, const int i_nData, Point *i_backupData,
         Point *i_splineFitData, const std::vector<Js> &i_data, std::vector<Js> i_pertFitData) const;
    void DetermineMinDUpMaxDDown(double& MinDUp, double& MaxDDown, int &o_nmUp, int &o_nmDown, const LocalPerturbation& i_localPert) const;
    void ConnectSplines();
	
	bool cycleDetected;
    int NumSplines, cycleIndex, iso, v, comp;
	Spline** splines;
	QList<SplinePoint*> *cycleDetectList;
    QString *stateName;
    const ElState *m_state;
    std::vector<LocalPerturbation*> m_localPerturbations;
    std::vector<ConnectSpline*> m_connectSplines;
};

#endif
