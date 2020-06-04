//
// C++ Interface: FCFTab
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2007 - 2016
//
// Copyright: See README file that comes with this source code
//
//


#ifndef FCFTAB_H
#define FCFTAB_H


#include <tablewindow.h>

class MainWindow;
class Molecule;


class FCFTab : public TableWindow
{
	Q_OBJECT
	
public:
    FCFTab(MainWindow* MW = 0, Molecule* M = 0);
	~FCFTab();
	void setData(QString Pot1, QString Pot2, int *NIso, int **NJ, int ***Nvs, int ***Nvss, 
				 float *****Data);
    bool calcFCF(ElState *State1, ElState *State2, int Maxvs, int Maxvss, int MaxJ, int NumWFPoints);
	void setStates(ElState *State1, ElState *State2);
	void setTransition(Transition *T);
	Transition *getTransition();
	ElState *getState1();
	ElState *getState2();
    virtual bool readData(QString Filename = "");
    virtual bool writeData(QString Filename = "");
	void getFCF(int PQR, int Iso, int Js, int vs, int &Nvss, float *&Data);
	void getFCF(int JD, int Iso, int Js, int vs, int &Nvss, double *Data);
	
private slots:
	void IsoChanged(int NIso);
	void JsChanged();
	void JssChanged();

private:
	void DestroyData();
	void cleanUp(double *E, double **WF, int *ri, int *ra, int *P, double *F, double *N, int Nv);
    void cleanUp(double *E, double ***WF, int *ri, int *ra, int *P, double *F, double *N, const int NChannels, const int Nv);
	void setData(int PQR, int Iso, int J, int Nvs, int Nvss, double ***FCF);
    void setData(const int PQR, const int Iso, const int J, const int Nvs, const int Nvss, double ****FCF);
	void shiftPar(double *&E1, double *&E2, double ***&WF1, double ***&WF2, int *&ri1, int *&ri2,
			      int *&ra1, int *&ra2, int *&P1, int *&P2, double *&F1, double *&F2, double *&N1,
				  double *&N2, int &Nv1, int &Nv2);
	
    int *NIso, **NJ, ***Nvs, ***Nvss, nIso, aPQR, m_NChannelsQ, m_NChannelsPR;
	float *****Data;
    double *m_ChannelFQ, *m_ChannelFPR;
	bool Block;
	ElState *St1, *St2;
	Transition *Trans;
};

#endif // FCFTAB_H
