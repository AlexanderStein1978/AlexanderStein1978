//
// C++ Interface: Datensatz
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2006 - 2020
//
// Copyright: See README file that comes with this source code
//
//

#ifndef DATENSATZ_H
#define DATENSATZ_H


struct Marker;
struct Element;


class Datensatz      //Einlesefunktion nur für Einlesen in einem Stück gedacht!
{
public:    
    Datensatz();
    ~Datensatz();
    void reinit();
    int GetDSL() const;
    void AddValue(const double &x, const double &y, const bool marked);
    void SubtractLine(const double Xstart, const double Xend, const int Ndata, const double *const Ydiff, const bool subtract);
    void InsertML(Datensatz &ML);
    void ReverseOrder();
    void SetMarker(int i, Marker *nmarker);
    Marker *GetMarker(int i);
    double GetValue(const int &i0, const int &i1);
    bool GetMarked(int i);
	double *GetPoint(int i);
private:
    Element *AE;
    int AIndex, G;
};

#endif
