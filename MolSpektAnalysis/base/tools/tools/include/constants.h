//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//

#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <QString>


const QString MAVersion("0.2.9.0");

#ifdef Q_WS_WIN
const QString DIRSEP("\\");
#else
const QString DIRSEP("/");
#endif

const int MaxDunCoefficients = 200;
const int MaxAtoms = 150;
const int MaxMolecules = 150;
const int MaxTermTables = 150;
const int MaxDunTables = 150;
const int MaxPotentials = 150;
const int MaxFitDataSets = 150;
const int MaxSpectSimulations = 150;
const int MaxIso = 50;
const int MaxStates = 50;
const int MaxTransitions = 100;
const int MaxLineTables = 100;
const int MaxFCFTables = 15;
const int MaxSpectra = 20;
const int MaxProgressions = 20000;
const int MaxTermRows = 700000;

const int cMaxv = 200;
const int cMaxJ = 500;
const double rmin = 1.8;
const double rmax = 200.0;
const double CMaxSearchDev = 0.02;
const int  NumPoints = 500000;
const int  NumFCF_WFPoints = 100000;
const int  MaxvFCF = 500;
const bool SaveMemory = false;

const double C_hq = 1.054571628e-34;
const double C_h = 6.62606896e-34;
const double C_c = 299792458.0;
const double C_u = 1.660538782e-27;
const double C_kB = 1.380658e-23;
const double C_K_C = 273.15;
const double hartree_cm = 219474.6313705;
const double a0_Angstrom = 0.529177249;

const double def_L_gsq = 0.025;

enum PotentialType {NoPotential, SplinePotential, analyticalPotential, 
	MorseLongRangePotential, TangToenniesPotential, ModifiedTangToenniesPotential};

#endif
