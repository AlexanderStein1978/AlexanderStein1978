//
// C++ Interface: constants
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2006 - 2017
//
// Copyright: See README file that comes with this source code
//
//

#ifndef MAVersion
#define MAVersion "0.2.8.0"

#ifdef Q_WS_WIN
#define DIRSEP "\\"
#else
#define DIRSEP "/"
#endif

#define MaxDunCoefficients 200
#define MaxAtoms 150
#define MaxMolecules 150
#define MaxTermTables 150
#define MaxDunTables 150
#define MaxPotentials 150
#define MaxFitDataSets 150
#define MaxSpectSimulations 150
#define MaxIso 50
#define MaxStates 50
#define MaxTransitions 100
#define MaxLineTables 100
#define MaxFCFTables 15
#define MaxSpectra 20
#define MaxProgressions 20000
#define MaxTermRows 700000

#define cMaxv 200
#define cMaxJ 500
#define rmin 1.8
#define rmax 200.0
#define CMaxSearchDev 0.02
#define NumPoints 500000
#define NumFCF_WFPoints 100000
#define MaxvFCF 500
#define SaveMemory false

#define C_hq 1.054571628e-34
#define C_h 6.62606896e-34
#define C_c 299792458.0
#define C_u 1.660538782e-27
#define C_kB 1.380658e-23
#define C_K_C 273.15
#define hartree_cm 219474.6313705
#define a0_Angstrom 0.529177249

#define def_L_gsq 0.025

enum PotentialType {NoPotential, SplinePotential, analyticalPotential, 
	MorseLongRangePotential, TangToenniesPotential, ModifiedTangToenniesPotential};

#endif
