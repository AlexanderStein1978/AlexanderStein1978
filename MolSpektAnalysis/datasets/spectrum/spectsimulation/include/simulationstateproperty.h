//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//

#ifndef SIMULATIONSTATEPROPERTY_H
#define SIMULATIONSTATEPROPERTY_H


class ElState;
class TermTable;
class FitData;
class FCFTab;


struct SimulationStateProperty
{
    enum LevelSource{UseOnlyTermTable, UseTermTableAndResidualFits, UseOnlyFitData, UseOnlyLevelsInFitDataButEnergiesFromTermTable,
                     UseOnlyLevelsInFitDataButEnergiesFromTermTableAndResidualFits};
    ElState* state;
    TermTable* termTable;
    FitData* fitData;
    FCFTab* fcfTab;
    int vMin, vMax;
    LevelSource levelSource;
};

#endif
