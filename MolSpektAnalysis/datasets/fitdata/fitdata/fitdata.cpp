//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//

#include "fitdata.h"
#include "fitdatasortfunctions.h"
#include "fitdatacoresortfunctor.h"
#include "sortforextractneworchangedwithoutsourcefunctor.h"
#include "sortbyltabandprogfunctor.h"
#include "linetable.h"
#include "elstate.h"
#include "molecule.h"
#include "Spektrum.h"
#include "termtable.h"
#include "utils.h"
#include "duntable.h"
#include "potential.h"
#include "naturalspline.h"
#include "tableline.h"
#include "ResidualFit.h"
#include "termenergy.h"
#include "heapsort.h"
#include "isotab.h"
#include "assign_vs_comptrans.h"
#include "progression.h"
#include "line.h"
#include "tlref.h"
#include "tsdialog.h"
#include "perturbation.h"
#include "espotfitinputelstateassigndialog.h" 
#include "mtable.h"
#include "basedata.h"

#include <QTableWidget>
#include <QStringList>
#include <QGridLayout>
#include <QListWidget>
#include <QComboBox>
#include <QLabel>
#include <QPushButton>
#include <QPainter>
#include <QFile>
#include <QTextStream>
#include <QCheckBox>
#include <QHeaderView>

#include <stdio.h>
#include <math.h>


FitData::FitData(ElState* nState, MainWindow* MW, Molecule* M): TableWindow(MDIChild::FitDataSet, MW, M)
{
    State = nState;
    setWindowTitle("Fit Data");
    fitDataCore = new FitDataCore();
	table->setModel(fitDataCore);
	table->horizontalHeader()->setVisible(true);
	table->verticalHeader()->setVisible(true);
	connect(table, SIGNAL(SelChanged()), this, SIGNAL(SelChanged()));
    Sources = 0;
    LineElStates = 0;
    NMarkedLevels = NSources = NSourceOffset = 0;
    FC = 0;
    lRow = -1;
    SourceOffset = 0;
    SourceOffsetNames = 0;
}

FitData::~FitData()
{
    if (Sources != 0) delete[] Sources;
    if (LineElStates != 0) delete[] LineElStates;
    if (FC != 0) delete[] FC;
    if (NSourceOffset > 0)
    {
        delete[] SourceOffset;
        delete[] SourceOffsetNames;
    }
    for (QList<ResidualFit*>::iterator it = residualFits.begin(); it != residualFits.end(); ++it) delete *it;
}

int FitData::addMarkedLevel(TermEnergy& TE, Spektrum *Source)
{
    NMarkedLevels++;
    table->blockSignals(true);
    int R = fitDataCore->addMarkedLevel(TE, Source);
    table->blockSignals(false);
    setBlockChangeSignal(true);
    Changed();
    setBlockChangeSignal(false);
    return R;
}

void FitData::AddRow()
{
    if (table == 0) return;
    table->blockSignals(true);
    int r, cr = table->currentIndex().row(), nr = fitDataCore->addRow(cr),  *FCB = new int[nr + 1];
    LineTable **SB = new LineTable*[nr+1];
    ElState **StB = new ElState*[nr+1];
    if (cr == -1) cr = 0;
    for (r=0; r < cr && r < NSources; r++)
    {
        SB[r] = Sources[r];
        FCB[r] = FC[r];
        StB[r] = LineElStates[r];
    }
    for (r = nr; r > cr; r--)
    {
        if (r <= NSources)
        {
            SB[r] = Sources[r-1];
            FCB[r] = FC[r-1];
            StB[r] = LineElStates[r-1];
        }
        else
        {
            SB[r] = 0;
            FCB[r] = -1;
            LineElStates[r] = 0;
        }
    }
    SB[cr] = 0;
    FCB[cr] = -1;
    LineElStates[cr] = 0;
    NSources = nr + 1;
    if (Sources != 0)
    {
        delete[] Sources;
        delete[] FC;
        delete[] LineElStates;
    }
    Sources = SB;
    FC = FCB;
    LineElStates = StB;
    table->blockSignals(false);
    Changed();
}

void FitData::clearMarkedLevels()
{
    if (NMarkedLevels > 0)
    {
        fitDataCore->setRowCount(fitDataCore->rowCount() - NMarkedLevels);
        NMarkedLevels = 0;
    }
}

void FitData::copyDataFromTable(const int i_numLines, int *const i_Lines, const FitData * const i_fitDataToCopyFrom)
{
    int currentRow = fitDataCore->rowCount(), NewRowCount = currentRow + i_numLines;
    table->blockSignals(true);
    if (currentRow > 0)
    {
        LineTable** newSources = new LineTable*[NewRowCount];
        int *newFC = new int[NewRowCount];
        ElState** newLineElStates = new ElState*[NewRowCount];
        for (int n=0; n < NSources; ++n)
        {
            newSources[n] = Sources[n];
            newFC[n] = FC[n];
            newLineElStates[n] = LineElStates[n];
        }
        for (int n = NSources; n < currentRow; ++n)
        {
            newSources[n] = 0;
            newFC[n] = -1;
            newLineElStates[n] = 0;
        }
        delete[] Sources;
        delete[] FC;
        delete[] LineElStates;
        Sources = newSources;
        FC = newFC;
        LineElStates = newLineElStates;
    }
    else
    {
        Sources = new LineTable*[NewRowCount];
        FC = new int[NewRowCount];
        LineElStates = new ElState*[NewRowCount];
    }
    NSources = NewRowCount;
    fitDataCore->addData(i_numLines, i_Lines, i_fitDataToCopyFrom->fitDataCore);
    for (int n=0; n < i_numLines; ++n)
    {
        if (i_Lines[n] < i_fitDataToCopyFrom->NSources)
        {
            Sources[currentRow] = i_fitDataToCopyFrom->Sources[i_Lines[n]];
            FC[currentRow] = i_fitDataToCopyFrom->FC[i_Lines[n]];
            LineElStates[currentRow++] = i_fitDataToCopyFrom->LineElStates[i_Lines[n]];
        }
        else
        {
            Sources[currentRow] = 0;
            FC[currentRow] = -1;
            LineElStates[currentRow++] = 0;
        }
    }
    table->blockSignals(false);
    bool different = (NSourceOffset != i_fitDataToCopyFrom->NSourceOffset);
    if (!different) for (int i=0; i < NSourceOffset; ++i) if (SourceOffset[i] != i_fitDataToCopyFrom->SourceOffset[i] || SourceOffsetNames[i] != i_fitDataToCopyFrom->SourceOffsetNames[i]) different = true;
    if (different)
    {
        if (NSourceOffset > 0)
        {
            delete[] SourceOffset;
            delete[] SourceOffsetNames;
            NSourceOffset = 0;
            SourceOffset = 0;
            SourceOffsetNames = 0;
        }
        if (i_fitDataToCopyFrom->NSourceOffset > 0)
        {
            NSourceOffset = i_fitDataToCopyFrom->NSourceOffset;
            SourceOffset = new double[NSourceOffset];
            SourceOffsetNames = new QString[NSourceOffset];
            for (int i=0; i < NSourceOffset; ++i)
            {
                SourceOffset[i] = i_fitDataToCopyFrom->SourceOffset[i];
                SourceOffsetNames[i] = i_fitDataToCopyFrom->SourceOffsetNames[i];
            }
        }
    }
    Changed();
    delete[] i_Lines;
}

bool FitData::containsState(ElState *State)
{
    int n;
    for (n=0; n < NSources; ++n) if (LineElStates[n] == State) return true;
    return false;
}

void FitData::DeleteRows()
{
    if (table == 0) return;
    table->blockSignals(true);
    int j, k, n, R = fitDataCore->rowCount(), *Rows = nullptr, N=0;
    table->getSelectedRows(Rows, N);
    if (NSources > 0) 
    {
        for (j=0; j < N && (Rows[j] >= NSources || Sources[Rows[j]] == nullptr); ++j) ;
        if (j < N)
        {
            QMessageBox::StandardButton B = QMessageBox::question(this, "MolSpektAnalysis", 
                "Delete the lines/progessions from the linetable, too?", QMessageBox::Yes | QMessageBox::No | QMessageBox::Cancel,
                QMessageBox::No);
            if (B == QMessageBox::Cancel) return;
            if (B == QMessageBox::Yes)
            {
                QList<LineTable*> SL;
                QList<int> CL;
                for (j=0; j < N; j++) if (Rows[j] < NSources && Sources[Rows[j]] != nullptr)
                {
                    for (n=0; (n < SL.count() ? SL[n] != Sources[Rows[j]] : false); n++) ;
                    if (n == SL.count())
                    {
                        SL.append(Sources[Rows[j]]);
                        CL.append(1);
                    }
                    else CL[n]++;
                }
                int **PN = new int*[SL.count()], **vss = new int*[SL.count()], **Jss = new int*[SL.count()];
                for (n=0; n < SL.count(); n++)
                {
                    PN[n] = new int[CL[n]];
                    if (SL[n]->getTransition()->getLowerState() == State)
                    {
                        vss[n] = new int[CL[n]];
                        Jss[n] = new int[CL[n]];
                    }
                    else vss[n] = Jss[n] = 0;
                    CL[n] = 0;
                }
                for (j=0; j<N; j++) if (Rows[j] < NSources)
                {
                    for (n=0; SL[n] != Sources[Rows[j]]; n++) ;
                    const BaseData* rowData = fitDataCore->getData(Rows[j]);
                    if (vss[n] != 0)
                    {
                        vss[n][CL[n]] = rowData->v;
                        Jss[n][CL[n]] = rowData->J;
                    }
                    PN[n][CL[n]++] = rowData->prog;
                }
                for (n=0; n < SL.count(); n++)
                {
                    SL[n]->DeleteRows(CL[n], PN[n], vss[n], Jss[n]);
                    delete[] PN[n];
                    if (vss[n] != 0)
                    {
                        delete[] vss[n];
                        delete[] Jss[n];
                    }
                }
                delete[] PN;
                delete[] vss;
                delete[] Jss;
            }
        }
    }
    bool* deleted = new bool[R];
    memset(deleted, 0, sizeof(bool) * R);
    for (n=0; n < N; ++n) deleted[Rows[n]] = true;
    for (j=0; j < N && !deleted[n]; j++) ;
    for (k=j+1; k < n; k++) if (!deleted[n])
    {
        if (k < NSources)
        {
            Sources[j] = Sources[k];
            FC[j] = FC[k];
            LineElStates[j] = LineElStates[k];
        }
        else if (j < NSources)
        {
            Sources[j] = 0;
            FC[j] = -1;
            LineElStates[j] = 0;
        }
        j++;
    }
    if (j < NSources) NSources = j;
    delete[] deleted;
    fitDataCore->deleteRows(Rows, N);
    delete[] Rows;
    table->blockSignals(false);
    Changed();
}

QList< LevelComb > FitData::getAvLevelComb(int iso, int comp, int vmin, int vmax)
{
    QList<LevelComb> R;
    LevelComb c;
    int n, m, N = fitDataCore->rowCount();
    for (n=0; n<N; n++) 
    {
        const BaseData* rowData = fitDataCore->getData(n);
        c.v = rowData->v;
        if (c.v < vmin || c.v > vmax) continue;
        c.iso = rowData->isotope;
        if (iso != -1 && iso != c.iso) continue;
        c.ef = 1 - abs(rowData->J - rowData->Js);
        if (comp < 2 && 1 - comp != c.ef) continue;
        for (m=0; m < R.size() && (R[m].iso != c.iso 
            || R[m].v != c.v || R[m].ef != c.ef); m++) ;
        if (m == R.size()) R.append(c);
    }
    return R;
}

void FitData::deleteRows(int* rows, int N)
{
    int n, m, NR = fitDataCore->rowCount();
    int mb = NR - NMarkedLevels;
    table->blockSignals(true);
    if (NSources > 0) 
    {
        for (n=0; (n<N ? (rows[n] < NSources ? Sources[rows[n]] == 0 : true) : false); n++) ;
        if (n<N ? Sources[rows[n]] != 0 : false)
        {
            QMessageBox::StandardButton B = QMessageBox::question(this, "MolSpektAnalysis", 
                "Delete the lines/progessions from the linetable, too?", QMessageBox::Yes | QMessageBox::No | QMessageBox::Cancel,
                QMessageBox::No);
            if (B == QMessageBox::Cancel) return;
            if (B == QMessageBox::Yes)
            {
                QList<LineTable*> SL;
                QList<int> CL;
                for (m=0; m<N; m++) if (rows[m] < NSources)
                {
                    for (n=0; (n < SL.count() ? SL[n] != Sources[rows[m]] : false); n++) ;
                    if (n == SL.count())
                    {
                        SL.append(Sources[rows[m]]);
                        CL.append(1);
                    }
                    else CL[n]++;
                }
                int **PN = new int*[SL.count()], **vss = new int*[SL.count()], **Jss = new int*[SL.count()];
                for (n=0; n < SL.count(); n++)
                {
                    PN[n] = new int[CL[n]];
                    if (SL[n]->getTransition()->getLowerState() == State)
                    {
                        vss[n] = new int[CL[n]];
                        Jss[n] = new int[CL[n]];
                    }
                    else vss[n] = Jss[n] = 0;
                    CL[n] = 0;
                }
                for (m=0; m<N; m++) if (rows[m] < NSources)
                {
                    for (n=0; SL[n] != Sources[rows[m]]; n++) ;
                    const BaseData* rowData = fitDataCore->getData(rows[m]);
                    if (vss[n] != 0)
                    {
                        vss[n][CL[n]] = rowData->v;
                        Jss[n][CL[n]] = rowData->J;
                    }
                    PN[n][CL[n]++] = rowData->prog;
                }
                for (n=0; n < SL.count(); n++)
                {
                    SL[n]->DeleteRows(CL[n], PN[n], vss[n], Jss[n]);
                    delete[] PN[n];
                    if (vss[n] != 0)
                    {
                        delete[] vss[n];
                        delete[] Jss[n];
                    }
                }
                delete[] PN;
                delete[] vss;
                delete[] Jss;
            }
        }
    }
    bool* deleted = new bool[NR];
    memset(deleted, 0, NR);
    for (n=0; n<N; n++)
    {
        deleted[rows[n]] = true;
        if (rows[n] >= mb) NMarkedLevels--;
    }
    for (n=0; n < NR && !deleted[n]; n++) ;
    for (m=n+1; m < NR; m++) if (deleted[n])
    {
        if (n < NSources)
        {
            if (m < NSources)
            {
                FC[n] = FC[m];
                Sources[n] = Sources[m];
                LineElStates[n] = LineElStates[m];
            }
            else NSources = n;
        }
        n++;
    }
    if (NSources > n) NSources = n;
    fitDataCore->deleteRows(rows, N);
    table->blockSignals(false);
    Changed();
}

void FitData::findWrongData()
{
    if (State == 0)
    {
        QMessageBox::information(this, "MolSpektAnalysis", 
                        "The data set has to be assigned to an electronic state first!");
        return;
    }
    int r, NR = fitDataCore->rowCount(), l = State->getLambda(), c, I, J;
    int NI = molecule->getNumIso(), NC = (l==0 ? 1 : 2);
    int JStep[NI], **JStart = CreateInt(NI, NC);
    Transition *Tr;
    ElState *St;
    std::string B;
    if (NR == 0) return;
    for (r=0; r < NI; r++) 
    {
        JStep[r] = molecule->getJStep(r);
        for (c=0; c < NC; c++) JStart[r][c] = State->getJStart(r, 1 - c);
    }
    if (lRow == -1) lRow = NR;
    const BaseData* rowData = fitDataCore->getData(r);
    int vs;
    for (r = lRow + 1; r != lRow; r++) 
    {
        if (r >= NR) if (lRow == (r=0)) break;
        rowData = fitDataCore->getData(r);
        J = rowData->J;
        I = rowData->isotope;
        c = (rowData->Js == J ? 1 : 0);
        vs = stoi(B = rowData->vs);
        if (c == 1 && l == 0 ) 
            if (B == "TE" || vs < -1 ? true : (r < NSources ?
                (Sources[r] != 0 ? ((Tr = Sources[r]->getTransition()) != 0 ? 
                  ((St = Tr->getUpperState()) != 0 ? St->getLambda() == 0 
                    : false) : false) : false) : false)) break;
        if (((J - JStart[I][c]) / JStep[I]) * JStep[I] + JStart[I][c] != J) break;
        if (fabs(rowData->devR) >= 3.5) break;
    }
    if (r != lRow) 
    {
        table->selectRow(lRow = r);
        if (r < NSources && rowData->vs != "TE" && Sources[r] != 0 && vs >= 0)
        {
            double WN = rowData->energy;
            I = rowData->isotope;
            J = rowData->Js;
            if (vs == 0 && (B == "TE" || B == "nA")) c = -1;
            Sources[r]->MarkLines(&I, &c, &J, &WN, 1);
        }
    }
    else QMessageBox::information(this, "MolSpektAnalysis", "No bad levels found!");
    Destroy(JStart, NI);
}

QList< int > FitData::getaFC()
{
    int i, j;
    QList<int> R;
    for (i=0; i < NSources; i++)
    {
        for (j=0; j < R.count() && FC[i] != R[j]; j++) ;
        if (j == R.count()) R.append(FC[i]);
    }
    for (i=0; i < R.count() - 1; i++) for (j=i+1; j < R.count(); j++) if (R[i] > R[j]) R.swapItemsAt(i, j);
    return R;
}

void FitData::getav(int &nv, int*& v)
{
    int n, b, N = getNumLines(), mv=0, Mv=0;
    for (n=0; n<N; n++)
    {
        if ((b = fitDataCore->get_v(n)) < mv) mv = -1;
        if (b > Mv) Mv = b;
    }
    bool *va = new bool[Mv - mv + 1];
    for (n=0; n <= Mv - mv; n++) va[n] = false;
    for (n=0; n<N; n++) va[((b = fitDataCore->get_v(n)) >= 0 ? b - mv : 0)] = true;
    for (n = nv = 0; n <= Mv - mv; n++) if (va[n]) nv++;
    v = new int[nv];
    for (n=mv, b=0; n <= Mv; n++) if (va[n]) v[b++] = n;
    delete[] va;
}

void FitData::getavIso(int*& Iso, int& NIso)
{
    int r, i, NR = fitDataCore->rowCount(), N = molecule->getNumIso();
    bool IsoAv[N];
    for (r=0; r<N; r++) IsoAv[r] = false;
    for (r=0; r < NR; r++) 
    {
        int iso = fitDataCore->getIso(r);
        if ((i = (iso != 0 ? iso : -1)) >= 0 ? i<N : false) IsoAv[i] = true;
    }
    for (i = NIso = 0; i<N; i++) if (IsoAv[i]) NIso++;
    if (NIso > 0) Iso = new int[NIso];
    for (i=r=0; i<N; i++) if (IsoAv[i]) Iso[r++] = i;
}

void FitData::getData(TableLine*& Lines, int& NLines, int JD, int F, int v, int mJ, int Iso, ElState *state)
{
    
    if ((NLines = getNumLines(JD, F, v, mJ, Iso)) > 0) 
    {
        Lines = new TableLine[NLines];
        int N = getNumLines();
        int *SA = new int[N], n, *S1 = heapSort(sortByProg);
        for (n=0; n < N; n++) SA[S1[n]] = n;
        delete[] S1;
        getData(Lines, SA, N, NLines, 0, v, 0, mJ, JD, F, Iso, true, state);
        delete[] SA;
    }
    else Lines = 0;
}

void FitData::getData(TableLine*& Lines, int& NLines, int *&RowN,
                      bool sortFuncs(const FitDataCore *const Tab, const int n, const int m), int *mv, int mJ)
{
    if ((NLines = (mv != 0 ? getNumLines(mv, mJ) : getNumLines())) != 0)
    {
        int n, N = getNumLines();
        int *SA = new int[N], *S1 = heapSort(sortFuncs);
        RowN = new int[NLines];
        for (n=0; n<N; n++) SA[S1[n]] = n;
        delete[] S1;
        Lines = new TableLine[NLines];
        getData(Lines, SA, N, NLines, RowN, -2, mv, mJ);
        delete[] SA;
    }
    else Lines = 0;
}

void FitData::getData(TableLine* Lines, int* SA, int i_SAL, int& NLines, int *RowN, int mv, int *Mv, int mJ, int JD, int F,
                      int Iso, bool OnlyAssignedVss, ElState* state)
{
    int n, b, m, MCT = 0, *CompTrans = 0, i, NR = fitDataCore->rowCount(), nSel, *Rows;
    ElState *XState = (molecule != 0 ? molecule->getStateP(0) : 0);
    TermTable *TT = (XState != 0 ? XState->getTermTable(false) : 0);
    Transition* Trans;
    bool *RV = new bool[NR], *RS = new bool[NR];
    getViewnRows(RV);
    table->getSelectedRows(Rows, nSel);
    for (n=0; n < NR; n++) RS[n] = false;
    for (n=0; n < nSel; n++) RS[Rows[n]] = true;
    if (TT != 0) TT->getCompT(MCT, CompTrans);
    for (n=m=0; m < NLines && n < i_SAL; n++)
    {
        if (n == NR) break;
        if (state != 0 && LineElStates != 0 && NSources > SA[n] && LineElStates[SA[n]] != state) continue;
        const BaseData* row = fitDataCore->getData(SA[n]);
        if (Iso >= 0 && row->isotope != Iso) continue;
        Lines[m].Iso = row->isotope;
        Lines[m].vss = row->v;
        if (OnlyAssignedVss && Lines[m].vss < 0) continue;
        Lines[m].Jss = row->J;
        if ((mJ > 0 ? Lines[m].Jss > mJ : false) || (Mv != 0 ? (mv >= 0 ? Lines[m].vss < mv : false) 
             || (Lines[m].Jss <= mJ ? Lines[m].vss > Mv[Lines[m].Jss] : true) 
                : (mv >= 0 ? Lines[m].vss != mv : mv == -1 && Lines[m].vss >= 0)))
            continue;
        Lines[m].Js = row->Js;
        if (JD != -1 ? fabs(Lines[m].Js - Lines[m].Jss) != JD : false) continue;
        if (F != -2 ? (SA[n] < NSources ? FC[SA[n]] != F : F != -1): false) continue;
        if (row->vs == "TE")
        {
            Lines[m].vs = -1;
            Lines[m].isTE = true;
        }
        else if ((b = stoi(row->vs)) < -1)
        {
            Lines[m].vs = b;
            Lines[m].isTE = true;
        }
        else
        {
            Lines[m].isTE = false;
            Lines[m].vs = (row->vs == "nA" ? -1 : b);
        }
        Lines[m].FC = ((SA[n] < NSources ? FC[SA[n]] >= 0 : false) ? 
            (FC[SA[n]] <= MCT && CompTrans != 0 ? CompTrans[FC[SA[n]]] : FC[SA[n]]) : 0);
        //printf("Lines[%d].FC=%d\n", m, Lines[m].FC);
        Lines[m].LTab = (SA[n] < NSources ? Sources[SA[n]] : 0);
        Lines[m].SourceName = row->source.c_str();
        Lines[m].State = (SA[n] < NSources ? LineElStates[SA[n]] : 0);
        if (Lines[m].LTab == 0 && molecule != 0 && row->source != "" && !SaveMemory)
        {
            Lines[m].LTab = molecule->getLineTable(row->source.c_str());
            if (Lines[m].LTab != 0 && (Trans = Lines[m].LTab->getTransition()) != 0) Lines[m].State = (Lines[m].isTE ? Trans->getUpperState() : Trans->getLowerState());
            if (SA[n] < NSources) Sources[SA[n]] = Lines[m].LTab;
        }
        Lines[m].PN = row->prog;
        Lines[m].File = row->file.c_str();
        Lines[m].WN = row->energy;
        if (Lines[m].WN == 0.0) continue;
        for (i=0; i < NSourceOffset; i++)
            if (Lines[m].File == SourceOffsetNames[i] || row->source.c_str() == SourceOffsetNames[i]) Lines[m].WN -= SourceOffset[i];
        Lines[m].err = row->uncert;
        Lines[m].dev = row->obs_calc;
        Lines[m].DevR = row->devR;
        if (RowN != 0) RowN[m] = SA[n];
        Lines[m].isViewn = RV[SA[n]];
        Lines[m].isSelected = RS[SA[n]];
        int sigDig = -int(floor(log10(Lines[n].err)));
        if (sigDig < 4) sigDig = 4;
        QString B = QString::number(row->uncert, 'f', sigDig);
        Lines[m].nDig = B.length() - B.indexOf('.') - 1;
        Lines[m].WeightFact = 1.0;
        Lines[m++].Row = SA[n];
    }
    delete[] RV;
    delete[] RS;
    delete[] CompTrans;
    delete[] Rows;
    NLines = m;
}

bool FitData::checkSourceConnections()
{
    int n, N = fitDataCore->rowCount(), CNSources = NSources;
    QMap<QString, LineTable*>& Translations = MW->getLineTableTranslations();
    bool result = true;
    if (Sources == 0)
    {
        Sources = new LineTable*[NSources = N];
        FC = new int[N];
        LineElStates = new ElState*[N];
    }
    else if (CNSources < N)
    {
        LineTable**SB = new LineTable*[NSources = N];
        int *FCB = new int[N];
        ElState **StB = new ElState*[N];
        for (n=0; n < CNSources; ++n)
        {
            SB[n] = Sources[n];
            FCB[n] = FC[n];
            StB[n] = LineElStates[n];
        }
        while (n < NSources)
        {
            SB[n] = 0;
            FCB[n] = -1;
            LineElStates[n++] = 0;
        }
        delete[] Sources;
        delete[] FC;
        delete[] LineElStates;
        Sources = SB;
        FC = FCB;
        LineElStates = StB;
    }
    for (n=0; n<N; ++n) if (Sources[n] == 0)
    {
        const QString Name = fitDataCore->getSource(n).c_str();
        const std::string& file = fitDataCore->getSourceFile(n);
        if (file.find(Name.toStdString()) == file.npos)
        {
            if (Translations.contains(Name)) Sources[n] = Translations[Name];
            else
            {
                LineTable* LTab = molecule->getLineTable(Name);
                if (LTab == 0) LTab = MW->getLineTable(Name, molecule);
                if (LTab != 0)
                {
                    Sources[n] = LTab;
                    Translations[Name] = LTab;
                    LineElStates[n] = LTab->getElState();
                }
                else result = false;
            }
            if (Sources[n] != 0) fitDataCore->setSource(n, Sources[n]->getName().toStdString());
        }
    }
    return result;
}

void FitData::extractChangedData(const FitData * const i_fitDataToCompare, FitData * const i_FitDataToAddChangedData) const
{
    if (0 == i_FitDataToAddChangedData) return;
    bool thisWithSources = AreSourcesAvailable(), toCompareWithSources = i_fitDataToCompare->AreSourcesAvailable(), withSources = (thisWithSources && toCompareWithSources);
    bool subtractSourceOffsets = (!withSources && (thisWithSources || toCompareWithSources) &&
                                  QMessageBox::question(0, "Extract changed levels",
                                                        "Only one of the fit data sets contains source offset information. Shall this be subtracted from the table?",
                                                        QMessageBox::Yes | QMessageBox::No, QMessageBox::Yes) == QMessageBox::Yes);
    int NRnew = fitDataCore->rowCount(), NRold = i_fitDataToCompare->fitDataCore->rowCount(), *Onew = new int[NRnew], *Oold = new int[NRold], rowOld, rowNew, *LinesToExtract = new int[NRnew];
    int NumLinesToExtract = 0;
    prepareForExtractNewOrChanged(i_fitDataToCompare, NRnew, NRold, withSources, subtractSourceOffsets, Onew, Oold);
    for (rowNew = rowOld = 0; rowNew < NRnew; ++rowNew)
    {
        sortForExtractNewOrChangedOrder order;
        while ((order = getForExtractNewOrChangedOrder(i_fitDataToCompare, Onew[rowNew], Oold[rowOld], withSources, subtractSourceOffsets)) <= SFENOCenergyIsSmaller) ++rowOld;
        if (order == SFENOCenergyIsBigger || (order != SFENOCisEqual && rowOld > 0 &&
                                              getForExtractNewOrChangedOrder(i_fitDataToCompare, Onew[rowNew], Oold[rowOld - 1], withSources, subtractSourceOffsets) == SFENOCenergyIsSmaller))
            LinesToExtract[NumLinesToExtract++] = Onew[rowNew];
    }
    delete[] Onew;
    delete[] Oold;
    i_FitDataToAddChangedData->copyDataFromTable(NumLinesToExtract, LinesToExtract, this);
    i_FitDataToAddChangedData->show();
}

void FitData::extractNewData(const FitData * const i_fitDataToCompare, FitData * const i_FitDataToAddNewData) const
{
    if (0 == i_FitDataToAddNewData) return;
    bool withSources = (AreSourcesAvailable() && i_fitDataToCompare->AreSourcesAvailable());
    int NRnew = fitDataCore->rowCount(), NRold = i_fitDataToCompare->fitDataCore->rowCount(), *Onew = new int[NRnew], *Oold = new int[NRold], rowOld, rowNew, *LinesToExtract = new int[NRnew];
    int NumLinesToExtract = 0;
    prepareForExtractNewOrChanged(i_fitDataToCompare, NRnew, NRold, withSources, false, Onew, Oold);
    for (rowNew = rowOld = 0; rowNew < NRnew; ++rowNew)
    {
        sortForExtractNewOrChangedOrder order;
        while ((order = getForExtractNewOrChangedOrder(i_fitDataToCompare, Onew[rowNew], Oold[rowOld], withSources, false)) == SFENOCisSmaller) ++rowOld;
        if (order == SFENOCisBigger) LinesToExtract[NumLinesToExtract++] = Onew[rowNew];
    }
    delete[] Onew;
    delete[] Oold;
    i_FitDataToAddNewData->copyDataFromTable(NumLinesToExtract, LinesToExtract, this);
    i_FitDataToAddNewData->show();
}

ElState* FitData::getElState()
{
    return State;
}

FitData::sortForExtractNewOrChangedOrder FitData::getForExtractNewOrChangedOrder(const FitData * const i_fitDataOld, const int i_RowNew, const int i_RowOld,
                                                                                 const bool i_withSources, const bool i_subtractSourceOffset) const
{
    const std::string Sold(i_fitDataOld->fitDataCore->getOtherState(i_RowOld)), Snew(fitDataCore->getOtherState(i_RowNew));
    if (Sold < Snew) return SFENOCisSmaller;
    if (Sold > Snew) return SFENOCisBigger;
    const int IsoOld = i_fitDataOld->fitDataCore->getIso(i_RowOld), IsoNew = fitDataCore->getIso(i_RowNew);
    if (IsoOld < IsoNew) return SFENOCisSmaller;
    if (IsoOld > IsoNew) return SFENOCisBigger;
    const int vOld = i_fitDataOld->fitDataCore->get_v(i_RowOld), vNew = fitDataCore->get_v(i_RowNew);
    if (vOld < vNew) return SFENOCisSmaller;
    if (vOld > vNew) return SFENOCisBigger;
    const int Jold = i_fitDataOld->fitDataCore->getJ(i_RowOld), Jnew = fitDataCore->getJ(i_RowNew);
    if (Jold < Jnew) return SFENOCisSmaller;
    if (Jold > Jnew) return SFENOCisBigger;
    const bool efOld = (i_fitDataOld->fitDataCore->getJs(i_RowOld) == Jold), efNew = (fitDataCore->getJs(i_RowNew) == Jnew);
    if (efOld && !efNew) return SFENOCisSmaller;
    if (!efOld && efNew) return SFENOCisBigger;
    double DeltaEold = 0.0, DeltaENew = 0.0;
    if (i_withSources)
    {
        const std::string& SourceOld(i_fitDataOld->fitDataCore->getSource(i_RowOld)), SourceNew(fitDataCore->getSource(i_RowNew));
        if (SourceOld < SourceNew) return SFENOCisSmaller;
        if (SourceOld > SourceNew) return SFENOCisBigger;
    }
    else
    {
        const bool absOld = QString::number(i_fitDataOld->fitDataCore->getUncertainty(i_RowOld), 'f', 4).contains(QRegExp("[234]01"));
        const bool absNew = QString::number(fitDataCore->getUncertainty(i_RowNew), 'f', 4).contains(QRegExp("[234]01"));
        if (absOld && !absNew) return SFENOCisSmaller;
        if (!absOld && absNew) return SFENOCisBigger;
        if (i_subtractSourceOffset)
        {
            if (i_fitDataOld->NSourceOffset > 0)
            {
                const std::string& SourceOld(i_fitDataOld->fitDataCore->getSource(i_RowOld));
                for (int i=0; i < i_fitDataOld->NSourceOffset; ++i) if (i_fitDataOld->SourceOffsetNames[i] == SourceOld.c_str()) DeltaEold = i_fitDataOld->SourceOffset[i];
            }
            if (NSourceOffset > 0)
            {
                const std::string& SourceNew(fitDataCore->getSource(i_RowNew));
                for (int i=0; i < NSourceOffset; ++i) if (SourceOffsetNames[i] == SourceNew.c_str()) DeltaENew = SourceOffset[i];
            }
        }
    }
    const double Eold = i_fitDataOld->fitDataCore->getEnergy(i_RowOld) - DeltaEold, Enew = fitDataCore->getEnergy(i_RowNew) - DeltaENew;
    if (Eold < Enew) return SFENOCenergyIsSmaller;
    if (Eold > Enew) return SFENOCenergyIsBigger;
    return SFENOCisEqual;
}

bool FitData::AreSourcesAvailable() const
{
    int N=0, NR = fitDataCore->rowCount();
    for (int n=0; n < NR; ++n) if (fitDataCore->getSource(n) != "") ++N;
    return (N > NR / 2);
}


int FitData::getMaxJ()
{
    return fitDataCore->getMaxJ();
}

int FitData::getMaxv()
{
    return fitDataCore->getMaxv();
}

bool FitData::isDataAvailable()
{
    if (fitDataCore->rowCount() > 0) return true;
    return false;
}

void FitData::prepareForExtractNewOrChanged(const FitData * const i_fitDataOld, const int i_NRnew, const int i_NRold, const bool i_withSources,
                                            const bool i_subtractSourceOffsets, int * const io_Onew, int * const io_Oold) const
{
    int *S0new, *S0old;
    if (i_withSources)
    {
        S0new = heapSort(sortForExtractNewOrChanged);
        S0old = i_fitDataOld->heapSort(sortForExtractNewOrChanged);
    }
    else
    {
        if (i_subtractSourceOffsets)
        {
            S0new = utils::heapSort(SortForExtractNewOrChangedWithoutSourceFunctor(fitDataCore, SourceOffsetNames, SourceOffset, NSourceOffset), fitDataCore->rowCount());
            S0old = utils::heapSort(SortForExtractNewOrChangedWithoutSourceFunctor(i_fitDataOld->fitDataCore, i_fitDataOld->SourceOffsetNames, i_fitDataOld->SourceOffset,
                                                                                   i_fitDataOld->NSourceOffset), i_fitDataOld->fitDataCore->rowCount());
        }
        else
        {
            S0new = utils::heapSort(SortForExtractNewOrChangedWithoutSourceFunctor(fitDataCore, 0, 0, 0), fitDataCore->rowCount());
            S0old = utils::heapSort(SortForExtractNewOrChangedWithoutSourceFunctor(i_fitDataOld->fitDataCore, 0, 0, 0), i_fitDataOld->fitDataCore->rowCount());
        }
    }
    for (int n=0; n < i_NRnew; ++n) io_Onew[S0new[n]] = n;
    for (int n=0; n < i_NRold; ++n) io_Oold[S0old[n]] = n;
    delete[] S0new;
    delete[] S0old;
}

bool FitData::readData(QTextStream& S)
{
    QString startSpecialPart = fitDataCore->readData(S);
    if (!startSpecialPart.isEmpty()) return ReadSpecialPart(S, startSpecialPart);
    return true;
}

bool FitData::readData(QString Filename)
{
    if (NSourceOffset > 0)
    {
        delete[] SourceOffset;
        delete[] SourceOffsetNames;
        NSourceOffset = 0;
        SourceOffsetNames = 0;
        SourceOffset = 0;
    }
    for (QList<ResidualFit*>::const_iterator it = residualFits.begin(); it != residualFits.end(); ++it)
    {
        disconnect(*it, SIGNAL(Changed()), this, SLOT(Changed()));
        delete *it;
    }
    residualFits.clear();

    if (!TableWindow::readData(Filename)) 
    {
        QFile F(Filename);
        F.open(QIODevice::ReadOnly);
        QTextStream S(&F);
        int r, c, NR = S.readLine().left(5).toInt(), mv = 0, mJ=0;
        bool Error = false;
        if (NR == 0)
        {
            QString B = S.readLine();
            NR = B.left(5).toInt();
            int n = B.mid(15, 5).toInt(), m = B.mid(20, 5).toInt();
            if (NR > 0 && n > 0 && m > 0)
            {
                if (MW == 0) return false;
                for (r=0; r <= n + m; ++r) S.readLine();
                QDialog *D = new QDialog;
                D->setWindowTitle("Import fit data");
                QGridLayout *L = new QGridLayout(D);
                L->addWidget(new QLabel("Please select the molecule the fit data belongs to.", D), 0, 0, 1, 2);
                L->addWidget(new QLabel("Selected molecule:", D), 1, 0);
                QComboBox *Box = new QComboBox(D);
                for (n=0; n < MW->getNumMolecules(); ++n) Box->addItem(MW->getMolecule(n)->getName());
                Box->setEditable(false);
                L->addWidget(Box, 1, 1);
                L->setRowMinimumHeight(2, 20);
                QPushButton *OK = new QPushButton("OK", D), *Cancel = new QPushButton("Cancel", D);
                L->addWidget(OK, 3, 0);
                L->addWidget(Cancel, 3, 1);
                connect(OK, SIGNAL(clicked()), D, SLOT(accept()));
                connect(Cancel, SIGNAL(clicked()), D, SLOT(reject()));
                if (D->exec() == QDialog::Accepted)
                {
                    molecule = MW->getMolecule(Box->currentIndex());
                    delete D;
                }
                else
                {
                    delete D;
                    return false;
                }
                IsoTab *Iso = molecule->getIso();
                fitDataCore->setRowCount(NR);
                int vs;
                QList<int> stateList;
                if (Sources != 0) delete[] Sources;
                Sources = new LineTable*[NSources = NR];
                if (FC == 0) delete[] FC;
                FC = new int[NR];
                if (LineElStates != 0) delete[] LineElStates;
                LineElStates = new ElState*[NR];
                table->blockSignals(true);
                for (r=0; r < NR && !S.atEnd(); ++r)
                {
                    BaseData data;
                    QStringList Buffer = S.readLine().split(' ', Qt::SkipEmptyParts);
                    if (Buffer.size() < 12) break;

                    data.isotope = Iso->getIsoIndex(Buffer[0].toInt(), Buffer[1].toInt());
                    data.v = Buffer[2].toInt();
                    if (data.v > mv) mv = data.v;
                    data.J = Buffer[3].toInt();
                    if (data.J > mJ) mJ = data.J;
                    data.vs = Buffer[4].toStdString();
                    data.Js = Buffer[5].toInt();
                    if (S.atEnd()) break;
                    vs = Buffer[4].toInt();
                    Buffer = S.readLine().split(' ', Qt::SkipEmptyParts);
                    if (Buffer.size() < 3) break;
                    data.energy = Buffer[0].toDouble();
                    data.uncert = Buffer[1].toDouble();
                    if (vs < 0 && !stateList.contains(vs)) stateList.push_back(vs);
                    FC[r] = -1;
                    Sources[r] = 0;
                    LineElStates[r] = 0;
                    fitDataCore->setRow(data, r);
                }
                if (r < NR) Error = true;
                else if (stateList.size() > 1)
                {
                    int N = stateList.size(), NSt = molecule->getNumStates();
                    QDialog *D = new QDialog(this);
                    QGridLayout *L = new QGridLayout(D);
                    QComboBox **BoxArray = new QComboBox*[N];
                    QPushButton *OK = new QPushButton("OK", D), *Cancel = new QPushButton("Cancel", D);
                    D->setWindowTitle("Different electronic states contained");
                    L->addWidget(new QLabel("Please select state assignment:", D), 0, 0, 1, 2);
                    L->addWidget(new QLabel("number:", D), 1, 0);
                    L->addWidget(new QLabel("assignment:", D), 1, 1);
                    for (int n=0; n<N; ++n)
                    {
                        L->addWidget(new QLabel(QString::number(stateList[n]), D), 2+n, 0);
                        L->addWidget(BoxArray[n] = new QComboBox(D), 2+n, 1);
                        BoxArray[n]->setEditable(false);
                        for (int s=0; s < NSt; ++s) BoxArray[n]->addItem(molecule->getStateP(s)->getName());
                    }
                    L->setRowMinimumHeight(2 + NSt, 20);
                    L->addWidget(OK, 3 + NSt, 0);
                    L->addWidget(Cancel, 3 + NSt, 1);
                    connect(OK, SIGNAL(clicked()), D, SLOT(accept()));
                    connect(Cancel, SIGNAL(clicked()), D, SLOT(reject()));
                    if (D->exec() == QDialog::Accepted)
                    {
                        ElState** StateArray = new ElState*[N];
                        for (int n=0; n<N; ++n) StateArray[n] = molecule->getStateP(BoxArray[n]->currentIndex());
                        for (r=0; r < NR; ++r) fitDataCore->setSecondState(r, (LineElStates[r] = StateArray[stateList.indexOf(stoi(fitDataCore->get_vs(r)))])->getName().toStdString());
                        delete[] StateArray;
                    }
                    delete D;
                    delete[] BoxArray;
                }
                table->blockSignals(false);
                molecule = 0;
            }
            else return false;
        }
        else
        {
            QStringList Buffer = S.readLine().split(' ', Qt::SkipEmptyParts);
            if (Buffer.count() >= 2)
            {
                setName("FitData" + ((r = Buffer[1].lastIndexOf('.')) != -1 ? Buffer[1].left(r) : Buffer[1]));
                setSource("Asen");
                setFileName(((r = Filename.lastIndexOf('.')) != -1 ? Filename.left(r) : Filename) + ".fdat");
            }
            else Error = true;
            if (Sources != 0) delete[] Sources;
            Sources = new LineTable*[NSources = NR];
            if (FC == 0) delete[] FC;
            FC = new int[NR];
            if (LineElStates != 0) delete[] LineElStates;
            LineElStates = new ElState*[NR];
            for (r=0; r<4; r++) S.readLine();
            table->blockSignals(true);
            fitDataCore->setRowCount(NR);
            for (r=0; r < NR; r++)
            {
                if (S.atEnd())
                {
                    Error = true;
                    break;
                }
                Buffer = S.readLine().split(' ', Qt::SkipEmptyParts);
                if (Buffer.count() < 7)
                {
                    Error = true;
                    continue;
                }
                BaseData data;
                data.isotope = (Buffer[6].toInt() - 1) / 10;
                data.v = Buffer[1].toInt();
                if (data.v > mv) mv = data.v;
                data.J = Buffer[2].toInt();
                if (data.J > mJ) mJ = data.J;
                data.vs = "TE";
                data.Js = (Buffer[4] == "f" ? data.J : data.J + 1);
                data.energy = Buffer[3].toDouble();
                data.uncert = Buffer[5].toDouble();
                Sources[r] = 0;
                LineElStates[r] = 0;
                FC[r] = ((c = Buffer[4].toInt()) > 0 ? c : -1);
            }
            if (Error) fitDataCore->setRowCount(r);
            table->blockSignals(false);
        }
        setvMax(mv);
        setJMax(mJ);
        emit propertiesChanged();
        return !Error;
    }
    table->blockSignals(true);
    int n, N = fitDataCore->rowCount(), b;
    if (Sources != 0) delete[] Sources;
    Sources = new LineTable*[NSources = N];
    if (FC == 0) delete[] FC;
    FC = new int[N];
    if (LineElStates != 0) delete[] LineElStates;
    LineElStates = new ElState*[N];
    for (n=0; n<N; n++)
    {
        if (n>0 && fitDataCore->getSource(n) == fitDataCore->getSource(n-1)) Sources[n] = Sources[n-1];
        else Sources[n] = (molecule != 0 && !SaveMemory ? molecule->getLineTable(fitDataCore->getSource(n).c_str()) : 0);
        if (Sources[n] != 0) LineElStates[n] = Sources[n]->getElState();
        if ((b = stoi(fitDataCore->get_vs(n))) == -1) fitDataCore->set_vs(n, "TE");
        if (b >= 990)
        {
            FC[n] = (b + 10) / 1000 - 1;
            fitDataCore->set_vs(n, (b -= 1000 * (FC[n] + 1)) == -1 ? "TE" : (b == -10 ? "nA" : std::to_string(b)));
        }
        else FC[n] = -1;
        if (molecule != 0)
        {
            QString Path = fitDataCore->getSourceFile(n).c_str(), MolPath = molecule->getFileName();
            fitDataCore->setSourceFile(n, getAbsolutePath(Path, MolPath).toStdString());
            if (n>0 && fitDataCore->getOtherState(n) == fitDataCore->getOtherState(n-1)) LineElStates[n] = LineElStates[n-1];
            else LineElStates[n] = molecule->getState(fitDataCore->getOtherState(n).c_str());
        }
    }
    table->blockSignals(false);
    emit propertiesChanged();
    return true;
}

bool FitData::ReadSpecialPart(QTextStream &i_stream, const QString &i_startString)
{
    bool Success = true;
    QString Buffer;
    QStringList SourceOffsetList;
    if (i_startString.indexOf("SourceOffsets:") >= 0)
    {
        while (!i_stream.atEnd())
        {
            Buffer = i_stream.readLine();
            if (Buffer.indexOf("Begin ResidualFit", 0, Qt::CaseInsensitive) >= 0) break;
            else SourceOffsetList << Buffer;
        }
        if (SourceOffsetList.size() > 0)
        {
            SourceOffset = new double[SourceOffsetList.size()];
            SourceOffsetNames = new QString[SourceOffsetList.size()];
            for (int n=0; n < SourceOffsetList.size(); ++n)
            {
                QStringList currentOffset = SourceOffsetList[n].split('\t');
                if (currentOffset.size() < 2) continue;
                SourceOffsetNames[NSourceOffset] = currentOffset[0];
                SourceOffset[NSourceOffset] = currentOffset[1].toDouble();
                if (!SourceOffsetNames[NSourceOffset].isEmpty() && 0.0 != SourceOffset[NSourceOffset]) ++NSourceOffset;
            }
        }
    }
    if (i_startString.indexOf("Begin ResidualFit", 0, Qt::CaseInsensitive) >= 0 || Buffer.indexOf("Begin ResidualFit", 0, Qt::CaseInsensitive) >= 0)
    {
        do
        {
            ResidualFit* newFit = new ResidualFit(0);
            Success = newFit->readData(&i_stream);
            if (!Success) break;
            if (molecule != 0) newFit->SetState(molecule->getState(*newFit->getStateName()));
            connect(newFit, SIGNAL(Changed()), this, SLOT(Changed()));
            residualFits.push_back(newFit);
            if (i_stream.atEnd()) break;
            Buffer = i_stream.readLine();
        } while(Buffer.indexOf("Begin ResidualFit", 0, Qt::CaseInsensitive) >= 0);
    }
    return Success;
}

void FitData::Assign_v(double**** TE, int NC, int NI, int NJ, int *Nv, int* IsoTrans, Assign_vs_CompTrans* CompTrans, double Tolerance,
                       double DoublAssTol)
{
    int ME, NE = 0, I, c, Jss, Js, n, m, i, j=0, NR = fitDataCore->rowCount();
    int MC = 0;
    for (n=0; n < NC; n++) if (CompTrans[n].FC > MC) MC = CompTrans[n].FC;
    ElState *StateTrans[MC + 1];
    for (n = m = ME = 0; n <= MC; n++)
    {
        if (n > CompTrans[m].FC) ++m;
        StateTrans[n] = CompTrans[m].State;
        ME += Nv[n];
    }
    double E[ME], LE, LLE = 0.0, Dev;
    int va[ME], CA[ME], *S1 = heapSort(sortIefJFreq), *SA = new int[NR], lv = -1, lc = -1, fc, k;
    for (n=0; n < NR; n++) SA[S1[n]] = n;
    int vc[MC+1];
    QString B;
    table->blockSignals(true);
    for (m=n=0, I = c = Jss = Js = -1; m <= NR; m++)
    {
        if (m < NR && I == fitDataCore->getIso(SA[m]) && Jss == fitDataCore->getJ(SA[m]) && Js == fitDataCore->getJs(SA[m])) continue;
        if (NE > 0) for (j=n, i=0; j<m; j++)
        {
            for (LE = fitDataCore->getEnergy(SA[j]); i < NE && E[i] < LE; i++) ;
            if (i>0 ? fabs(E[i] - LE) > fabs(E[i-1] - LE) : false) i--;
            if (fabs(Dev = LE - E[i]) <= Tolerance)
            {
                if (lv == va[i] && lc == CA[i])
                {
                    if (LE - LLE < DoublAssTol && fitDataCore->getSourceFile(SA[j]) != fitDataCore->getSourceFile(SA[j-1])) fc = lc;
                    else
                    {
                        if (fabs(E[i] - LE) < fabs(E[i] - LLE))
                        {
                            if (j-2 >= n ? lv == fitDataCore->get_v(SA[j-2]) && (k = stoi(fitDataCore->get_vs(SA[j-2]))) < -1 - MC : false)
                            {
                                fitDataCore->set_vs(SA[j-1], std::to_string(k-1));
                                FC[SA[j-1]] = -k;
                            }
                            else 
                            {
                                fitDataCore->set_vs(SA[j-1], std::to_string(-2 - MC));
                                FC[SA[j-1]] = 1 + MC;
                            }
                            fc = lc;
                        }
                        else
                        {
                            if ((k = stoi(fitDataCore->get_vs(SA[j-1]))) < -1 - MC) fc = -k;
                            else
                            {
                                if (j-2 >= n && lv == fitDataCore->get_v(SA[j-2]) && (k = stoi(fitDataCore->get_vs(SA[j-2]))) < -1 - MC) fc = -k;
                                else fc = MC + 1;
                            }
                        }
                    }
                }
                else fc = CA[i];
            }
            else
            {
                if (j>0 && lv == va[i] && (k = stoi(fitDataCore->get_vs(SA[j-1]))) < -1 - MC) fc = -k;
                else fc = MC + 1;
                va[i] = -1;
            }
            fitDataCore->set_v(SA[j], (lv = va[i]));
            //printf("SA[%d]=%d\n", j, SA[j]);
            if (fc == 0)
            {
                fitDataCore->set_vs(SA[j], "TE");
                FC[SA[j]] = -1;
            }
            else fitDataCore->set_vs(SA[j], std::to_string(-1 - (FC[SA[j]] = fc)));
            fitDataCore->setObsCalc(SA[j], Dev);
            fitDataCore->setDevRatio(SA[j], Dev / fitDataCore->getUncertainty(SA[j]));
            LineElStates[SA[j]] = StateTrans[CA[i]];
            if (LineElStates[SA[j]] != 0) fitDataCore->setSecondState(SA[j], LineElStates[SA[j]]->getName().toStdString());
            lc = CA[i];
            LLE = LE;
        }
        if (m < NR)
        {
            I = fitDataCore->getIso(SA[n=m]);
            Jss = fitDataCore->getJ(SA[n]);
            if ((I < NI && IsoTrans[I] >= 0) && Jss < NJ)
            {
                Js = fitDataCore->getJs(SA[n]);
                c = (MC > 0 ? 1 - fabs(Jss - Js) : 0);
                lv = lc = -1;
                for (i=0; i <= MC; i++) vc[i] = 0;
                for (NE = 0; true; )
                {
                    for (i=0, E[NE] = 1e10; i < NC; ++i)
                        if (CompTrans[i].ef == c && vc[k = CompTrans[i].FC] < Nv[i] ? (LE = TE[i][IsoTrans[I]][vc[k]][Jss]) < E[NE] : false)
                            if (LE != 0.0)
                    {
                        E[NE] = LE;
                        j=k;
                    }
                    if (E[NE] != 1e10)
                    {
                        va[NE] = vc[j]++;
                        CA[NE++] = j;
                    }
                    else break;
                }
            }
            else NE = 0;
        }
    }
    table->blockSignals(false);
    delete[] S1;
    delete[] SA;
    Changed();
}

void FitData::getNumLevels(int**& LNum, int& NumIso, int& NumComp, int type, int maxv, bool ef, bool disComp, bool disElSt, QList<ElState*>* statesList)
{
    int c, NR = fitDataCore->rowCount(), n, NC = 0, st;
    QList<ElState*> StatesList;
    if (type <= 2)
    {
        bool ****L;
        int I, J, v, Nv = 0, NJ = 0;
        int **b = CreateInt(NR, 4);
        for (n = 0, NumIso = NumComp = 1; n < NR; n++)
        {
            if ((b[n][0] = fitDataCore->getIso(n)) >= NumIso) NumIso = b[n][0] + 1;
            if (disComp && (b[n][3] = (FC[n] >= 0 ? FC[n] : -1 - stoi(fitDataCore->get_vs(n)))) >= NumComp) NumComp = b[n][3] + 1;
            if (!disComp || b[n][3] < 0) b[n][3] = 0;
            if ((b[n][1] = fitDataCore->get_v(n)) >= Nv) Nv = b[n][1] + 1;
            if ((b[n][2] = fitDataCore->getJ(n)) >= NJ) NJ = b[n][2] + 1;
            if (disElSt)
            {
                QList<ElState*>::const_iterator it;
                for (it = StatesList.begin(); it != StatesList.end() && *it != LineElStates[n]; ++it) ;
                if (it == StatesList.end()) StatesList.push_back(LineElStates[n]);
            }
        }
        if (maxv >= 0 && maxv < Nv) Nv = maxv + 1;
        if (ef) NumComp *= 2;
        if (disElSt)
        {
            NC = NumComp;
            NumComp *= StatesList.count();
        }
        L = CreateBool(NumIso, NumComp, Nv, NJ);
        for (n=0; n < NumComp; n++) for (I=0; I < NumIso; I++) for (v=0; v < Nv; v++)
            for (J=0; J < NJ; J++) L[I][n][v][J] = false;
        for (n=0; n < NR; n++)
            if ((type == 0 || (type == 1 && Sources[n] != 0) || (type == 2 && Sources[n] == 0)) && b[n][1] < Nv && b[n][1] >= 0)
        {
            for (st = 0; st < StatesList.size() && StatesList[st] != LineElStates[n]; ++st) ;
            c = st * NC + (ef ? (b[n][2] == fitDataCore->getJs(n) ? 2 * b[n][3] + 1 : 2 * b[n][3]) : b[n][3]);
            L[b[n][0]][c][b[n][1]][b[n][2]] = true;
        }
        Destroy(b, NR);
        LNum = CreateInt(NumIso, NumComp);
        for (I=0; I < NumIso; I++) for (c=0; c < NumComp; c++)
            for (LNum[I][c] = v = 0; v < Nv; v++) for (J=0; J < NJ; J++)
                if (L[I][c][v][J]) LNum[I][c]++;
        Destroy(L, NumIso, NumComp, Nv);
    }
    else
    {
        SortByLTabAndProgFunctor SortOperator(fitDataCore, Sources, NSources, MW, molecule);
        int *fsa = utils::heapSort(SortOperator, NR), *Pos = new int[NR], m = 0, NP, PNum, F=0, v;
        for (n=0; n < NR; ++n) Pos[fsa[n]]=n;
        delete[] fsa;
        fsa = 0;
        for (n=0, NumIso = NumComp = 1; n < NR; ++n)
        {
            if ((maxv >= 0 && fitDataCore->get_v(n) > maxv) || (n < NSources && Sources[n] == 0)) continue;
            if ((m = fitDataCore->getIso(n)) >= NumIso) NumIso = m+1;
            if (disComp && (m = (FC != 0 && n < NSources && FC[n] >= 0 ? FC[n] : -1 - stoi(fitDataCore->get_vs(n)))) >= NumComp) NumComp = m+1;
            if (disElSt)
            {
                QList<ElState*>::const_iterator it;
                for (it = StatesList.begin(); it != StatesList.end() && *it != LineElStates[n]; ++it) ;
                if (it == StatesList.end()) StatesList.push_back(LineElStates[n]);
            }
        }
        if (ef) NumComp *= 2;
        if (disElSt)
        {
            NC = NumComp;
            NumComp *= StatesList.count();
        }
        LNum = CreateInt(NumIso, NumComp);
        for (n=0; n < NumIso; ++n) for (m=0; m < NumComp; ++m) LNum[n][m] = 0;
        LineTable* currentLTab = 0;
        Progression* Progs;
        for (n=0; n < NR; ++n)
            if ((maxv < 0 || ((v = fitDataCore->get_v(Pos[n])) <= maxv && v >= 0)) && ((Pos[n] < NSources && Sources[Pos[n]] != 0)
                || fitDataCore->getProgression(Pos[n]) > 0))
        {
            if (Sources != 0 && Pos[n] < NSources && Sources[Pos[n]] != 0 ? currentLTab != Sources[Pos[n]] : currentLTab == 0
                || currentLTab->getName().toStdString() != fitDataCore->getSource(Pos[n]))
            {
                currentLTab = (Sources != 0 && Pos[n] < NSources && Sources[Pos[n]] != 0 ? Sources[Pos[n]] : MW->getLineTable(fitDataCore->getSource(Pos[n]).c_str(), molecule));
                if (currentLTab == 0)
                {
                    printf("GetNumLevels: Pos[%d]=%d, currentLTab==0!\n", n, Pos[n]);
                    break;
                }
                currentLTab->getProgressions(NP, Progs, sortByProgression);
                m=0;
            }
            PNum = fitDataCore->getProgression(Pos[n]);
            while(m < NP && Progs[m].PNum < PNum) ++m;
            if (PNum != Progs[m].PNum)
            {
                printf("GetNumLevels: PNum=%d, Progs[%d].PNum=%d\n", PNum, m, Progs[m].PNum);
                break;
            }
            if (disComp)
            {
                F = (FC != 0 && Pos[n] < NSources && FC[Pos[n]] >= 0 ? FC[Pos[n]] : -1 - stoi(fitDataCore->get_vs(Pos[n])));
                if (F<0) F=0;
            }
            for (st = 0; st < StatesList.size() && StatesList[st] != LineElStates[Pos[n]]; ++st) ;
            c = st * NC + (ef ? (Progs[m].Js == Progs[m].L[0].Jss ? 2*F+1 : 2*F) : F);
            if (c<0 || c >= NumComp || Progs[m].Iso < 0 || Progs[m].Iso >= NumIso)
            {
                printf("GetNumLevels: c=%d, NumComp=%d, Progs[%d].Iso=%d, NumIso=%d\n", c, NumComp, m, Progs[m].Iso, NumIso);
                break;
            }
            LNum[Progs[m].Iso][c] += Progs[m].N;
        }
        if (n < NR) for (n=0; n < NumIso; ++n) for (m=0; m < NumComp; ++m) LNum[n][m] = 0;
    }
    if (statesList != 0) *statesList = StatesList;
}

int FitData::getNumLines(int* mv, int mJ)
{
    int n, R, J, N = fitDataCore->rowCount() - NMarkedLevels;
    for (n=R=0; n<N; n++) if ((J = fitDataCore->getJ(n)) <= mJ && fitDataCore->get_v(n) < mv[J]) R++;
    return R;
}

int FitData::getNumLines(int JD, int F, int v, int mJ, int Iso)
{
    int N = getNumLines(), n, R;
    for (n=R=0; n<N; n++) if (fitDataCore->getEnergy(n) != 0.0 && (mJ == 0 || fitDataCore->getJ(n) <= mJ)
            && (JD == 0 || JD == abs(fitDataCore->getJ(n) - fitDataCore->getJs(n)))
            && (Iso == -1 || fitDataCore->getIso(n) == Iso)
            && (F == -2 || (n < NSources ? FC[n] == F : F == -1))
             && (v >= 0 ? (v == fitDataCore->get_v(n)) : (v != -1 || fitDataCore->get_v(n) >= 0))) R++;
    return R;
}

int FitData::getNumLines() const
{
    int n, N=0, NR = fitDataCore->rowCount() - NMarkedLevels;
    for (n=0; n < NR; n++) if (fitDataCore->getEnergy(n) != 0.0) N++;
    return N;
}

int FitData::getNumProgressions(int *mv, int mJ)
{
    int n, m, l, b, R, J, N = getNumLines();
    if (N==0) return 0;
    int *S1 = heapSort(sortByProg), *SA = new int[N];
    for (n=0; n<N; n++) SA[S1[n]] = n;
    delete[] S1;
    LineTable *LTab = (NSources > 0 ? Sources[SA[0]] : 0);
    for (l = fitDataCore->getProgression(n), n=1, R=m=0; n<N; n++)
    {
        if (mJ > 0 && ((J = fitDataCore->getJ(SA[n])) > mJ || fitDataCore->get_v(SA[n]) > mv[J])) continue;
        if ((b = fitDataCore->getProgression(SA[n])) == l && (SA[n] >= NSources || LTab == Sources[SA[n]])) m++;
        else
        {
            if (m>0 && l>0)
            {
                R++;
                //printf("%d, ", n-m-1);
            }
            l=b;
            m=0;
            LTab = (SA[n] < NSources ? Sources[SA[n]] : 0);
        }
    }
    if (m>0) R++;
    //printf("\n");
    delete[] SA;
    return R;
}

void FitData::getSourceOffset(QStringList& Names, double*& Offsets)
{
    int n, m, N = fitDataCore->rowCount();
    QString Buffer;
    for (n=0; n<N; n++)
    {
        for (m=0, Buffer = fitDataCore->getSource(n).c_str(); m < Names.count() && Buffer != Names[m]; m++) ;
        if (m == Names.count() && !Buffer.isEmpty()) Names.append(Buffer);
        for (m=0, Buffer = fitDataCore->getSourceFile(n).c_str(); (m < Names.count() && Buffer != Names[m]); m++) ;
        if (m == Names.count() && !Buffer.isEmpty()) Names.append(Buffer);
    }
    if (Names.count() > 0) Offsets = new double[Names.count()];
    for (m=0; m < Names.count(); m++) Offsets[m] = 0.0;
    for (n=0; n < NSourceOffset; n++)
    {
        for (m=0; m < Names.count() && SourceOffsetNames[n] != Names[m]; m++) ;
        if (m < Names.count()) Offsets[m] = SourceOffset[n];
    }
}

QStringList FitData::getSources()
{
    int n, N = fitDataCore->rowCount();
    QStringList RList;
    QString Buffer;
    for (n=0; n<N; n++) if (!RList.contains(Buffer = fitDataCore->getSource(n).c_str())) RList.append(Buffer);
    return RList;
}

void FitData::getUncertaintyStats(QList< double >& Uncert, QList< int >& Numbers)
{
    int i, n, N = fitDataCore->rowCount(), vM = getvMax(), JM = getJMax();
    double U;
    for (n=0; n<N; n++) 
        if (fitDataCore->get_v(n) <= vM
            && fitDataCore->getJ(n) <= JM)
    {
        U = fitDataCore->getUncertainty(n);
        for (i=0; (i < Uncert.count() ? U > Uncert[i] : false); i++) ;
        if (i < Uncert.count() ? U == Uncert[i] : false) Numbers[i]++;
        else
        {
            Uncert.insert(i, U);
            Numbers.insert(i, 1);
        }
    }
}

double FitData::getUncertaintyOfvibLevel(int v, int I, Spektrum* Source)
{
    int cv = 1000, cI = 1000, av, aI, n, m=0, M, N = fitDataCore->rowCount(), dI, dv;
    double U;
    QList<double> Unc;
    QList<int> Num;
    for (n=0; n<N; n++) if (fitDataCore->getSource(n) == Source->getName().toStdString())
    {
        av = fitDataCore->get_v(n);
        aI = fitDataCore->getIso(n);
        if (abs(aI - I) < abs(cI - I) || (abs(aI - I) == abs(cI - I) 
            && abs(av - v) < abs(cv - v)))
        {
            cv = av;
            cI = aI;
        }
    }
    if (cv == 1000) return 0.0201;
    dI = abs(cI - I);
    dv = abs(cv - v);
    for (n=0; n<N; n++) 
        if (fitDataCore->getSource(n) == Source->getName().toStdString()
            && abs(fitDataCore->getIso(n) - I) <= dI
            && abs(fitDataCore->get_v(n) - v) <= dv)
    {
        for (m=0, U = fitDataCore->getUncertainty(n);
             m < Unc.count() && Unc[m] != U; m++) ;
        if (m < Unc.count()) Num[m]++;
        else
        {
            Unc.append(U);
            Num.append(1);
        }
    }
    for (n=M=0; n < Num.count(); n++) if (Num[n] > M) M = Num[m=n];
    return Unc[m];
}

void FitData::removeDataFSource()
{
    int n, m, N = fitDataCore->rowCount();
    QStringList S;
    QDialog *D;
    QGridLayout *L;
    QPushButton *O, *C;
    QComboBox *B;
    QString qb;
    std::string b;
    std::vector<int> list;
    int c;
    for (n=0; n<N; n++) if (!S.contains(qb = fitDataCore->getSource(n).c_str())) S.append(qb);
    D = new QDialog(this);
    D->setWindowTitle("MolSpektAnalysis");
    L = new QGridLayout(D);
    L->addWidget(new QLabel("Please select the source of which the data shall be removed from the table:", this), 0, 0, 1, 2);
    L->addWidget(new QLabel("Source:", this), 1, 0);
    L->addWidget(B = new QComboBox(this), 1, 1);
    L->setRowMinimumHeight(2, 20);
    L->addWidget(O = new QPushButton("OK", this), 3, 0);
    L->addWidget(C = new QPushButton("Cancel", this), 3, 1);
    B->addItems(S);
    B->setEditable(false);
    connect(O, SIGNAL(clicked()), D, SLOT(accept()));
    connect(C, SIGNAL(clicked()), D, SLOT(reject()));
    c = D->exec();
    b = B->currentText().toStdString();
    delete D;
    if (c == QDialog::Rejected) return;
    for (n=0; fitDataCore->getSource(n) != b; n++) ;
    table->blockSignals(true);
    for (m=n; m<N; m++)
    {
        if (fitDataCore->getSource(n) != b)
        {
            if (n < NSources)
            {
                if (m < NSources)
                {
                    Sources[n] = Sources[m];
                    FC[n] = FC[m];
                    LineElStates[n] = LineElStates[m];
                }
                else
                {
                    Sources[n] = 0;
                    FC[n] = -1;
                    LineElStates[n] = 0;
                }
            }
            n++;
        }
        else list.push_back(n);
    }
    NSources = n;
    fitDataCore->deleteRows(list.data(), list.size());
    fitDataCore->blockSignals(false);
}

void FitData::RemoveDoubled()
{
    int NR = fitDataCore->rowCount(), n, m, c;
    int *SA = new int[NR], *S1 = heapSort(sortIefJFreqv), *CompT, MCT;
    bool* toDel = new bool[NR];
    std::vector<int> delList;
    TermTable *TT = (State != 0 ? State->getTermTable() : 0);
    if (TT != 0) TT->getCompT(MCT, CompT);
    table->blockSignals(true);
    for (n=0; n < NR; n++)
    {
        SA[S1[n]] = n;
        toDel[n] = false;
    }
    for (m=0, n=1; n < NR; n++) 
    {
        if (fabs(fitDataCore->getEnergy(SA[m]) - fitDataCore->getEnergy(SA[n])) < 3e-3
            && fitDataCore->getIso(SA[m]) == fitDataCore->getIso(SA[n])
             && (fitDataCore->get_v(SA[m]) == fitDataCore->get_v(SA[n])
             || (stoi(fitDataCore->get_vs(SA[m])) < 0 && stoi(fitDataCore->get_vs(SA[n])) < 0))
             && fitDataCore->getJ(SA[m]) == fitDataCore->getJ(SA[n])
             && fitDataCore->getJs(SA[m]) == fitDataCore->getJs(SA[n])
             && (n >= NSources || FC[m] == FC[n] || FC[m] < 0 || FC[n] < 0 || FC[m] > MCT || FC[n] > MCT || CompT[FC[m]] < 0 || CompT[FC[n]] < 0))
        {
            if (fitDataCore->getSourceFile(SA[m]) == fitDataCore->getSourceFile(SA[n]) || fitDataCore->getSourceFile(SA[n]) == "")
            {
                toDel[n] = true;
                delList.push_back(n);
            }
            else 
            {
                if (fitDataCore->getSourceFile(SA[m]) == "")
                {
                    toDel[n] = true;
                    delList.push_back(n);
                }
                m=n;
            }
        }
        else m=n;
    }
    for (n=0; (n < NR && toDel[n] == false); n++) ;
    for (m=n+1; m < NR; m++) if (toDel[n] == false)
    {
        if (n < NSources)
        {
            if (m < NSources)
            {
                Sources[n] = Sources[m];
                FC[n] = FC[m];
                LineElStates[n] = LineElStates[m];
            }
            else
            {
                Sources[n] = 0;
                FC[n] = -1;
                LineElStates[n] = 0;
            }
        }
        n++;
    }
    if (n < NSources) NSources = n;
    fitDataCore->deleteRows(delList.data(), delList.size());
    table->blockSignals(false);
    Changed();
    delete[] SA;
    delete[] S1;
    delete[] toDel;
}

void FitData::removeMarkedLevel(TermEnergy& TE, Spektrum* Source)
{
    int n, N = fitDataCore->rowCount();
    table->blockSignals(true);
    for (n = N - NMarkedLevels; n<N; n++)
    {
        if (fitDataCore->getIso(n) != TE.Iso) continue;
        if (fitDataCore->get_v(n) != TE.v) continue;
        if (fitDataCore->getJ(n) != TE.J) continue;
        if (fitDataCore->getSource(n) != Source->getName().toStdString()) continue;
        if (fitDataCore->getSourceFile(n) != Source->getFileName().toStdString()) continue;
        if (fabs(fitDataCore->getEnergy(n) - TE.E) < 1e-4) break;
    }
    if (n==N)
    {
        printf("FitData::removeMarkedLevel error: Level cannot be found!");
        return;
    }
    fitDataCore->deleteRow(n);
    for (n+=1 ; n < NSources; n++)
    {
        FC[n-1] = FC[n];
        Sources[n-1] = Sources[n];
        LineElStates[n-1] = LineElStates[n];
    }
    if (NSources == N) NSources = N-1;
    else 
    {
        Sources[NSources - 1] = 0;
        FC[NSources - 1] = -1;
        LineElStates[NSources - 1] = 0;
    }
    NMarkedLevels--;
    table->blockSignals(false);
}

void FitData::removeSingleLines()
{
    int n, m, c, N = fitDataCore->rowCount(), lc;
    int *SA = new int[N], *S1 = heapSort(sortByProg);
    bool* del = new bool[N];
    memset(del, 0, sizeof(bool) * N);
    for (n=0; n<N; n++) SA[S1[n]] = n;
    delete[] S1;
    table->blockSignals(true);
    for (n = lc = 1; n < N - NMarkedLevels; n++, lc++) if (fitDataCore->getProgression(SA[n]) != fitDataCore->getProgression(SA[n-1]) || Sources[SA[n-1]] != Sources[SA[n]])
    {
        if (lc == 1 && fitDataCore->get_vs(SA[n-1]) != "TE" && stoi(fitDataCore->get_vs(SA[n-1])) >= -1) del[SA[n]] = true;
        lc = 0;
    }
    if (lc == 1) del[SA[n-1]] = true;
    for (n=0; n<N && !del[n]; n++) ;
    for (m=n+1; m<N; m++)
    {
        if (!del[m])
        {
            if (m < NSources)
            {
                LineElStates[n] = LineElStates[m];
                FC[n] = FC[m];
                Sources[n++] = Sources[m];
            }
            else if (n < NSources)
            {
                LineElStates[n] = 0;
                FC[n] = -1;
                Sources[n++] = 0;
            }
            else n++;
        }
        else fitDataCore->deleteRow(n);
    }
    table->blockSignals(false);
    if (n < NSources) NSources = n;
    Changed();
}

void FitData::selectDataFSource(QString Source)
{
    int n, t, N = fitDataCore->rowCount(), C = fitDataCore->columnCount();
    std::string source = Source.toStdString();
    QItemSelectionModel* model = new QItemSelectionModel;
    table->blockSignals(true);
    for (n=0; n < N; n++)
    {
        while (n<N && fitDataCore->getSource(n) != source) n++;
        t=n;
        while (n<N && fitDataCore->getSource(n) == source) n++;
        if (n>t)
        {
            QModelIndex topLeft = fitDataCore->getIndex(t, 0), bottomRight = fitDataCore->getIndex(n-1, C);
            QItemSelection newSelection(topLeft, bottomRight);
            model->select(newSelection, QItemSelectionModel::Select);
        }
    }
    table->setSelectionModel(model);
    table->blockSignals(false);
    emit SelChanged();
}

void FitData::addData(TableLine* Lines, int NLines)
{
    int n, m = fitDataCore->rowCount(), sigDig, i, j;
    TermTable *TT = (State != 0 ? State->getTermTable() : 0);
    int *CompZ = (TT != 0 ? TT->getCompZ() : 0), NFC = (TT != 0 ? TT->getNumComp() : 1);
    if (Sources != 0)
    {
        int *FCB = new int[m + NLines];
        LineTable **SB = new LineTable*[m + NLines];
        ElState **StB = new ElState*[m + NLines];
        for (n=0; n < NSources; n++)
        {
            FCB[n] = FC[n];
            SB[n] = Sources[n];
            StB[n] = LineElStates[n];
        }
        while (n<m)
        {
            FCB[n] = -1;
            StB[n] = 0;
            SB[n++] = 0;
        }
        delete[] Sources;
        delete[] FC;
        FC = FCB;
        Sources = SB;
        LineElStates = StB;
        NSources = m + NLines;
    }
    else
    {
        Sources = new LineTable*[NSources = NLines + m];
        FC = new int[NSources];
        LineElStates = new ElState*[NSources];
        for (n=0; n<m; n++)
        {
            FC[n] = -1;
            Sources[n] = 0;
            LineElStates[n] = 0;
        }
    }
    table->blockSignals(true);
    fitDataCore->setRowCount(NSources);
    for (n=0; n < NLines; n++, m++)
    {
        FC[m] = (CompZ != 0 && Lines[n].FC >= 0 && Lines[n].FC < NFC ? CompZ[Lines[n].FC] : Lines[n].FC);
        fitDataCore->setIso(m, Lines[n].Iso);
        fitDataCore->set_v(m, Lines[n].vss);
        fitDataCore->setJ(m, Lines[n].Jss);
        fitDataCore->set_vs(m, (Lines[n].vs != -1 ? std::to_string(Lines[n].vs) : (Lines[n].isTE ? "TE" : "nA")));
        fitDataCore->setJs(m, Lines[n].Js);
        if ((Sources[n] = Lines[n].LTab) != 0) fitDataCore->setSource(m, Sources[n]->getName().toStdString());
        else if (!Lines[n].SourceName.isEmpty()) fitDataCore->setSource(m, Lines[n].SourceName.toStdString());
        else
        {
            i = Lines[n].File.lastIndexOf(QRegExp("[\\/]")) + 1;
            j = Lines[n].File.indexOf('.', i);
            fitDataCore->setSource(m, Lines[n].File.mid(i, j-i).toStdString());
        }
        LineElStates[m] = Lines[n].State;
        if (LineElStates[m] == 0 && Sources[m] != 0) LineElStates[m] = Sources[m]->getElState();
        fitDataCore->setProgression(m, Lines[n].PN);
        fitDataCore->setSourceFile(m, Lines[n].File.toStdString());
        sigDig = -int(floor(log10(Lines[n].err)));
        if (sigDig < 4) sigDig = 4;
        fitDataCore->setEnergy(m, Lines[n].WN);
        fitDataCore->setUncertainty(m, Lines[n].err);
        fitDataCore->setObsCalc(m, Lines[n].dev);
        fitDataCore->setDevRatio(m, Lines[n].DevR);
    }
    table->blockSignals(false);
    delete[] Lines;
    Changed();
}

void FitData::setData(TableLine* Lines, int NLines)
{
    table->blockSignals(true);
    if (Sources != 0)
    {
        delete[] Sources;
        delete[] FC;
        delete[] LineElStates;
        NSources = 0;
        Sources = 0;
        LineElStates = 0;
        FC = 0;
    }
    fitDataCore->setRowCount(0);
    addData(Lines, NLines);
    setNewCreated();
}

void FitData::setDev(double* dev, double* DevR)
{
    int n, m, N = fitDataCore->rowCount(), mJ = getJMax(), mv = getvMax();
    QString Buffer;
    table->blockSignals(true);
    for (n=m=0; n<N; n++) if (fitDataCore->get_v(n) <= mv && fitDataCore->getJ(n) <= mJ)
    {
        fitDataCore->setObsCalc(n, dev[m]);
        fitDataCore->setDevRatio(n, DevR[m++]);
    }
    table->blockSignals(false);
    Changed();
}

void FitData::setDev(double* dev, int* RowN, int N)
{
    int n;
    QString Buffer;
    table->blockSignals(true);
    for (n=0; n<N; n++)
    {
        fitDataCore->setObsCalc(RowN[n], dev[n]);
        fitDataCore->setDevRatio(RowN[n], dev[n] / fitDataCore->getUncertainty(RowN[n]));
    }
    table->blockSignals(false);
    Changed();
}

void FitData::setDev(TableLine **TL, TLRef *SortArray, int NE, double tol)
{
    int *inSort = heapSort(sortIefJvFreq), i, e, N = getNumLines(), *intSort = new int[N], I, DJ, J, v;
    double E, ObsCalc;
    table->blockSignals(true);
    for (i=0; i<N; ++i)
    {
        fitDataCore->setObsCalc(i, 0.0);
        fitDataCore->setDevRatio(i, 0.0);
        intSort[inSort[i]] = i;
    }
    for (i=e=0; i<N; ++i)
    {
        I = fitDataCore->getIso(intSort[i]);
        J = fitDataCore->getJ(intSort[i]);
        DJ = J - fitDataCore->getJs(intSort[i]);
        v = fitDataCore->get_v(intSort[i]);
        E = fitDataCore->getEnergy(intSort[i]);
        while (e < NE - 1 && (TL[SortArray[e].state][SortArray[e].line].Iso < I
               || (TL[SortArray[e].state][SortArray[e].line].Iso == I && TL[SortArray[e].state][SortArray[e].line].Jss - TL[SortArray[e].state][SortArray[e].line].Js > DJ)
               || (TL[SortArray[e].state][SortArray[e].line].Iso == I && TL[SortArray[e].state][SortArray[e].line].Jss - TL[SortArray[e].state][SortArray[e].line].Js == DJ
                   && TL[SortArray[e].state][SortArray[e].line].Jss < J)
               || (TL[SortArray[e].state][SortArray[e].line].Iso == I && TL[SortArray[e].state][SortArray[e].line].Jss - TL[SortArray[e].state][SortArray[e].line].Js == DJ
                   && TL[SortArray[e].state][SortArray[e].line].Jss == J && TL[SortArray[e].state][SortArray[e].line].vss < v)
               || (TL[SortArray[e].state][SortArray[e].line].Iso == I && TL[SortArray[e].state][SortArray[e].line].Jss - TL[SortArray[e].state][SortArray[e].line].Js == DJ
                   && TL[SortArray[e].state][SortArray[e].line].Jss == J && TL[SortArray[e].state][SortArray[e].line].vss == v
                   && E - TL[SortArray[e].state][SortArray[e].line].WN > tol))) ++e;
        if (TL[SortArray[e].state][SortArray[e].line].Iso == I && TL[SortArray[e].state][SortArray[e].line].Jss - TL[SortArray[e].state][SortArray[e].line].Js == DJ
                && TL[SortArray[e].state][SortArray[e].line].Jss == J
                && abs(TL[SortArray[e].state][SortArray[e].line].WN - E) < tol)
        {
            ObsCalc = E - TL[SortArray[e].state][SortArray[e].line].WN + TL[SortArray[e].state][SortArray[e].line].dev;
            fitDataCore->setObsCalc(intSort[i], ObsCalc);
            fitDataCore->setDevRatio(intSort[i], ObsCalc / fitDataCore->getUncertainty(intSort[i]));
            TL[SortArray[e].state][SortArray[e].line].isSelected = true;
        }
    }
    table->blockSignals(false);
    delete[] inSort;
    delete[] intSort;
    Changed();
}

void FitData::setElState(ElState* nState)
{
    State = nState;
    for (int n=0; n < NSources; ++n) if (0 == LineElStates[n])
    {
        LineElStates[n] = nState;
        fitDataCore->setSecondState(n, nState->getName().toStdString());
    }
}

void FitData::setMolecule(Molecule *Mol)
{
    TableWindow::setMolecule(Mol);
    int n;
    for (n=0; n < NSources; n++) LineElStates[n] = Mol->getState(fitDataCore->getOtherState(n).c_str());
    for (QList<ResidualFit*>::iterator it = residualFits.begin(); it != residualFits.end(); ++it)
        (*it)->SetState(Mol->getState(*(*it)->getStateName()));
}

void FitData::setFC(int nFC)
{
    int n, N, *Rows;
    table->getSelectedRows(Rows, N);
    if (N > 0)
    {
        for (n=0; n<N; n++) if (Rows[n] < NSources) FC[Rows[n]] = nFC;
    }
    else for (n=0; n < NSources; n++) FC[n] = nFC;
}

void FitData::setRWErr(double* RWErr)
{
    table->blockSignals(true);
    int n, N = fitDataCore->rowCount();
    QString Buffer;
    for (n=0; n<N; n++)
    {
        Buffer = QString::number(fitDataCore->getUncertainty(n), 'f', 4);
        fitDataCore->setSecondState(n, QString::number(RWErr[n], 'f', Buffer.length() - Buffer.indexOf('.') - 1).toStdString());
    }
    fitDataCore->setRWError("RWErr");
    table->blockSignals(false);
}

void FitData::setSourceOffset(QStringList& Names, double* Offsets)
{
    int n, m;
    if (NSourceOffset > 0)
    {
        delete[] SourceOffset;
        delete[] SourceOffsetNames;
    }
    for (n=m=0; n < Names.count(); n++) if (Offsets[n] != 0.0) m++;
    if ((NSourceOffset = m) > 0)
    {
        SourceOffset = new double[m];
        SourceOffsetNames = new QString[m];
    }
    else
    {
        SourceOffset = 0;
        SourceOffsetNames = 0;
    }
    for (n=m=0; n < Names.count(); n++) if (Offsets[n] != 0.0)
    {
        SourceOffset[m] = Offsets[n];
        SourceOffsetNames[m++] = Names[n];
    }
    if (Names.count() > 0) delete[] Offsets;
    Changed();
}

void FitData::setUncertainty(double Uncertainty, bool Min)
{
    table->blockSignals(true);
    int n, m, i, N, *Rows;
    double cU;
    table->getSelectedRows(Rows, N);
    double U = Uncertainty;
    double mU = 9.0 + Uncertainty, MU = 99.0 + Uncertainty, nU;
    QList<LineTable*> SL;
    QList<int> CL;
    if (N > 0)
    {
        for (n=0; n<N; n++) if (Rows[n] < NSources && Sources[Rows[n]] != 0 && Sources[Rows[n]]->getTransition()->getLowerState() == State)
        {
            for (i=0; (i < SL.count() ? SL[i] != Sources[Rows[n]] : false); i++) ;
            if (i == SL.count())
            {
                SL.append(Sources[Rows[n]]);
                CL.append(1);
            }
            else CL[i]++;
        }
        int *PN[SL.count()], *vss[SL.count()], *Jss[SL.count()];
        QString *Err[SL.count()];
        for (n=0; n < SL.count(); n++)
        {
            PN[n] = new int[CL[n]];
            vss[n] = new int[CL[n]];
            Jss[n] = new int[CL[n]];
            Err[n] = new QString[CL[n]];
            CL[n] = 0;
        }
        for (n=0; n<N; n++)
        {
            
            if (!Min || (cU = fitDataCore->getUncertainty(Rows[n])) < Uncertainty) fitDataCore->setUncertainty(Rows[n], nU = U);
            else if (cU >= 99.0) fitDataCore->setUncertainty(Rows[n], nU = MU);
            else if (cU >= 9.0) fitDataCore->setUncertainty(Rows[n], nU = mU);
            if (Rows[n] < NSources && Sources[Rows[n]] != nullptr && Sources[Rows[n]]->getTransition()->getLowerState() == State)
            {
                for (i=0; SL[i] != Sources[Rows[n]]; i++) ;
                PN[i][CL[i]] = fitDataCore->getProgression(Rows[n]);
                vss[i][CL[i]] = fitDataCore->get_v(Rows[n]);
                Jss[i][CL[i]] = fitDataCore->getJ(Rows[n]);
                Err[i][CL[i]++] = nU;
            }
        }
        for (n=0; n < SL.count(); n++)
        {
            SL[n]->SetError(CL[n], PN[n], vss[n], Jss[n], Err[n]);
            delete[] PN[n];
            delete[] vss[n];
            delete[] Jss[n];
            delete[] Err[n];
        }
    }
    else 
    {
        for (m=0; m < NSources; m++) if (Sources[m]->getTransition()->getLowerState() == State)
        {
            for (i=0; (i < SL.count() ? SL[i] != Sources[m] : false); i++) ;
            if (i == SL.count())
            {
                SL.append(Sources[m]);
                CL.append(1);
            }
            else CL[i]++;
        }
        int *PN[SL.count()], *vss[SL.count()], *Jss[SL.count()];
        QString *Err[SL.count()];
        for (n=0; n < SL.count(); n++)
        {
            PN[n] = new int[CL[n]];
            vss[n] = new int[CL[n]];
            Jss[n] = new int[CL[n]];
            Err[n] = new QString[CL[n]];
            CL[n] = 0;
        }
        for (m=0; m < fitDataCore->rowCount(); m++)
        {
            if (!Min || (cU = fitDataCore->getUncertainty(m)) < Uncertainty) fitDataCore->setUncertainty(m, nU = U);
            else if (cU >= 99.0) fitDataCore->setUncertainty(m, nU = MU);
            else if (cU >= 9.0) fitDataCore->setUncertainty(m, nU = mU);
            if (m < NSources) if (Sources[m]->getTransition()->getLowerState() == State)
            {
                for (i=0; SL[i] != Sources[m]; i++) ;
                PN[i][CL[i]] = fitDataCore->getProgression(m);
                vss[i][CL[i]] = fitDataCore->get_v(m);
                Jss[i][CL[i]] = fitDataCore->getJ(m);
                Err[i][CL[i]++] = nU;
            }
        }
        for (n=0; n < SL.count(); n++)
        {
            SL[n]->SetError(CL[n], PN[n], vss[n], Jss[n], Err[n]);
            delete[] PN[n];
            delete[] vss[n];
            delete[] Jss[n];
            delete[] Err[n];
        }
    }
    table->blockSignals(false);
    Changed();
}

void FitData::setUncertainty(int* Rows, double* Uncertainties, int N)
{
    int n, i;
    table->blockSignals(true);
    QList<LineTable*> SL;
    QList<int> CL;
    for (n=0; n<N; n++) if (Rows[n] < NSources && Sources[Rows[n]] != 0)
        if (Sources[Rows[n]]->getTransition()->getLowerState() == State)
    {
        for (i=0; i < SL.count() && SL[i] != Sources[Rows[n]]; i++) ;
        if (i == SL.count())
        {
            SL.append(Sources[Rows[n]]);
            CL.append(1);
        }
        else CL[i]++;
    }
    int *PN[SL.count()], *vss[SL.count()], *Jss[SL.count()];
    double *Err[SL.count()];
    for (n=0; n < SL.count(); n++)
    {
        PN[n] = new int[CL[n]];
        vss[n] = new int[CL[n]];
        Jss[n] = new int[CL[n]];
        Err[n] = new double[CL[n]];
        CL[n] = 0;
    }
    for (n=0; n<N; n++) 
    {
        fitDataCore->setUncertainty(Rows[n], Uncertainties[n]);
        fitDataCore->setDevRatio(Rows[n], fitDataCore->getObsCalc(Rows[n]) / Uncertainties[n]);
        if (Rows[n] < NSources && Sources[Rows[n]] != 0 && Sources[Rows[n]]->getTransition()->getLowerState() == State)
        {
            for (i=0; SL[i] != Sources[Rows[n]]; i++) ;
            PN[i][CL[i]] = fitDataCore->getProgression(Rows[n]);
            vss[i][CL[i]] = fitDataCore->get_v(Rows[n]);
            Jss[i][CL[i]] = fitDataCore->getJ(Rows[n]);
            Err[i][CL[i]++] = Uncertainties[n];
        }
    }
    for (n=0; n < SL.count(); n++)
    {
        SL[n]->SetError(CL[n], PN[n], vss[n], Jss[n], Err[n]);
        delete[] PN[n];
        delete[] vss[n];
        delete[] Jss[n];
        delete[] Err[n];
    }
    table->blockSignals(false);
    Changed();
}

void FitData::setWeightFactors(int N, int* Rows, double* Factors)
{
    table->blockSignals(true);
    int n, NR = fitDataCore->rowCount();
    std::string Buffer = "1.0";
    fitDataCore->setRWError("WFact");
    for (n=0; n < NR; n++) fitDataCore->setSecondState(n, Buffer);
    for (n=0; n<N; n++) fitDataCore->setUncertainty(Rows[n], Factors[n]);
    table->blockSignals(false);
    delete[] Rows;
    delete[] Factors;
    Changed();
}

void FitData::setMCSEnergies(double* MCSE, int* Rows, int NR)
{
    int n;
    table->blockSignals(true);
    for (n=0; n < NR; n++) fitDataCore->setEnergy(Rows[n], MCSE[n]);
    table->blockSignals(false);
    delete[] MCSE;
    delete[] Rows;
}

void FitData::sortTab(int* S2)
{
    int i, P1=0, n, N = Tab->rowCount(), AFC[2];
    BaseData *AIt[2];
    LineTable *AS[2];
    ElState *ASt[2];
    Tab->blockSignals(true);
    for (i=0; i<N; i++) if (S2[i] != i)
    {
        AIt[P1] = fitDataCore->getData(i);
        if (i < NSources)
        {
            AFC[P1] = FC[i];
            AS[P1] = Sources[i];
            ASt[P1] = LineElStates[i];
        }
        else
        {
            AFC[P1] = 0;
            AS[P1] = 0;
            ASt[P1] = 0;
        }
        while (S2[i] != i) 
        {
            AIt[1-P1] = fitDataCore->getData(S2[i]);
            fitDataCore->setData(S2[i], AIt[P1]);
            if (S2[i] < NSources)
            {
                AFC[1-P1] = FC[S2[i]];
                AS[1-P1] = Sources[S2[i]];
                ASt[1-P1] = LineElStates[S2[i]];
                FC[S2[i]] = AFC[P1];
                Sources[S2[i]] = AS[P1];
                LineElStates[S2[i]] = ASt[P1];
            }
            else
            {
                AFC[1-P1] = 0;
                AS[1-P1] = 0;
                ASt[1-P1] = 0;
            }
            P1 = 1 - P1;
            n = S2[S2[i]];
            S2[S2[i]] = S2[i];
            S2[i] = n;
        }
        fitDataCore->setData(i, AIt[P1]);
        if (i < NSources)
        {
            FC[i] = AFC[P1];
            Sources[i] = AS[P1];
            LineElStates[i] = ASt[P1];
        }    
    }
    table->blockSignals(false);
    delete[] S2;
    Changed();
}

void FitData::sortByLTabAndProg()
{
    SortByLTabAndProgFunctor sortFunct(fitDataCore, Sources, NSources, MW, molecule);
    sortTab(utils::heapSort(sortFunct, fitDataCore->rowCount()));
}

void FitData::updateData()
{
    if (State == 0 || molecule == 0)
    { 
        printf("FitData::updateData error: This object is not assigned to a molecule!\n");
        return;
    }
    TSDialog Diag(this, molecule, State);
    if (Diag.exec() == QDialog::Rejected) return;
    int rr, r = 0, l, n, m, N = Diag.getNumSelected(), NL[N], MJ = 0, Mv = 0, *FCB, NR = fitDataCore->rowCount();
    if (N == 0 && !Diag.getUpdateSA()) return;
    LineTable *L[N], **SB;
    ElState **StB;
    TermEnergy *E;
    TableLine *TL;
    QString CSource;
    bool *Del = new bool[NR];
    memset(Del, 0, NR * sizeof(bool));
    for (n=0; n<N; n++) L[n] = Diag.getSelected(n);
    table->blockSignals(true);
    if (Diag.getKeepOther())
    {
        for (n=0; n < NSources; n++)
        {
            CSource = fitDataCore->getSource(n).c_str();
            for (l=0; l<N && Sources[n] != L[l] && CSource != L[l]->getName(); l++) ;
            if (l<N) Del[n] = true;
        }
        if (Diag.getUpdateSA() && molecule != 0 && State != molecule->getStateP(0))
        {
            int NS = 0, NumStates = molecule->getNumStates() - 2;
            TermTable *TT = State->getTermTable();
            int NP = (TT != 0 ? TT->getNPerturbations() : 0), Numv = TT->getMaxv() + 1, nC = TT->getNumComp();
            int NumJ = TT->getMaxJ() + 1, NumI = TT->getNumIso();
            double D, sD, En, ****TData = TT->getData();
            if (NP > 0)
            {
                int NComp[NumStates], Nv[NumStates], s, c, I=0, J=0, v=0, bv = 0, bs, R;
                bool *PComp[NumStates];
                TermTable *PTT[NumStates];
                FitData *PFD[NumStates];
                Perturbation *Pert;
                double ****PTData[NumStates];
                QPixmap Pix(10, 10);
                QPainter P(&Pix);
                P.setPen(QColor(255, 0, 0));
                P.setFont(QFont("Arial", 10));
                P.drawText(0, 10, "N");
                for (s=-1; s < NS; s++)
                {
                    if (s>=0) NP = PTT[s]->getNPerturbations();
                    for (n=0; n < NP; n++)
                    {
                        if (s==-1) Pert = TT->getPerturbation(n);
                        else Pert = PTT[s]->getPerturbation(n);
                        if (Pert->Perturber == TT || Pert->Perturber == 0 || !PComp[Pert->Comp]) continue;
                        for (m=0; m < NS && Pert->Perturber != PTT[m]; m++) ;
                        if (m == NS)
                        {
                            PTT[m] = Pert->Perturber;
                            if (PTT[m]->getElState() == 0) continue;
                            NComp[m] = PTT[m]->getNumComp();
                            Nv[m] = PTT[m]->getMaxv() + 1;
                            PComp[m] = new bool[NComp[m]];
                            for (c=0; c < NComp[m]; c++) PComp[m][c] = false;
                            PTData[m] = PTT[m]->getData();
                            PFD[NS++] = 0;
                        }
                        PComp[m][Pert->PComp] = true;
                    }
                }
                for (r=0; r < NR; r++) if (fitDataCore->getSource(r).empty() && (I = fitDataCore->getIso(r)) < NumI && (v = fitDataCore->get_v(r)) < Numv && (J = fitDataCore->getJ(r)) < NumJ)
                {
                    c = (fitDataCore->getJs(r) != J ? 0 : (nC > 1 ? 1 : 0));
                    for (sD = fabs((En = fitDataCore->getEnergy(r)) - TData[c][I][v][J]), s=0, bs = -1; s < NS; s++)
                        for (c=0; c < NComp[s]; c++) if (PComp[s][c]) for (v=0; v < Nv[s]; v++)
                            if ((D = fabs(En - PTData[s][c][I][v][J])) < sD)
                    {
                        sD = D;
                        bs = s;
                        bv = v;
                    }
                    if (bs >= 0)
                    {
                        if (PFD[bs] == 0)
                        {
                            if ((PFD[bs] = PTT[bs]->getElState()->getFitData()) == 0)
                            {
                                if ((PFD[bs] = MW->CreateFitData()) == 0) continue;
                                PTT[bs]->getElState()->addFitData(PFD[bs]);
                            }
                            if (!PFD[bs]->isVisible()) PFD[bs]->show();
                            PFD[bs]->activateWindow();
                            PFD[bs]->setFocus();
                        }
                        R = PFD[bs]->fitDataCore->rowCount();
                        PFD[bs]->NMarkedLevels++;
                        PFD[bs]->table->blockSignals(true);
                        PFD[bs]->fitDataCore->setRowCount(R+1);
                        PFD[bs]->fitDataCore->setData(R, fitDataCore->getData(r));
                        PFD[bs]->fitDataCore->set_v(R, bv);
                        PFD[bs]->table->blockSignals(false);
                        PFD[bs]->setBlockChangeSignal(true);
                        PFD[bs]->Changed();
                        PFD[bs]->setBlockChangeSignal(false);
                        Del[r] = true;
                    }
                }
                for (s=0; s < NS; s++) delete[] PComp[s];
            }
        }
        for (r=0, NR = fitDataCore->rowCount(); r < NR && !Del[r]; r++) ;
        for (l=r+1; l < NR; l++) if (!Del[l])
        {
            fitDataCore->setData(r, fitDataCore->getData(l));
            if (l < NSources)
            {
                Sources[r] = Sources[l];
                FC[r] = FC[l];
                LineElStates[r] = LineElStates[l];
            }
            else if (r < NSources)
            {
                Sources[r] = 0;
                FC[r] = -1;
                LineElStates[r] = 0;
            }
            r++;
        }
        if (r < NSources) NSources = r;
    }
    else
    {
        NSources = 0;
        delete[] Sources;
        Sources = 0;
        delete[] FC;
        FC = 0;
        delete[] LineElStates;
        LineElStates = 0;
    }
    fitDataCore->setRowCount(rr = r);
    for (n=0; n<N; n++) 
    {
        if (L[n]->getTransition()->getLowerState() != State)
        {
            L[n]->getgoodTE(NL[n], E);
            fitDataCore->setRowCount(r + NL[n]);
            if (FC != 0)
            {
                FCB = new int[r + NL[n]];
                StB = new ElState*[r + NL[n]];
                for (m=0; m<r; m++)
                {
                    FCB[m] = FC[m];
                    StB[m] = LineElStates[m];
                }
                delete[] FC;
                delete[] LineElStates;
                FC = FCB;
                LineElStates = StB;
            }
            else FC = new int[r + NL[n]];
            for (l=0; l < NL[n]; l++) 
            {
                FC[r] = E[l].FC;
                LineElStates[r] = L[n]->getTransition()->getUpperState();
                fitDataCore->setIso(r, E[l].Iso);
                fitDataCore->set_v(r, E[l].v);
                fitDataCore->setJ(r, E[l].J);
                fitDataCore->set_vs(r, "TE");
                fitDataCore->setJs(r, (E[l].ef ? E[l].J + 1 : fitDataCore->getJ(r)));
                fitDataCore->setSource(r, L[n]->getName().toStdString());
                fitDataCore->setProgression(r, E[l].PN);
                fitDataCore->setSourceFile(r, E[l].File.toStdString());
                fitDataCore->setEnergy(r, E[l].E);
                fitDataCore->setUncertainty(r, E[l].err);
                fitDataCore->setObsCalc(r, E[l].dev);
                fitDataCore->setDevRatio(r, E[l].DevR);
                fitDataCore->setSecondState(r, LineElStates[r]->getName().toStdString());
                ++r;
                if (E[l].v > Mv) Mv = E[l].v;
                if (E[l].J > MJ) MJ = E[l].J;
            }
            delete[] E;
        }
        else
        {
            L[n]->getgoodLines(NL[n], TL);
            fitDataCore->setRowCount(r + NL[n]);
            if (FC != 0)
            {
                FCB = new int[r + NL[n]];
                StB = new ElState*[r + NL[n]];
                for (m=0; m<r; m++)
                {
                    FCB[m] = FC[m];
                    StB[m] = LineElStates[m];
                }
                delete[] FC;
                delete[] LineElStates;
                FC = FCB;
                LineElStates = StB;
            }
            else FC = new int[r + NL[n]];
            for (l=0; l < NL[n]; l++)
            {
                FC[r] = TL[l].FC;
                LineElStates[r] = L[n]->getTransition()->getLowerState();
                fitDataCore->setIso(r, TL[l].Iso);
                fitDataCore->set_v(r, TL[l].vss);
                fitDataCore->setJ(r, TL[l].Jss);
                fitDataCore->set_vs(r, TL[l].vs != -1 ? to_string(TL[l].vs) : "nA");
                fitDataCore->setJs(r, TL[l].Js);
                fitDataCore->setSource(r, L[n]->getName().toStdString());
                fitDataCore->setProgression(r, TL[l].PN);
                fitDataCore->setSourceFile(r, TL[l].File.toStdString());
                fitDataCore->setEnergy(r, TL[l].WN);
                fitDataCore->setUncertainty(r, TL[l].err);
                fitDataCore->setObsCalc(r, TL[l].dev);
                fitDataCore->setDevRatio(r++, TL[l].DevR);
                if (TL[l].vss > Mv) Mv = TL[l].vss;
                if (TL[l].Jss > MJ) MJ = TL[l].Jss;
            }
            delete[] TL;
        }
    }
    table->blockSignals(false);
    SB = Sources;
    Sources = new LineTable*[r];
    for (n=0; n < NSources; n++) Sources[n] = SB[n];
    if (SB != 0) delete[] SB;
    while (n < rr) Sources[n++] = 0;
    for (n=0, NSources = r; n<N; n++) for (l=0; l < NL[n]; l++) Sources[rr++] = L[n];
    if (getvMax() < Mv) setvMax(Mv);
    if (getJMax() < MJ) setJMax(MJ);
    Changed();
}

void FitData::updateEnergy(int N, int *r, double *Energy)
{
    int n;
    table->blockSignals(true);
    for (n=0; n<N; n++) fitDataCore->setEnergy(r[n], Energy[n]);
    table->blockSignals(false);
    Changed();
}

void FitData::updateLevels(int N, int *r, double *Energy, int *vs, int *vss)
{
    int n;
    table->blockSignals(true);
    for (n=0; n<N; n++) 
    {
        fitDataCore->setEnergy(r[n], Energy[n]);
        fitDataCore->set_v(r[n], vss[n]);
        fitDataCore->set_vs(r[n], (vs[n] == -1 ? "TE" : QString::number(vs[n]).toStdString()));
    }
    table->blockSignals(false);
    Changed();
}

void FitData::updateRow(TableLine* Line)
{
    int n, N = fitDataCore->rowCount();
    for (n=0; n<N; n++) if (fitDataCore->getEnergy(n) == Line->WN && fitDataCore->getIso(n) == Line->Iso && fitDataCore->get_v(n) == Line->vss && fitDataCore->getJ(n) == Line->Jss)
        fitDataCore->setUncertainty(n, Line->err);
}

void FitData::writeData(QTextStream& S)
{
    fitDataCore->writeData(S);
}

bool FitData::writeData(QString Filename)
{
    bool Success;
    int n, N = fitDataCore->rowCount(), NSources = fitDataCore->getNSources();
    std::string B;
    if (NMarkedLevels > 0)
    {
        QPixmap P;
        LineTable **NewSources = new LineTable*[N];
        int *NFC = new int[N];
        ElState **NSt = new ElState*[N];
        for (n=0; n < NSources; n++)
        {
            NewSources[n] = Sources[n];
            NFC[n] = FC[n];
            NSt[n] = LineElStates[n];
        }
        delete[] Sources;
        delete[] FC;
        delete[] LineElStates;
        Sources = NewSources;
        FC = NFC;
        LineElStates = NSt;
        table->blockSignals(true);
        for (n = NSources; n<N; n++)
        {
            Sources[n] = nullptr;
            LineElStates[n] = (molecule != 0 && C > fdcLineElState ? molecule->getState(fitDataCore->getOtherState(n).c_str()) : nullptr);
            B = fitDataCore->get_vs(n);
            FC[n] = (B == "TE" ? -1 : -1 - stoi(B));
        }
        table->blockSignals(false);
        NMarkedLevels = 0;
        Changed();
        emit AssignmentsAccepted(this);
    }
    table->blockSignals(true);
    for (n=0; n < NSources; n++) if (FC[n] >= 0)
    {
        B = fitDataCore->get_vs(n), ch = to_string(B == "TE" ? -1 : (B == "nA" ? -10 : stoi(B))) + 1000 * (FC[n] + 1);
        fitDataCore->set_vs(n, ch);
    }
    QString MolPath = molecule->getFileName();
    if (molecule != 0) for (n=0; n<N; ++n)
    {
        QString CurPath = fitDataCore->getSourceFile(n).c_str();
        std::string relPath = getRelativePath(CurPath, MolPath).toStdString();
        fitDataCore->setSourceFile(n, relPath);
    }
    Success = TableWindow::writeData(Filename);
    for (n=0; n < NSources; n++) if (FC[n] >= 0)
    {
        std::string vs = (c = stoi(fitDataCore->get_vs(n)) - 1000 * (FC[n] + 1)) == -1 ? "TE" : (c == -10 ? "nA" : to_string(c));
        fitDataCore->set_vs(n, vs);
    }
    if (molecule != 0) for (n=0; n<N; ++n)
    {
        QString CurPath = fitDataCore->getSourceFile(n).c_str();
        std::string absPath = getAbsolutePath(CurPath, MolPath);
        fitDataCore->setSourceFile(n, absPath);
    }
    table->blockSignals(false);
    if ((NSourceOffset > 0 || residualFits.size() > 0) && Success)
    {
        QFile F(getFileName());
        F.open(QIODevice::Append);
        QTextStream S(&F);
        if (NSourceOffset > 0)
        {
            S << "SourceOffsets:\n";
            for (n=0; n < NSourceOffset; n++) S << SourceOffsetNames[n] << '\t' << QString::number(SourceOffset[n], 'f', 12) << '\n';
        }
        for (QList<ResidualFit*>::const_iterator it = residualFits.begin(); it != residualFits.end(); ++it) (*it)->writeData(&S);
    }
    return Success;
}

bool FitData::writeExPotFitInput(QString Filename)
{
    if (Filename.isEmpty()) return false;
    if (molecule == 0) return false;
    IsoTab *Iso = molecule->getIso();
    if (Iso == 0) return false;
    QFile File(Filename);
    int nIso = Iso->numIso, n, N = Tab->rowCount(), I, sI, MsI = 0, *StA = 0, SourceI = 0;
    QString IsoStr[nIso], B;
    ElState *stateArray[molecule->getNumStates()];
    for (sI=0; sI < molecule->getNumStates(); ++sI) stateArray[sI] = 0;
    if (LineElStates != 0) for (n=0; n < NSources; ++n)
    {
        for (sI = 0; stateArray[sI] != LineElStates[n] && stateArray[sI] != 0; ++sI) ;
        if (stateArray[sI] == 0) stateArray[sI] = LineElStates[n];
        if (sI > MsI) MsI = sI;
    }
    if (MsI > 0)
    {
        EsPotFitInputElStateAssignDialog D(this, stateArray, MsI + 1);
        if (D.exec() == QDialog::Rejected)
        {
            delete Iso;
            return true;
        }
        if (!D.IsAssignmentFromTabSelected()) StA = D.GetAssignment();
    }
    if (!File.open(QIODevice::WriteOnly))
    {
        delete Iso;
        return false;
    }
    QTextStream WS(&File);
    for (n=0; n < nIso; n++) 
        IsoStr[n] = ("    " + QString::number(Iso->mNumIso1[n])).right(5)
                  + ("    " + QString::number(Iso->mNumIso2[n])).right(5);
    QString FileName;
    for (n=0; n<N; n++)
    {
        I = Tab->item(n, 0)->text().toInt();
        B = Tab->item(n, 3)->text();
        if (NSourceOffset > 0)
        {
            FileName = Tab->item(n, fdcFile)->text();
            for (SourceI = 0; SourceI < NSourceOffset && FileName.indexOf(SourceOffsetNames) == -1; ++SourceI) ;
        }
        if (StA != 0 && n < NSources) for (sI = 0; stateArray[sI] != LineElStates[n]; ++sI) ;
        WS << IsoStr[I] + ("     " + Tab->item(n, 1)->text()).right(5) 
              + ("     " + Tab->item(n, 2)->text()).right(5) 
              + ("     " + (StA != 0 && n < NSources ? QString::number(StA[sI]) :
                                                       (B.indexOf("TE") >= 0 || B.indexOf("nA") >= 0 ? "-1" : B))).right(5)
              + ("     " + Tab->item(n, 4)->text()).right(5) + IsoStr[I]
              + "    0    0    0    0\n"
              + ("               "
                        + (SourceI < NSourceOffset ? QString::number(Tab->item(n, 8)->text().toDouble() - SourceOffset[SourceI], 'f', 4) : Tab->item(n, 8)->text())).right(15)
              + ("               " + Tab->item(n, 9)->text()).right(15) + "    2\n";
    }
    File.close();
    delete Iso;
    if (StA != 0) delete[] StA;
    return true;
}

bool FitData::writeTFGS(QString Filename)
{
    if (Filename.isEmpty()) return false;
    if (molecule == 0) return false;
    IsoTab *Iso = molecule->getIso();
    if (Iso == 0) return false;
    Atom *atom1 = molecule->getAtom1(), *atom2 = molecule->getAtom2();
    if (atom1 == 0 || atom2 == 0) return false;
    QFile File(Filename);
    if (!File.open(QIODevice::WriteOnly))
    {
        delete Iso;
        return false;
    }
    int nIso = Iso->numIso, n, N = Tab->rowCount(), I, vo = 0, *S1 = heapSort(sortforTFGS), *SA = new int[N], vs, Js, lvs = -2, lJs = -2;
    int lI = -1, SI = 0, lSI = -1, PN, lPN = -1, wv, nDig;
    double WN, err;
    QString IsoStr[nIso], B;
    QStringList SourceList;
    QTextStream S(&File);
    for (n=0; n<N; n++) SA[S1[n]] = n;
    delete[] S1;
    S << ("    " + QString::number(N)).right(5) << "   0" << ("    " + QString::number(atom1->getnIso())).right(5)
      << ("    " + QString::number(atom2->getnIso())).right(5) << ("    " + QString::number(State->getOmega())).right(5) 
      << "    0    1     levels, Dunham parameters, isotope A and B, omega, select dublicates (=0 no, =1 best, = 2 average)\n";
    for (n=0; n < atom1->getnIso(); n++) 
        S << ("    " + QString::number(atom1->getnNuc(n))).right(5) 
          << ("            " + QString::number(atom1->getIsoMass(n), 'f', 8)).right(13) << "           0.0\n";
    for (n=0; n < atom2->getnIso(); n++) 
        S << ("    " + QString::number(atom2->getnNuc(n))).right(5) 
          << ("            " + QString::number(atom2->getIsoMass(n), 'f', 8)).right(13) << "           0.0\n";
    S << "    0" 
      << ("       " + (State->getDunTable() != 0 ? QString::number(State->getDunTable()->getwe(), 'f', 1) : QString("60.0"))).right(8)
      << "  -0.0          0" << ("    " + QString::number(getvMax())).right(5) << "    0"
      << ("    " + QString::number(getJMax())).right(5) + "      v Verschiebung, vibr. spacing, level limits\n"
      << "  0.0              1           constant zero, switch of quantum numbers\n";
    for (n=0; n < nIso; n++)
        IsoStr[n] = ("    " + QString::number(Iso->mNumIso1[n])).right(5) + ("    " + QString::number(Iso->mNumIso2[n])).right(5);
    delete Iso;
    for (n=0; n<N; n++)
    {
        WN = Tab->item(SA[n], 8)->text().toDouble();
        B = Tab->item(SA[n], 9)->text();
        nDig = B.length() - B.indexOf('.') - 1;
        err = B.toDouble();
        PN = Tab->item(SA[n], 6)->text().toInt();
        Js = Tab->item(SA[n], 4)->text().toInt();
        I = Tab->item(SA[n], 0)->text().toInt();
        B = Tab->item(SA[n], 3)->text();
        vs = (B != "nA" ? B.toInt() : -1);
        if (B == "TE")
        {
            wv = 9999;
            WN *= -1.0;
            SI = -1;
        }
        else 
        {
            if (PN != lPN)
            {
                if (SI == lSI && lvs == vs && Js == lJs && lI == I) vo += 100;
                else vo = 0;
            }
            B = Tab->item(SA[n], 5)->text();
            SI = SourceList.indexOf(B);
            if (SI == -1)
            {
                SI = SourceList.count();
                SourceList << B;
            }
            wv = 1000 * (SI + 1) + vo + vs;
        }
        B = Tab->item(SA[n], 2)->text();
        S << IsoStr[I] << ("    " + QString::number(wv)).right(5) << ("    " + (wv != 9999 ? QString::number(Js) : B)).right(5)
          << ("     " + Tab->item(SA[n], 1)->text()).right(5)  << ("     " + B).right(5) + IsoStr[I] << "    0    0    0    0\n";
        S << ("              " + QString::number(WN, 'f', nDig)).right(15) 
          << ("              " + QString::number(err, 'f', nDig)).right(15) << "    2\n";
        lvs = vs;
        lJs = Js;
        lI = I;
        lSI = SI;
        lPN = PN;
    }
    delete[] SA;
    return true;
}

void FitData::setResidualFit(ResidualFit *i_residualFit)
{
    connect(i_residualFit, SIGNAL(Changed()), this, SLOT(Changed()));
    for (QList<ResidualFit*>::Iterator it = residualFits.begin(); it != residualFits.end(); ++it)
        if ((*it)->getIso() == i_residualFit->getIso() && (*it)->getv() == i_residualFit->getv() && (*it)->getComp() == i_residualFit->getComp()
                && (*it)->getStateName() == i_residualFit->getStateName())
    {
        disconnect(*it, SLOT(Changed()), this, SLOT(Changed()));
        delete *it;
        *it = i_residualFit;
        return;
    }
    residualFits.push_back(i_residualFit);
    Changed();
}

ResidualFit* FitData::getResidualFit(ElState * const i_state, const int i_Iso, const int i_v, const int i_comp)
{
    for (QList<ResidualFit*>::iterator it = residualFits.begin(); it != residualFits.end(); ++it)
        if ((*it)->getIso() == i_Iso && (*it)->getv() == i_v && (*it)->getComp() == i_comp && (*it)->getStateName() == i_state->getName()) return *it;
    return 0;
}

int *FitData::heapSort(bool sortFuncs(const FitDataCore *const, const int, const int)) const
{
    return utils::heapSort(FitDataCoreSortFunctor(Tab, sortFuncs), getNumLines());
}
