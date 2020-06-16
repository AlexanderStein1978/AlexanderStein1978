//
// C++ Implementation: ThreadControl
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2006 - 2020
//
// Copyright: See README file that comes with this source code
//
//

#include "ThreadControl.h"
#include "MainWindow.h"
#include "tablewindow.h"
#include "elstate.h"
#include "potential.h"

#include <QGridLayout>
#include <QLineEdit>
#include <QLabel>
#include <QProgressBar>
#include <QIntValidator>
#include <QMessageBox>
#include <QCloseEvent>


ThreadControl::ThreadControl(MainWindow* mw, ElState* nSt, int NumIterations, int NumParFits)
{
    int n;
    MW = mw;
    St = nSt;
    MFD = 0;
    NFC = 0;
    AdTabRows = 0;
    NumIt = NumIterations;
    MMParF = MaxParF = NumParFits;
    NumParF = ItFinished = NumTermF = 0;
    NumResultTabColumns = 0;
    PotI = new int[MaxParF];
    Pots = new Potential*[MaxParF];
    for (n=0; n < MaxParF; n++)
    {
        PotI[n] = -1;
        Pots[n] = 0;
    }
    ItI = new int[NumIterations];
    for (n=0; n < NumIterations; n++) ItI[n] = -1;
    QGridLayout *L = new QGridLayout(this);
    L->addWidget(new QLabel("Number of iterations:", this), 0, 0);
    L->addWidget(NumItE = new QLineEdit(QString::number(NumIt), this), 0, 1);
    L->addWidget(new QLabel("Number of parallel fits:", this), 1, 0);
    L->addWidget(NumParFE = new QLineEdit(QString::number(MaxParF), this), 1, 1);
    L->setRowMinimumHeight(2, 20);
    L->addWidget(Progress = new QProgressBar(this), 3, 0, 1, 2);
    Progress->setTextVisible(true);
    Progress->setRange(0, NumIt);
    Progress->setValue(0);
    Progress->setFormat("%v iterations finished");
    NumItE->setValidator(NumItEValid = new QIntValidator(0, 100000, NumItE));
    NumParFE->setValidator(NumParFEValid = new QIntValidator(0, 1000, NumParFE));
    connect(NumItE, SIGNAL(editingFinished()), this, SLOT(NItChanged()));
    connect(NumParFE, SIGNAL(editingFinished()), this, SLOT(NParFChanged()));
    ResultTab = MW->createTableWindow();
    ResultTab->show();
}

ThreadControl::~ThreadControl()
{
    delete[] Pots;
    delete[] PotI;
    delete[] ItI;
    if (NFC > 0) delete[] FCRows;
}

void ThreadControl::CalcFinished(int ThreadNum)
{
    PotI[ThreadNum] = -1;
    NumParF--;
    NumTermF++;
    Progress->setValue(++ItFinished);
    if (ItFinished + NumParF < NumIt && NumParF < MaxParF) StartCalcs();
    else if (ItFinished == NumIt) AllCalculationsFinished();
}

void ThreadControl::AllCalculationsFinished()
{
    QMessageBox::information(this, "MolSpektAnalysis", "The calculations are finished.");
}

void ThreadControl::CalcTerminated(int ThreadNum)
{
    MaxParF--;
    NumParF--;
    NumTermF++;
    ItI[PotI[ThreadNum]] = -1;
    PotI[ThreadNum] = -1;
    NumParFE->setText(QString::number(MaxParF));
}

void ThreadControl::closeEvent(QCloseEvent* event)
{
    if (NumParF > 0 ? QMessageBox::information(this, "MolSpektAnalysis", 
            "Closing this window will make you unable to directly control the running calculations", 
            QMessageBox::Ok | QMessageBox::Cancel, QMessageBox::Cancel) == QMessageBox::Cancel : false) 
        event->ignore();
    else 
    {
        QString Name;
        int i, j, N;
        for (i=0; i < MMParF; i++)
        {
            if (Pots[i] == 0) continue;
            for (j=0, N = St->getNumPotentials(), Name = Pots[i]->getName(); j<N; j++)
                if (St->getPotentialName(j) == Name) St->removePotential(j);
        }
        event->accept();
    }
}

void ThreadControl::NItChanged()
{
    int n, nIt = NumItE->text().toInt(), *nItI;
    ResultTab->setTabDimensions(nIt + AdTabRows, NumResultTabColumns);
    if (nIt > NumIt)
    {
        nItI = new int[nIt];
        for (n=0; n < NumIt; n++) nItI[n] = ItI[n];
        while (n < nIt) nItI[n++] = -1;
        delete[] ItI;
        ItI = nItI;
    }
    NumIt = nIt;
    Progress->setMaximum(NumIt);
    if (NumParF < MaxParF && ItFinished + NumParF < NumIt) StartCalcs(); 
}

void ThreadControl::NParFChanged()
{
    MaxParF = NumParFE->text().toInt();
    if (MaxParF > MMParF)
    {
        Potential** nPots = new Potential*[MaxParF];
        int n, *nPotI = new int[MaxParF];
        for (n=0; n < MMParF; n++)
        {
            nPots[n] = Pots[n];
            nPotI[n] = PotI[n];
        }
        while (n < MaxParF)
        {
            nPots[n] = 0;
            nPotI[n++] = -1;
        }
        delete[] Pots;
        delete[] PotI;
        Pots = nPots;
        PotI = nPotI;
        MMParF = MaxParF;
    }
    if (NumParF < MaxParF && ItFinished + NumParF < NumIt) StartCalcs(); 
}

void ThreadControl::StartCalcs()
{
    int n, p;
    while (NumParF < MaxParF && ItFinished + NumParF < NumIt)
    {
        for (p=0; PotI[p] != -1; p++) ;
        for (n=0; ItI[n] != -1; n++) ;
        ItI[n] = p;
        PotI[p] = n;
        StartCalc(p, n);
        if (Pots[p] == 0) 
        {
            MaxParF = p;
            NumParFE->setText(QString::number(p));
            break;
        }
        NumParF++;
    }
    NumItEValid->setBottom(ItFinished + NumParF);
}

void ThreadControl::StartCalc(int, int)
{

}
