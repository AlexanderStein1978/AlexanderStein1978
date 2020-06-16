//
// C++ Implementation: PotFitLogWindow
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2006 - 2019
//
// Copyright: See README file that comes with this source code
//
//


#include "PotFitLogWindow.h"

#include <QGridLayout>
#include <QPushButton>
#include <QPlainTextEdit>
#include <QMessageBox>


PotFitLogWindow::PotFitLogWindow(QWidget*) : MDIChild()
{
    QGridLayout *L = new QGridLayout(this);
    fitFinished = fitRunning = fitStopped = false;
    Stop = new QPushButton("Stop fit", this);
    Cancel = new QPushButton("Cancel fit", this);
    Close = new QPushButton("Close window", this);
    Stop->setEnabled(false);
    Cancel->setEnabled(false);
    Close->setEnabled(false);
    L->addWidget(Stop, 0, 0);
    L->addWidget(Cancel, 0, 1);
    L->addWidget(Close, 0, 2);
    L->addWidget(Log = new QPlainTextEdit(this), 1, 0, 1, 3);
    Log->setReadOnly(true);
    L->setRowStretch(1, 1);
    connect(Stop, SIGNAL(clicked()), this, SLOT(stopFit()));
    connect(Cancel, SIGNAL(clicked()), this, SLOT(cancelFit()));
    connect(Close, SIGNAL(clicked()), this, SLOT(close()));
}

void PotFitLogWindow::addTextRow(QString Text)
{
    Log->appendPlainText(Text);
}

void PotFitLogWindow::cancelFit()
{
    if (QMessageBox::question(this, "MolSpektAnalysis", "Do you really want to kill the running fit?", 
                              QMessageBox::Yes | QMessageBox::No, QMessageBox::No) == QMessageBox::Yes)
    {
        emit CancelFit(false);
        Stop->setEnabled(false);
        Cancel->setEnabled(false);
    }
}

bool PotFitLogWindow::askForQuit()
{
    if (!fitFinished && !fitRunning && !fitStopped)
    {
        QMessageBox::information(this, "MolSpektAnalysis", "The window cannot be closed while a fit is starting!");
        return false;
    }
    if (!fitRunning) return true;
    if (QMessageBox::information(this, "MolSpektAnalysis", "Closing this window kills the running fit!",
                                 QMessageBox::Ok | QMessageBox::Cancel, QMessageBox::Cancel) == QMessageBox::Ok)
    {
        emit CancelFit(true);
        Stop->setEnabled(false);
        Cancel->setEnabled(false);
        Close->setEnabled(false);
        Stop->setText("Stop fit");
        return true;
    }
    else return false;
}

void PotFitLogWindow::FitFinished()
{
    fitRunning = false;
    fitFinished = true;
    Stop->setText("Start fit");
    Cancel->setEnabled(false);
    Log->appendPlainText("Fit finished");
}

void PotFitLogWindow::FitStopped()
{
    fitRunning = false;
    fitStopped = true;
    Cancel->setEnabled(false);
    Stop->setEnabled(true);
    Log->appendPlainText("Fit stopped");
}

void PotFitLogWindow::FitRunning()
{
    if (fitStopped)
    {
        Log->appendPlainText("Fit restarted");
        fitStopped = false;
    }
    fitFinished = false;
    fitRunning = true;
    Stop->setEnabled(true);
    Cancel->setEnabled(true);
    Close->setEnabled(true);
}

void PotFitLogWindow::FitTerminated()
{
    fitRunning = false;
    fitFinished = true;
    Stop->setText("Start fit");
    Cancel->setEnabled(false);
    Log->appendPlainText("Fit terminated!");
}

void PotFitLogWindow::stopFit()
{
    if (Stop->text() == "Stop fit")
    {
        emit StopFit();
        Stop->setText("Restart fit");
    }
    else if (Stop->text() == "Restart fit")
    {
        emit RestartFit();
        Stop->setEnabled(false);
    }
    else 
    {
        emit StartFit();
        fitFinished = fitRunning = fitStopped = false;
        Stop->setEnabled(false);
        Cancel->setEnabled(false);
        Close->setEnabled(false);
        Stop->setText("Stop fit");
    }
}
