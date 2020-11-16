#include "particlewatchtable.h"
#include "window.h"
#include "watchpoint.h"

#include <QGridLayout>


ParticleWatchTable::ParticleWatchTable(Window *iwindow, MainWindow *MW) : TableWindow(External, MW), mWatchPoint(iwindow->getNumSteps(), iwindow->getNumParticles()),
    CoordBox(new QComboBox(this)), StepBox(new QComboBox(this)), ParticleBox(new QComboBox(this)), SumLabel(new QLabel(this)), NextButton(new QPushButton("Next iteration", this)),
    window(iwindow)
{
    QGridLayout *L = new QGridLayout(this);
    int numParticles(window->getNumParticles()), xDim(window->getXDim()), yDim(numParticles / xDim);

    window->setParticleWatchPoint(&mWatchPoint);

    L->addWidget(new QLabel("Particle:", this), 0, 0);
    L->addWidget(ParticleBox, 0, 1);
    L->addWidget(SumLabel, 0, 2);
    L->addWidget(NextButton, 0, 3);
    L->addWidget(new QLabel("Step:", this), 0, 4);
    L->addWidget(StepBox, 0, 5);
    L->addWidget(new QLabel("Coordinate:", this), 0, 6);
    L->addWidget(CoordBox, 0, 7);
    L->addWidget(Tab = new QTableWidget(2 * yDim, xDim, this), 1, 0, 1, 8);

    for (int i=0; i < xDim; ++i)
    {
        Tab->setHorizontalHeaderItem(i, new QTableWidgetItem(QString::number(2*i)));
        for (int j=0; j < 2 * yDim; ++j) Tab->setItem(j, i, new QTableWidgetItem);
    }
    for (int i=0; i < yDim; ++i)
    {
        Tab->setVerticalHeaderItem(2*i, new QTableWidgetItem(QString("%1 - %2").arg(2*i * xDim).arg(2*(i+1) * xDim - 2)));
        Tab->setVerticalHeaderItem(2*i+1, new QTableWidgetItem(QString("%1 - %2").arg(2*i * xDim + 1).arg(2*(i+1) * xDim - 1)));
    }

    ParticleBox->setValidator(new QIntValidator(0, numParticles - 1, ParticleBox));
    for (int i=0; i < numParticles; ++i) ParticleBox->addItem(QString::number(i));

    StepBox->setEditable(false);
    for (int i=0; i < window->getNumSteps(); ++i) StepBox->addItem(QString::number(i));

    CoordBox->setEditable(false);
    CoordBox->addItems(QStringList() << "x" << "y" << "z");

    connect(NextButton, SIGNAL(clicked()), this, SLOT(nextClicked()));
    connect(StepBox, SIGNAL(currentIndexChanged(int)), this, SLOT(stepChanged(int)));
    connect(CoordBox, SIGNAL(currentIndexChanged(int)), this, SLOT(coordChanged(int)));
}

void ParticleWatchTable::coordChanged(const int i)
{
    fillTable(StepBox->currentIndex(), i);
}

void ParticleWatchTable::stepChanged(const int i)
{
    SumLabel->setText(QString("Sum a: (%1;%2;%3)").arg(mWatchPoint.getSumX(i)).arg(mWatchPoint.getSumY(i)).arg(mWatchPoint.getSumZ(i)));
    fillTable(i, CoordBox->currentIndex());
}

void ParticleWatchTable::nextClicked()
{
    mWatchPoint.reset();
    window->setParticleWatch(ParticleBox->currentText().toInt());
    stepChanged(StepBox->currentIndex());
}

void ParticleWatchTable::fillTable(const int step, const int coord)
{
    int xDim(window->getXDim()), yDim(window->getNumParticles() / xDim);
    for (int row = 0; row < yDim; ++row) for (int column = 0; column < xDim; ++column)
    {
        Tab->item(2 * row, column)->setText(QString::number(mWatchPoint.get(step, 2 * (row * xDim + column), coord)));
        Tab->item(2 * row + 1, column)->setText(QString::number(mWatchPoint.get(step, 2 * (row * xDim + column) + 1, coord)));
    }
}
