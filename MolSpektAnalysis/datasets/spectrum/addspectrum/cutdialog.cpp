//
// C++ Implementation: cutdialog
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2009 - 2016
//
// Copyright: See README file that comes with this source code
//
//

#include "cutdialog.h"
#include "MainWindow.h"
#include "Spektrum.h"

#include <stdio.h>

#include <QComboBox>
#include <QLabel>
#include <QGridLayout>
#include <QTableWidget>
#include <QPushButton>
#include <QRadioButton>
#include <QStringList>
#include <QMessageBox>

CutDialog::CutDialog(MainWindow *mw, Spektrum *S)
{
	setAttribute(Qt::WA_DeleteOnClose);
	if (mw == 0)
	{
		printf("Error: mw == 0, closing now!");
		close();
		return;
	}
	int n, N = mw->getNumSpectra();
	Spektrum *Sp;
	MW = mw;
	setWindowTitle("Cut spectrum");
	QGridLayout *Layout = new QGridLayout(this);
	QLabel *L1 = new QLabel("Spektrum to cut:", this);
	Layout->addWidget(L1, 0, 0);
	SpektrumBox = new QComboBox(this);
	SpektrumBox->setEditable(false);
	Layout->addWidget(SpektrumBox, 0, 1);
	for (n=0; n<N; n++) 
	{
		SpektrumBox->addItem((Sp = MW->getSpectrum(n))->getName());
		//printf("n=%d, N=%d, Sp=%s\n", n, N, Sp->getName().ascii());
		if (Sp == S) SpektrumBox->setCurrentIndex(n);
	}
	connect(SpektrumBox, SIGNAL(currentIndexChanged(QString)), this, SLOT(SpektrumChanged(QString)));
	connect(MW, SIGNAL(SpektrumDestroyed(int)), this, SLOT(SpektrumDeleted(int)));
	out = new QRadioButton("use outer areas", this);
	out->setChecked(true);
	Layout->addWidget(out, 1, 0);
	in = new QRadioButton("use inner areas", this);
	Layout->addWidget(in, 1, 1);
	Layout->setRowMinimumHeight(2, 10);
	QLabel *L2 = new QLabel("Ranges:", this);
	Layout->addWidget(L2, 3, 0);
	Tab = new QTableWidget(100, 2, this);
	Tab->setHorizontalHeaderLabels(QStringList() << "Left border" << "Right border");
	Layout->addWidget(Tab, 4, 0, 1, 2);
	Layout->setRowMinimumHeight(5, 20);
	QPushButton *apply = new QPushButton("Apply", this);
	Layout->addWidget(apply, 6, 0);
	connect(apply, SIGNAL(pressed()), this, SLOT(Apply()));
	QPushButton *Close = new QPushButton("Close", this);
	Layout->addWidget(Close, 6, 1);
	connect(Close, SIGNAL(pressed()), this, SLOT(close()));
	SpektrumChanged(SpektrumBox->currentText());
}

CutDialog::~CutDialog()
{
	disconnect(0, 0, this, SLOT(setValue(double, double)));
}

void CutDialog::SpektrumChanged(QString Name)
{
	disconnect(0, 0, this, SLOT(setValue(double, double)));
	int n, N = MW->getNumSpectra();
	for (n=0; (n<N ? (Spekt = MW->getSpectrum(n))->getName() != Name : false); n++) ;
	if (n<N) connect(Spekt, SIGNAL(clicked(double, double)), this, SLOT(setValue(double, double)));
}

void CutDialog::SpektrumDeleted(int Index)
{
	int n = SpektrumBox->currentIndex();
	SpektrumBox->removeItem(Index);
	if (n > Index && n != SpektrumBox->currentIndex() - 1) SpektrumBox->setCurrentIndex(n-1);
}

void CutDialog::setValue(double x, double /*y*/)
{
	int r = Tab->currentRow(), c = Tab->currentColumn();
	//printf("CutDialog::setValue, r=%d, c=%d\n", r, c);
	if (r==-1 || c==-1) Tab->setCurrentCell(r=0, c=0);
	if (Tab->currentItem() == 0) Tab->setItem(r, c, new QTableWidgetItem(QString::number(x, 'g', 9)));
	else Tab->currentItem()->setText(QString::number(x, 'g', 9));
	if (c==0) c++;
	else
	{
		c=0;
		r++;
	}
	if (r < 100) Tab->setCurrentCell(r, c);
}

void CutDialog::Apply()
{
	if (Spekt == 0)
	{
		QMessageBox::information(this, "MolSpektAnalysis", 
			"Please select an existing spectrum first!", QMessageBox::Ok);
		return;
	}
	int n, N;
	QTableWidgetItem *IB;
	for (N=0; (N < 100 ? Tab->item(N, 0) != 0 && Tab->item(N, 1) != 0 : false); N++) ;
	for (n=0; n<N; n++) if (Tab->item(n, 0)->text().toDouble() > Tab->item(n, 1)->text().toDouble())
	{
		IB = Tab->takeItem(n, 0);
		Tab->setItem(n, 0, Tab->takeItem(n, 1));
		Tab->setItem(n, 1, IB);
	}
	for (IB = Tab->item(0, 0); IB != 0; IB = 0) for (n=1; n<N; n++)
			if (Tab->item(n, 0)->text().toDouble() < Tab->item(n-1, 0)->text().toDouble())
	{
		IB = Tab->takeItem(n-1, 0);
		Tab->setItem(n-1, 0, Tab->takeItem(n, 0));
		Tab->setItem(n, 0, IB);
		IB = Tab->takeItem(n-1, 1);
		Tab->setItem(n-1, 1, Tab->takeItem(n, 1));
		Tab->setItem(n, 1, IB);
	}
	if (in->isChecked()) N++;
	double Start[N], Stop[N];
	if (in->isChecked())
	{
		Start[0] = 0.0;
		for (n=0; n<N-1; n++) 
		{
			Stop[n] = Tab->item(n, 0)->text().toDouble();
			Start[n+1] = Tab->item(n, 1)->text().toDouble();
		}
		Stop[n] = Spekt->getMaxR() + 1.0;
	}
	else for (n=0; n<N; n++)
	{
		Start[n] = Tab->item(n, 0)->text().toDouble();
		Stop[n] = Tab->item(n, 1)->text().toDouble();
	}
	Spekt->cut(Start, Stop, N);
}
