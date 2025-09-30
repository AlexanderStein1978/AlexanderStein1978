//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//

#include "spectlist.h"
#include "intprog.h"
#include "Spektrum.h"
#include "MainWindow.h"

#include <QFile>
#include <QFileDialog>
#include <QTextStream>


SpectList::SpectList(MainWindow *MW) : TableWindow(MDIChild::TextTable1, MW)
{
	setWindowTitle("List of spectra");
	Tab->setColumnCount(2);
	Tab->setHorizontalHeaderLabels(QStringList() << "File name" << "Laser frequency");
	setFilter("SpectList (*.slist)");
	setFileExt(".slist");
}

void SpectList::AutoSLP()
{
	if (MW == 0) return;
	int n, m, N = Tab->rowCount(), NL;
	QFile F, Log("AutoSLP.log");
	QStringList Dirs;
	QString FN;
	IntProg **P;
	for (n=0; n<N; n++)
	{
		if (Tab->item(n, 0) == 0) continue;
		if ((FN = Tab->item(n, 0)->text()).isEmpty()) continue;
		F.setFileName(FN);
		if (!F.exists())
		{
			FN = FN.right(FN.length() - FN.lastIndexOf(QRegExp("[\\/]"))  - 1);
			for (m=0; m < Dirs.count(); m++)
			{
				F.setFileName(Dirs[m] + FN);
				if (F.exists()) break;
			}
			if (m < Dirs.count()) Tab->item(n, 0)->setText(Dirs[m] + FN);
			else
			{
				FN = QFileDialog::getOpenFileName(this, "Show file " + FN, FN, "Scpectra (*.spect *.dat)");
				if (!FN.isEmpty())
				{
					Tab->item(n, 0)->setText(FN);
					FN = FN.left(FN.lastIndexOf(QRegExp("[\\/]")) + 1);
					for (m=0; m < Dirs.count(); m++) if (Dirs[m] == FN) break;
					if (m == Dirs.count()) Dirs << FN;
				}
			}
		}
	}
	Log.open(QIODevice::WriteOnly);
	QTextStream LS(&Log);
	Spektrum *S = MW->CreateSpectrum();
	for (n=0; n<N; n++)
	{
		if (Tab->item(n, 0) == 0) continue;
		if ((FN = Tab->item(n, 0)->text()).isEmpty()) continue;
		if (!S->readData(FN))
		{
			printf("Error reading file %s\n", FN.toLatin1().data());
			LS << "Error reading file " + FN + "\n";
			continue;
		}
		printf("File %s:\n", FN.toLatin1().data());
		LS << "File " + FN + ":\n";
        if (S->getNPoints() == 0)
		{
			printf("Error, file is not a spectrum!\n");
			LS << "Error, file is not a spectrum!\n";
			continue;
		}
		S->ShowMarker();
		S->MSLP(P, NL);
		printf("%d new progressions found", NL);
		LS << QString::number(NL) + " new progressions found\n";
		for (m=0; m < NL; m++)
		{
			S->ClearMarked();
			S->MarkProgression(*P[m]);
			S->addMarker(true);
			delete P[m];
		}
	}
}
