//
// C++ Implementation: main
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2006 - 2021
//
// Copyright: See README file that comes with this source code
//
//

#include <qapplication.h>

#include "MainWindow.h"

int main( int argc, char ** argv)
{
	QString LF;
	QApplication a( argc, argv );
	QFont F = a.font();
	F.setPixelSize(12);
	a.setFont(F);
	MainWindow *w = new MainWindow();
	/*while(n < argc)
	{
		LF = argv[n];
		if (LF.right(4) == ".log") w->setResumeComputationLogfile(LF);
	}*/
	w->show();
    a.connect( &a, SIGNAL( lastWindowClosed() ), &a, SLOT( quit() ) );
    int res = a.exec();
	delete w;
	return res;
}

