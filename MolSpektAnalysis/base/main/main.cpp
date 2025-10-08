//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
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

