#include "soundmainwindow.h"

#include <QMdiArea>


SoundMainWindow::SoundMainWindow()
{
    workspace = new QMdiArea;
	setCentralWidget(workspace);
}

void SoundMainWindow::showMDIChild(QWidget* C)
{
    workspace->addSubWindow(C);
	C->show();
}
