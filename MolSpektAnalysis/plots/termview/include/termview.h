//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//

#ifndef TERMVIEW_H
#define TERMVIEW_H

#include <tablewindow.h>

class TermTable;


class TermView : public TableWindow
{
	Q_OBJECT
public:
    TermView(TermTable *Term);
    ~TermView();
private slots:
	void ShowData();
private:
	double ****Data;
	int vMax, JMax, NComp, NIso;
	TermTable *TT;
};

#endif
