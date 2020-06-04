//
// C++ Interface: SHOWACTION
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2006 - 2019
//
// Copyright: See README file that comes with this source code
//
//


#ifndef SHOWACTION_H
#define SHOWACTION_H


#include <QAction>


class ShowAction : public QAction
{
	Q_OBJECT
			
	public:
		ShowAction(const QString &text, QObject *parent) : QAction(text, parent)
		{
		}
		
	public slots:
	
		inline void setName(QString name)
		{
			setText(name);
		}
		
};

#endif
