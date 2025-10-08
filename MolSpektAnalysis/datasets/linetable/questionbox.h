//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//

#ifndef QUESTIONBOX_H
#define QUESTIONBOX_H


#include <QMessageBox>


class QuestionBox : public QMessageBox
{
	Q_OBJECT
	
	public:
		QuestionBox(QString Caption, QString Question, QWidget* parent = 0);
		void addButton(QString Text);
		int getResult();
	
	private slots:
		void b0();
		void b1();
		void b2();
		void b3();
		void b4();
				
	private:
		int Result, NButtons;
};

#endif
