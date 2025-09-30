//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//

#ifndef COMBOBOX_H
#define COMBOBOX_H 


#include <QComboBox>


class ComboBox : public QComboBox
{
	Q_OBJECT
	public:
		ComboBox();
		ComboBox(QWidget *parent);
		void setNumber(int Number);
		
	protected:
		virtual void mouseDoubleClickEvent(QMouseEvent *E);
		virtual void mousePressEvent(QMouseEvent *event);
		
	signals:
		void DoubleClicked(int Number);
		void RightClicked(QPoint P);
		
	private:
		int Number;
};

#endif
