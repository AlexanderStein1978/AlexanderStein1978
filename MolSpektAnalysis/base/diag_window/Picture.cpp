//
// C++ Implementation: Picture
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2006 - 2016
//
// Copyright: See README file that comes with this source code
//
//

#include "Picture.h"
#include "math.h"

#include <qpainter.h>
#include <QPixmap>
#include <QResizeEvent>
#include <QMouseEvent>
#include <QKeyEvent>
#include <QPaintEvent>

Picture::Picture(QWidget *Parent, const char */*name*/) : QWidget(Parent)
{
	Pixmap = 0;
    QSize S = size();
    Resize(S);
    fdrawn = false;
    mPressed = false;
	Zoom = true;
	splitScreen = false;
	setMouseTracking(true);
	crossLeft = crossTop = 0;
	crossRight = width();
	crossBottom = height();
	mouseCross = false;	
	crossX = crossY = -1;
	rmx = rmy = 0;
}

Picture::~Picture() 
{
	if (Pixmap != 0) delete Pixmap;
}

void Picture::setSplitScreen(const bool E)
{
	splitScreen = E;
	if (E) Zoom = false;
}

void Picture::setZoom(const bool E)
{
	Zoom = E;
	if (E) splitScreen = false;
}

bool Picture::splitScreenEnabled() const
{
	return splitScreen;
}

bool Picture::ZoomEnabled() const
{
	return Zoom;
}

void Picture::GetMW(const int &x1, const int &x2, int &m, int &w)
{
    if (x1 > x2)
    {
	m = x2;
	w = x1 - x2;
    }
    else
    {
	m = x1;
	w = x2 - x1;
    }
    return;
}

void Picture::Line(const int &x1, const int &y1, const int &x2, const int &y2, QColor *Color)
{
    //printf("Beginn von Line\n");
    mPressed = false;
    mx1 = mx2 = my1 = my2 = 0;
    QPainter P(Pixmap);
    if (Color != 0) P.setPen(*Color);
    P.drawLine(x1, y1, x2, y2);
    update(x1, y1, x2, y2);
    //printf("Ende von Line\n");
}

void Picture::clear()
{
    QSize S = size();
    int w = S.width();
    int h = S.height();
    QRect A(1, 1, w - 2, h - 2);
    QPainter P(Pixmap);
    P.setBackground(QColor(255, 255, 255));
    P.eraseRect(A);
    update(A);
}

void Picture::drawVertText(const int &x, const int &y, const QString &T)
{
    //printf("DrawVertText\n");
    int X = -1 * y;
    int Y = x;
    QPainter P(Pixmap);
    P.rotate(270);
    P.drawText(X, Y, T);
    update();
}

QRect Picture::contentsRect()
{
    QSize S = size();
    int w = S.width();
    int h = S.height();
    QRect A(1, 1, w - 2, h - 2);
    return A;
}

void Picture::paintEvent(QPaintEvent *P) 
{
    //printf("PaintEvent\n");
    QRect R = P->rect();
    QPoint D(R.left(), R.top());
    QPainter Paint(this);
    if (Pixmap!= 0) Paint.drawPixmap(D, *Pixmap, R);
	if (mouseCross && crossX > crossLeft && crossX < crossRight && crossY > crossTop 
		   && crossY < crossBottom) 
	{
		Paint.setPen(QColor(170, 170, 170));
		Paint.drawLine(crossX, crossTop, crossX, crossBottom);
		Paint.drawLine(crossLeft, crossY, crossRight, crossY);
		Paint.setPen(QColor(0, 0, 0));
	}
    if (mPressed && Zoom) Paint.drawRect(SRect);
	if (splitScreen) 
	{
		Paint.setPen(QColor(0, 0, 255));
		Paint.drawLine(mx1, my1, mx2, my2);
	}
    fdrawn = true;
}

void Picture::resizeEvent(QResizeEvent *R)
{
    QSize S = R->size();
    Resize(S);
    emit Resized(&S);
}

void Picture::mouseOutside()
{
	if (crossX > 0 || crossY > 0)
	{
		crossX = crossY = -1;
		if (mouseCross) update();
	}
}

void Picture::mouseMoveEvent(QMouseEvent *e) 
{
    //printf("MouseMoveEvent\n");
    e->setAccepted(false);
	emit mouseMoved(e);
	if (mouseCross) 
	{
		crossX = e->x();
		crossY = e->y();
		update();
	}
	if (mPressed && !e->isAccepted())
    {
		int xm, ym, xw, yw;
		GetMW(mx1, mx2, xm, xw);
		GetMW(my1, my2, ym, yw);
		QRect RA(xm, ym, xw + 1, yw + 1);
		mx2 = e->x();
		my2 = e->y();
		GetMW(mx1, mx2, xm, xw);
		GetMW(my1, my2, ym, yw);
		SRect.setRect(xm, ym, xw, yw);
		update(MaxRect(RA, SRect));
	}
	emit MouseMoved(e->x(), e->y());
	e->setAccepted(true);
}

void Picture::mousePressEvent(QMouseEvent *e)
{
    //printf("MousePressEvent\n");
	e->setAccepted(false);
	emit mousePressed(e);
	if (e->isAccepted()) return;
    if (e->button() == Qt::LeftButton)
    {
		mx2 = mx1 = e->x();
		my2 = my1 = e->y();
		if (e->modifiers() == Qt::ShiftModifier)
		{
			mx1 = rmx;
			my1 = rmy;
		}
		SRect.setRect(mx1, my1, 1, 1);
		mPressed = true;
		e->setAccepted(true);
    }
	if (splitScreen) 
	{
		update();
		e->setAccepted(true);
	}
}

void Picture::mouseReleaseEvent(QMouseEvent *e)
{
    //printf("MouseReleaseEvent\n");
    e->setAccepted(false);
	emit mouseReleased(e);
	if (!e->isAccepted())
	{
		if (mPressed && Zoom && mx1 != mx2 && my1 != my2)
    	{
			int xm, ym, xw, yw;
			rmx = mx2;
			rmy = my2;
			GetMW(mx1, mx2, xm, xw);
			GetMW(my1, my2, ym, yw);
			QRect RA(xm, ym, xw + 1, yw + 1);
			QRect R(xm, ym, xw, yw);
			if (xw > 10 && yw > 10)
			{
				emit SelectionChanged(&R);
				emit SelectionChanged(&R, e->modifiers() == Qt::ControlModifier);
				update(RA);
			}
    	}
		if (mPressed && splitScreen && (mx1 != mx2 || my1 != my2))
		{
			//printf("mx1=%d, my1=%d, mx2=%d, my2=%d\n", mx1, my1, mx2, my2);
			if (mx1 == mx2)
			{
				my1 = 0;
				my2 = height();
			}
			else 
			{
				if (mx1 > mx2)
				{
					int buff = mx1;
					mx1 = mx2;
					mx2 = buff;
					buff = my1;
					my1 = my2;
					my2 = buff;
				}
				double St = (double)(my2 - my1) / (mx2 - mx1);
				//printf("St=%f\n", St);
				if (St <= 1.0)
				{
					my1 = my1 - (int)(St * (double)mx1);
					mx1 = 0;
					my2 = my2 + (int)(St * (double)(width() - mx2));
					mx2 = width();
				}
				else
				{
					mx1 = mx1 - (int)((double)my1 / St);
					my1 = 0;
					mx2 = mx2 + (int)((double)(height() - my2) / St);
					my2 = height();
				}
			}
			//printf("mx1=%d, my1=%d, mx2=%d, my2=%d\n", mx1, my1, mx2, my2);
			update();
			emit ScreenSplitted(mx1, my1, mx2, my2);
		}
		mPressed = false;
		if (fabs(mx1 - mx2) <= 10.0 && fabs(my1 - my2) <= 10.0)
		{
			if (e->button() == Qt::LeftButton) 
			{
				QPoint P(e->x(), e->y());
				emit LeftClicked(&P);
				emit LeftClicked(&P, e->modifiers() == Qt::ControlModifier);
			}
    		else if (e->button() == Qt::RightButton) 
			{
				QPoint P(e->globalPos());
				emit RightClicked();
				emit RightClicked(&P);
			}
		}
    	if (!splitScreen) mx1 = mx2 = my1 = my2 = 0;
	}
	QWidget::mouseReleaseEvent(e);
}

void Picture::Resize(QSize &Size)
{
    int w = Size.width(), h = Size.height();
	if (Pixmap != 0) 
	{
		delete Pixmap;
		Pixmap = new QPixmap(w, h);
	}
	/*
	QPainter P(&Pixmap);
    P.setBackground(QColor(255, 255, 255));
    P.eraseRect(0, 0, --w, --h); 
    P.drawLine(0, 0, w, 0);
    P.drawLine(0, 0, 0, h);
    P.drawLine(w, 0, w, h);
    P.drawLine(0, h, w, h);
    update();*/
}

void Picture::keyPressEvent(QKeyEvent *K)
{
    emit KeyPressed(K);
}

QPixmap *Picture::getPixmap()
{
    if (Pixmap == 0) Pixmap = new QPixmap(width(), height());
	return Pixmap;
}

void Picture::UpdatePicture()
{
    if (fdrawn) update();
}

void Picture::GetMW(const int l1, const int l2, const int r1, const int r2, 
					int &l, int &w)
{
	int r;
	if (l1 < l2) l = l1 - 1;
	else l = l2 - 1;
	if (r1 > r2) r = r1 + 1;
	else r = r2 + 1;
	w = r - l + 1;
}

QRect Picture::MaxRect(const QRect &Rect1, const QRect &Rect2)
{
	int ML, MT, MW, MH;
	GetMW(Rect1.left(), Rect2.left(), Rect1.right(), Rect2.right(), ML, MW);
	GetMW(Rect1.top(), Rect2.top(), Rect1.bottom(), Rect2.bottom(), MT, MH);
	return QRect(ML, MT, MW, MH);
}

bool Picture::getMouseCross() const
{
	return mouseCross;
}

void Picture::setMouseCross(const bool &show)
{
	mouseCross = show;
}

void Picture::getMouseCrossBorders(int &left, int &top, int &right, int &bottom) const
{
	left = crossLeft;
	top = crossTop;
	right = crossRight;
	bottom = crossBottom;
}

void Picture::setMouseCrossBorders(const int &left, const int &top, const int &right, const int &bottom)
{
	crossLeft = left;
	crossTop = top;
	crossRight = right;
	crossBottom = bottom;
}

