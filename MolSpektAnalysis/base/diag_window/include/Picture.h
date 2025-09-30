//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//

#ifndef PICTURE_H
#define PICTURE_H

#include <qwidget.h>
#include <qpixmap.h>
#include <QResizeEvent>
#include <QMouseEvent>
#include <QKeyEvent>
#include <QPaintEvent>

class Picture : public QWidget
{
    Q_OBJECT
    
public:
    Picture(QWidget *Parent = 0, const char *name = 0);
    ~Picture();
    void Line(const int &x1, const int &y1, const int &x2, const int &y2, QColor *Color = 0);
    void clear();
    void drawVertText(const int &x, const int &y, const QString &S);
    QRect contentsRect();
    QPixmap *getPixmap();
    void UpdatePicture();
	void setZoom(const bool E);
	bool ZoomEnabled() const;
	void setSplitScreen(const bool E);
	bool splitScreenEnabled() const;
	void setMouseCross(const bool &show);
	bool getMouseCross() const;
	void setMouseCrossBorders(const int &left, const int &top, const int &right, const int &bottom);
	void getMouseCrossBorders(int &left, int &top, int &right, int &bottom) const;
	void mouseOutside();
	
	inline bool savePicture(QString FileName, QString Format)
	{
		if (Pixmap != 0) return Pixmap->save(FileName, Format.toLatin1());
		return false;
	}
	
signals:
    void SelectionChanged(QRect *MarkedArea);
	void SelectionChanged(QRect *MarkedArea, bool ControlPressed);
    void Resized(QSize *Size);
    void RightClicked();
	void RightClicked(QPoint *GlobalPosition);
    void LeftClicked(QPoint *Position);
	void LeftClicked(QPoint *Position, bool ControlPressed);
    void KeyPressed(QKeyEvent *K);
	void ScreenSplitted(int x1, int y1, int x2, int y2);
	void MouseMoved(int x, int y);
	void mouseMoved(QMouseEvent *e);
	void mousePressed(QMouseEvent *e);
	void mouseReleased(QMouseEvent *e);
protected:
    void paintEvent(QPaintEvent *P);
    void resizeEvent(QResizeEvent *R);
    void mouseMoveEvent(QMouseEvent *e);
    void mousePressEvent(QMouseEvent *e);
    void mouseReleaseEvent(QMouseEvent *e);
    void keyPressEvent(QKeyEvent *K);
private:
    QPixmap *Pixmap;
    int mx1, my1, mx2, my2, crossLeft, crossTop, crossBottom, crossRight, crossX, crossY;
    int rmx, rmy;
	bool mPressed, Zoom, splitScreen;
    QRect SRect;
    bool fdrawn, mouseCross;
    QRect MaxRect(const QRect &Rect1, const QRect &Rect2);
    void GetMW(const int l1, const int l2, const int r1, const int r2, int &l, int &w);
    void Resize(QSize &Size);
    void GetMW(const int &x1, const int &x2, int &m, int &w);
};

#endif 
