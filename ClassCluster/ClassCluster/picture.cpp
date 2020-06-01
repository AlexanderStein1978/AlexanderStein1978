#include "picture.h"

#include <QPixmap>
#include <QPainter>
#include <QRect>
#include <QPaintEvent>


Picture::Picture(QWidget* parent): QWidget(parent)
{
    Map = 0;
}

Picture::~Picture()
{
    if (Map != 0) delete Map;
}

QPixmap *Picture::getPixmap()
{
    if (Map == 0) Map = new QPixmap(width(), height());
    return Map;
}

void Picture::paintEvent(QPaintEvent *event)
{
    if (Map == 0) return;
    QPainter P(this);
    QRect R = event->rect();
    P.drawPixmap(R, *Map, R);
}
