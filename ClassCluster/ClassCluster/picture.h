#ifndef PICTURE_H
#define PICTURE_H


#include <QWidget>


class QPixmap;
class QPaintEvent;


class Picture : public QWidget
{
    public:
        Picture(QWidget *parent = 0);
        ~Picture();
        QPixmap *getPixmap();

    protected:
        void paintEvent(QPaintEvent *event);

    private:
        QPixmap *Map;
};

#endif // PICTURE_H
