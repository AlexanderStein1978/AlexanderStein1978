#pragma once

#include <QMainWindow>

class QMdiArea;


class SoundMainWindow : public QMainWindow
{
public:
    SoundMainWindow();

    void showMDIChild(QWidget *C);

private:
    QMdiArea *workspace;
};
