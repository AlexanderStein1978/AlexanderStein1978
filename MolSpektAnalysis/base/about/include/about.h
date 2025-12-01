//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//

#ifndef ABOUT_H
#define ABOUT_H

#include <QWidget>

class QPushButton;
class QTextBrowser;

class About : public QWidget
{
    Q_OBJECT

public:
    About();

signals:
    void closeThis();

private:
    QPushButton *close;
    QTextBrowser *browser;
};

#endif // ABOUT_H
