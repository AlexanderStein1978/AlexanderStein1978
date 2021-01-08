//
// C++ Interface: LineWindowBase
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2020 - 2020
//
// Copyright: See README file that comes with this source code
//
//

#ifndef LINEWINDOWBASE_H
#define LINEWINDOWBASE_H


#include <QWidget>

class QComboBox;
class QFocusEvent;

class MainWindow;
class Spektrum;
class Gaussian;


class LineWindowBase : public QWidget
{
    Q_OBJECT
public:
    LineWindowBase(MainWindow *mw, Spektrum* spect = nullptr, Gaussian* line = nullptr);

protected slots:
    void SpektrumChanged(const QString& SpectName);
    void LineChanged(const int index);
    void LineRemoved();

protected:
    void focusInEvent(QFocusEvent *event) override;
    virtual void disconnectSpectrum() = 0;
    virtual void connectSpectrum() = 0;
    virtual void lineChanged() = 0;

    QComboBox *SpektrumBox, *LineBox;
    MainWindow *MW;
    Spektrum *mSpektrum;
    Gaussian* mLine;

private:
    void ConnectSpectrum();
    void DisconnectSpectrum();
};

#endif // LINEWINDOWBASE_H
