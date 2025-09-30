//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//

#ifndef LINEWINDOWBASE_H
#define LINEWINDOWBASE_H


#include <QWidget>

#include <functional>

class QComboBox;
class QFocusEvent;

class MainWindow;
class Spektrum;
class LineProfile;


class LineWindowBase : public QWidget
{
    Q_OBJECT
public:
    LineWindowBase(MainWindow *mw, Spektrum* spect = nullptr, LineProfile *line = nullptr);

protected slots:
    void SpektrumChanged(const QString& SpectName);
    void LineChanged(const int index);
    void NumberOfLinesChanged();

protected:
    void focusInEvent(QFocusEvent *event) override;
    virtual void disconnectSpectrum() = 0;
    virtual void connectSpectrum() = 0;
    virtual void lineChanged() = 0;
    void modifyLine(const std::function<void()>& func);

    QComboBox *SpektrumBox, *LineBox;
    MainWindow *MW;
    Spektrum *mSpektrum;
    LineProfile* mLine;
};

#endif // LINEWINDOWBASE_H
