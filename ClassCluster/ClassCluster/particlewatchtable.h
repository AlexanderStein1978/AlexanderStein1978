#ifndef PARTICLEWATCHTABLE_H
#define PARTICLEWATCHTABLE_H


#include "tablewindow.h"
#include "watchpoint.h"

class QComboBox;
class QPushButton;
class QLabel;

class Window;
class MainWindow;


class ParticleWatchTable : public TableWindow
{
    Q_OBJECT

public:

    ParticleWatchTable(Window* window, MainWindow* MW);

private slots:

    void stepChanged(const int i);
    void coordChanged(const int i);
    void nextClicked();

private:
    void fillTable(const int step, const int coord);

    WatchPoint mWatchPoint;

    QComboBox *CoordBox, *StepBox, *ParticleBox;
    QLabel *SumLabel;
    QPushButton *NextButton;
    Window* window;
};

#endif // PARTICLEWATCHTABLE_H
