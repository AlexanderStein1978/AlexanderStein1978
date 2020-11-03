#ifndef POTCONTROL_H
#define POTCONTROL_H


#include <QObject>

class ControlWindow;
class MainWindow;
class Potential;
class PotentialPlot;
class PotStruct;

class QLineEdit;
class QPushButton;
class QCheckBox;
class QGridLayout;
class QTextStream;
class QComboBox;


class PotControl : QObject
{
    Q_OBJECT

public:
    PotControl(ControlWindow* parent, MainWindow* mw);
    virtual ~PotControl();

    void Init(const QString& data);
    void Serialize(QTextStream& stream) const;
    void FillLayout(QGridLayout* layout, const int row) const;
    void FillStruct(PotStruct& potStruct) const;
    void UpdatePotentialBox();

    inline bool isChanged() const
    {
        return changed;
    }

private slots:
    void Plot(const bool show);
    void RecalcExtensions();
    void adjustRe();
    void PotentialBoxIndexChanged(const int newIndex);

    void Changed()
    {
        changed = true;
    }

private:
    void exchangePotential(Potential *const newPot);

    ControlWindow* parent;
    Potential* pot;
    QLineEdit *VScale, *RScale;
    QComboBox *PotentialBox;
    QPushButton *adjustReB;
    QCheckBox *showBox;
    MainWindow* MW;
    bool changed, changing;
};

#endif // POTCONTROL_H
