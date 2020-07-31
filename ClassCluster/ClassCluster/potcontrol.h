#ifndef POTCONTROL_H
#define POTCONTROL_H


#include <QObject>

class ControlWindow;
class Potential;
class PotentialPlot;
class PotStruct;

class QLineEdit;
class QPushButton;
class QCheckBox;
class QGridLayout;
class QTextStream;


class PotControl : QObject
{
    Q_OBJECT

public:
    PotControl(ControlWindow* parent);
    virtual ~PotControl();

    void Init(const QString& data);
    void Serialize(QTextStream& stream) const;
    void FillLayout(QGridLayout* layout, const int row) const;
    void FillStruct(PotStruct& potStruct) const;
    bool canPotBeClosed() const;

    inline bool isChanged() const
    {
        return changed;
    }

    inline void setPlot(PotentialPlot *const i_plot)
    {
        plot = i_plot;
    }

private slots:

    void Open();
    void ShowOpenDialog();
    void Save();
    void SaveAs();
    void Plot(const bool show);
    void RecalcExtensions();

    void Changed()
    {
        changed = true;
    }

private:
    void openPotential();

    ControlWindow* parent;
    Potential* pot;
    PotentialPlot* plot;
    QLineEdit *fileName, *VScale, *RScale;
    QPushButton *openB, *saveB, *saveAsB;
    QCheckBox *showBox;
    bool changed, changing;
};

#endif // POTCONTROL_H
