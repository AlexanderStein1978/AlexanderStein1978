#ifndef LOGWINDOW_H
#define LOGWINDOW_H


#include "tablewindow.h"
#include "logmodel.h"

class QFile;
class QTextStream;
class QLineEdit;
class QPushButton;
class QCheckBox;
class QMutex;


class LogWindow : public TableWindow
{
    Q_OBJECT
public:
    explicit LogWindow(QWidget *parent = 0);

    void LogMessage(QStringList& message);
    void SetMessageBuffer(QMutex& loggerMutex, QList<QStringList> buffer);

private slots:
    void ShowFileDialog();
    void WriteLogFileChanged();

private:
    QLineEdit* mFilenameEdit;
    QPushButton* mFileDialogButton;
    QCheckBox* mWriteLogFileCheckBox;

    LogModel mModel;

    QMutex* mLoggerMutex;
    QTextStream *mLogStream;
    QFile *mLogFile;
    bool mWriteToFile;
};

#endif // LOGWINDOW_H
