#include "logwindow.h"

#include <QGridLayout>
#include <QCheckBox>
#include <QFileDialog>
#include <QLineEdit>
#include <QMessageBox>
#include <QPushButton>
#include <QTextStream>
#include <QMutex>


LogWindow::LogWindow(MainWindow *parent) : TableWindow(External, parent, nullptr), mFilenameEdit(new QLineEdit(this)), mMaxTableSize(new QLineEdit(this)), mFileDialogButton(new QPushButton("...", this)),
    mWriteLogFileCheckBox(new QCheckBox("Write messages to file", this)), mModel(), mLoggerMutex(nullptr), mLogStream(nullptr), mLogFile(nullptr), mWriteToFile(false)
{
    QGridLayout* L = new QGridLayout(this);
    L->addWidget(mWriteLogFileCheckBox, 0, 0);
    L->addWidget(new QLabel("Filename:", this), 0, 1);
    L->addWidget(mFilenameEdit, 0, 2);
    L->addWidget(mFileDialogButton, 0, 3);
    L->addWidget(new QLabel("Max table rows: ", this), 0, 4);
    L->addWidget(mMaxTableSize, 0, 5);
    L->setColumnMinimumWidth(3, 15);
    L->setColumnStretch(0, 10);
    L->setColumnStretch(1, 10);
    L->setColumnStretch(2, 10);
    L->addWidget(table = new MTable(this), 1, 0, 1, 4);
    L->setRowMinimumHeight(0, 20);
    L->setRowStretch(1, 10);
    mMaxTableSize->setValidator(new QIntValidator(0, 2147483647, mMaxTableSize));
    table->setModel(&mModel);
    connect(mWriteLogFileCheckBox, SIGNAL(stateChanged(int)), this, SLOT(WriteLogFileChanged(int)));
    connect(mFileDialogButton, SIGNAL(clicked()), this, SLOT(ShowFileDialog()));
    connect(mMaxTableSize, SIGNAL(editingFinished()), this, SLOT(MaxRowsCanged()));
}

LogWindow::~LogWindow()
{
    if (nullptr != mLogStream) delete mLogStream;
    if (nullptr != mLogFile) delete mLogFile;
}

void LogWindow::LogMessage(QStringList &message)
{
    mModel.AddLogMessage(message);
    if (mWriteToFile)
    {
        (*mLogStream) << message.join('\t') << '\n';
        mLogStream->flush();
    }
}

void LogWindow::MaxRowsCanged()
{
    mModel.SetMaxRows(mMaxTableSize->text().toInt());
}


void LogWindow::SetMessageBuffer(QMutex &loggerMutex, QList<QStringList> buffer)
{
    mModel.SetMessageBuffer(buffer);
    mLoggerMutex = &loggerMutex;
}

void LogWindow::ShowFileDialog()
{
    QString filename = QFileDialog::getSaveFileName(this, "Select filename for log file", "", "*.dat;;*.log;;*.txt;;*.csv");
    if (!filename.isEmpty()) mFilenameEdit->setText(filename);
}

void LogWindow::WriteLogFileChanged(int state)
{
    if (state == Qt::Unchecked)
    {
        mWriteToFile = false;
        mFilenameEdit->setReadOnly(false);
        mFileDialogButton->setEnabled(true);

    }
    else if (state == Qt::Checked)
    {
        QString filename = mFilenameEdit->text();
        if (filename.isEmpty())
        {
            QMessageBox::information(this, "ClassCluster logging", "Please choose a filename to write the log file.");
            mWriteLogFileCheckBox->setChecked(false);
            return;
        }
        if (nullptr == mLogFile) mLogFile = new QFile(filename);
        else mLogFile->setFileName(filename);
        mLogFile->open(QIODevice::WriteOnly);
        if (!mLogFile->isWritable())
        {
            QMessageBox::warning(this, "ClassCluster logging", "Error: unable to write to the selected file!");
            mWriteLogFileCheckBox->setChecked(false);
            return;
        }
        if (nullptr == mLogStream) mLogStream = new QTextStream(mLogFile);
        (*mLogStream) << "Time\tSeverity\tFunction\tFile\tMessage\n";
        const QList<QStringList>& messageBuffer = mModel.GetMessageBuffer();
        int currentMessageCount = messageBuffer.size();
        for (int i=0; i < currentMessageCount; ++i) (*mLogStream) << messageBuffer[i].join('\t') << '\n';
        mLoggerMutex->lock();
        for (int i = currentMessageCount; i < messageBuffer.size(); ++i) (*mLogStream) << messageBuffer[i].join('\t') << '\n';
        mLogStream->flush();
        mWriteToFile = true;
        mLoggerMutex->unlock();
        mFilenameEdit->setReadOnly(true);
        mFileDialogButton->setEnabled(false);
    }
}
