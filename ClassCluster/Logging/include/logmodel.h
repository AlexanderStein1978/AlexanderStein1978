//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//

#ifndef LOGMODEL_H
#define LOGMODEL_H


#include <QAbstractTableModel>

class LogList;


class LogModel : public QAbstractTableModel
{
public:
    LogModel(QObject *parent = 0);
    virtual ~LogModel();

    int rowCount(const QModelIndex &parent = QModelIndex()) const override;
    int columnCount(const QModelIndex &parent = QModelIndex()) const override;
    QVariant data (const QModelIndex &index, int role = Qt::DisplayRole) const override;
    QVariant headerData(int section, Qt::Orientation orientation, int role = Qt::DisplayRole) const override;

    void AddLogMessage(const QStringList& message);
    void SetMaxRows(const int maxValue);
    void SetMessageBuffer(LogList& messageBuffer);

    inline const LogList& GetMessageBuffer() const
    {
        return *mMessageBuffer;
    }

private:
    LogList* mMessageBuffer;
};

#endif // LOGMODEL_H
