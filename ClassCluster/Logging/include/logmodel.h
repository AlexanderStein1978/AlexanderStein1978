#ifndef LOGMODEL_H
#define LOGMODEL_H


#include <QAbstractTableModel>


class LogModel : public QAbstractTableModel
{
public:
    LogModel(QObject *parent = 0);

    int rowCount(const QModelIndex &parent = QModelIndex()) const override;
    int columnCount(const QModelIndex &parent = QModelIndex()) const override;
    QVariant data (const QModelIndex &index, int role = Qt::DisplayRole) const override;
    QVariant headerData(int section, Qt::Orientation orientation, int role = Qt::DisplayRole) const override;

    void AddLogMessage(QStringList& message);
    void SetMessageBuffer(QList<QStringList>& messageBuffer);

    inline const QList<QStringList>& GetMessageBuffer() const
    {
        return mMessageBuffer;
    }

private:
    QList<QStringList> mMessageBuffer;
};

#endif // LOGMODEL_H
