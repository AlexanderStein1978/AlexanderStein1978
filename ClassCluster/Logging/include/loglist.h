#pragma once


#include <QMutex>
#include "loglistitem.h"


class LogList
{
public:
    LogList();
    ~LogList();

    LogList(const LogList& copy) = delete;
    LogList& operator=(const LogList& right) = delete;

    void clear();
    void addElement(const QStringList& newElement);
    const QStringList getElement(const int index) const;
    void SetMaxSize(const int maxSize);

    inline int size() const
    {
        return mSize;
    }

    inline bool isMaxReached() const
    {
       return (mSize == mMaxSize);
    }

private:
    QMutex *mMutex;
    LogListItem *mFirst, *mLast, **mCurrent;
    int *mCurrentPos, mMaxSize, mSize;
};
