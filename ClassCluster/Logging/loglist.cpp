//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//

#include "loglist.h"
#include <limits>


LogList::LogList() : mMutex(new QMutex), mFirst(nullptr), mLast(nullptr), mCurrent(new LogListItem*), mCurrentPos(new int), mMaxSize(std::numeric_limits<int>::max()), mSize(0)
{
    *mCurrentPos = -1;
}

LogList::~LogList()
{
    clear();
    delete mCurrent;
    delete mCurrentPos;
    delete mMutex;
}

void LogList::clear()
{
    while (nullptr != mFirst)
    {
        LogListItem *toDel = mFirst;
        mFirst = mFirst->next;
        delete toDel;
    }
    (*mCurrent) = mFirst = mLast = nullptr;
    mSize = 0;
    *mCurrentPos = -1;
}

void LogList::addElement(const QStringList& newElement)
{
    QMutexLocker lock(mMutex);
    if (nullptr == mFirst)
    {
        mFirst = new LogListItem;
        mLast = mFirst;
        mFirst->next = nullptr;
        mFirst->prev = nullptr;
        mSize = 1;
    }
    else if (mSize < mMaxSize)
    {
        LogListItem *newElement = new LogListItem;
        newElement->next = mFirst;
        mFirst->prev = newElement;
        mFirst = newElement;
        ++mSize;
        if (-1 != *mCurrentPos) ++(*mCurrentPos);
    }
    else
    {
        mLast->next = mFirst;
        mFirst->prev = mLast;
        mLast->message.clear();
        mFirst = mLast;
        mLast = mLast->prev;
        mFirst->prev = nullptr;
        mLast->next = nullptr;
        if (*mCurrentPos == mSize - 1) *mCurrent = mLast;
    }
    mFirst->message = newElement;
}

const QStringList LogList::getElement(const int index) const
{
    if (nullptr == mFirst)
    {
        return QStringList();
    }
    QMutexLocker lock(mMutex);
    if (-1 == *mCurrentPos)
    {
        if (index > mSize / 2)
        {
            *mCurrentPos = mSize - 1;
            *mCurrent = mLast;
        }
        else
        {
            *mCurrentPos = 0;
            *mCurrent = mFirst;
        }
    }
    while (*mCurrentPos > index && nullptr != (*mCurrent)->prev)
    {
        *mCurrent = (*mCurrent)->prev;
        --(*mCurrentPos);
    }
    while (*mCurrentPos < index && nullptr != (*mCurrent)->next)
    {
        *mCurrent = (*mCurrent)->next;
        ++(*mCurrentPos);
    }
    return (*mCurrent)->message;
}

void LogList::SetMaxSize(const int maxSize)
{
    mMaxSize = maxSize;
    QMutexLocker lock(mMutex);
    while (mSize > maxSize)
    {
        LogListItem* toDel = mLast;
        mLast = mLast->prev;
        mLast->next = nullptr;
        delete toDel;
        --mSize;
    }
    if (*mCurrentPos >= mSize)
    {
        *mCurrentPos = mSize - 1;
        *mCurrent = mLast;
    }
}


