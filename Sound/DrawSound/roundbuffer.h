//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//

#pragma once


class RoundBuffer
{
public:
    struct BufferElement
    {
        double frequency, amplitude;
    };

    RoundBuffer(const int size);
    RoundBuffer(const RoundBuffer& other);
    virtual ~RoundBuffer();

    RoundBuffer& operator=(const RoundBuffer& right);

    virtual void setNextElement(const double frequency, const double amplitude);
    const BufferElement& element(const int index) const;

    inline int getNumElements() const
    {
        return mNumElements;
    }

private:
    BufferElement* m_elements;
    int m_size, m_currentIndex, mNumElements;
};
