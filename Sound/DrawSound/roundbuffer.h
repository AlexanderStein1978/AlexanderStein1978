#pragma once


class RoundBuffer
{
public:
    struct BufferElement
    {
        double frequency, amplitude;
    };

    RoundBuffer(const int size) : m_elements(new BufferElement[size]), m_size(size), m_currentIndex(0u)
    {
    }

    ~RoundBuffer()
    {
        delete[] m_elements;
    }

    inline void setNextElement(const double frequency, const double amplitude)
    {
        if (++m_currentIndex >= m_size) m_currentIndex = 0;
        m_elements[m_currentIndex].frequency = frequency;
        m_elements[m_currentIndex].amplitude = amplitude;
    }

    inline const BufferElement& element(const int index) const
    {
        int i = m_currentIndex - index;
        if (i < 0) i += m_size;
        return m_elements[i];
    }

private:
    BufferElement* m_elements;
    int m_size, m_currentIndex;
};
