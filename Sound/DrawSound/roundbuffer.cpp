#include "roundbuffer.h"

#include <cstring>


RoundBuffer::RoundBuffer(const int size) : m_elements(new BufferElement[size]), m_size(size), m_currentIndex(0u), mNumElements(0u)
{
}

RoundBuffer::RoundBuffer(const RoundBuffer& other) : m_elements(new BufferElement[other.m_size]), m_size(other.m_size), m_currentIndex(other.m_currentIndex), mNumElements(other.mNumElements)
{
    memcpy(m_elements, other.m_elements, m_size * sizeof(BufferElement));
}

RoundBuffer::~RoundBuffer()
{
    delete[] m_elements;
}

RoundBuffer & RoundBuffer::operator=(const RoundBuffer& right)
{
    m_elements = new BufferElement[right.m_size];
    memcpy(m_elements, right.m_elements, right.m_size * sizeof(BufferElement));
    m_size = right.m_size;
    m_currentIndex = right.m_currentIndex;
    mNumElements = right.mNumElements;

    return *this;
}

void RoundBuffer::setNextElement(const double frequency, const double amplitude)
{
    if (++m_currentIndex >= m_size) m_currentIndex = 0;
    m_elements[m_currentIndex].frequency = frequency;
    m_elements[m_currentIndex].amplitude = amplitude;
    if (mNumElements < m_size) ++mNumElements;
}

const RoundBuffer::BufferElement& RoundBuffer::element(const int index) const
{
    int i = m_currentIndex - index;
    if (i < 0) i += m_size;
    return m_elements[i];
}
