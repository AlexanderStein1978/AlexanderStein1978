#ifndef RANDOM_SORTFUNCTOR_H
#define RANDOM_SORTFUNCTOR_H


class RandomSortFunctor
{
public:

    RandomSortFunctor(const double *const randomArray) : m_randomArray(randomArray)
    {
    }

    bool operator()(const int n1, const int n2)
    {
        if (n1 == -1) return false;
        if (n2 == -1) return true;
        return (m_randomArray[n1] > m_randomArray[n2]);
    }

private:
    const double *const m_randomArray;
};

#endif // RANDOM_SORTFUNCTOR_H
