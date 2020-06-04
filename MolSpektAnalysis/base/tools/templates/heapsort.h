//
// C++ Interface: heapSort
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2006 - 2019
//
// Copyright: See README file that comes with this source code
//
//

#ifndef HEAPSORT
#define HEAPSORT



namespace utils
{
    template <class T> int* heapSort(T sortOperator, int N)
    {
        int i, n, m;
        int S1[N], *S2 = new int[N], B=0;
        for (i=0; i<N; i++) S1[i] = i;
        while (B>=0)
        {
            B=-1;
            for (n=0, m=2; m<=N; n++, m+=2)
            {
                if (m == N || sortOperator(S1[m-1], S1[m]))
                    i=m-1;
                else i=m;
                if (sortOperator(S1[i], S1[n]))
                {
                    B = S1[i];
                    S1[i] = S1[n];
                    S1[n] = B;
                }
            }
        }
        for (i=0; i<N; i++)
        {
            S2[S1[0]]=i;
            for (n=0, m=2; S1[n] != -1; m=2*((n=m)+1))
            {
                if (m < N)
                {
                    if (sortOperator(S1[m-1], S1[m])) m--;
                }
                else if (m==N) m--;
                else break;
                S1[n] = S1[m];
            }
            if (m >= N) S1[n] = -1;
        }
        return S2;
    }
}

#endif
