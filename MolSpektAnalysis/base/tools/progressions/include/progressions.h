//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//

#ifndef PROGRESSIONS_H
#define PROGRESSIONS_H


struct IntProg;
struct FLI;
struct Marker;


class Progressions
{
public:
    Progressions();
    ~Progressions();
	void Insert(const IntProg *Prog, const int &G);
    void Insert(Marker **marker, const int &N, const int &G);
    void RotateB();
    void RotateF();
    void GetMarker(Marker **&marker, int &N);
    void Clear();

private:
	void Insert(FLI *newElement);
	
    FLI *P;
};

#endif
