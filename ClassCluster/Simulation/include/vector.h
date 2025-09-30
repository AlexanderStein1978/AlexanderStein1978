//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//

#ifndef VECTOR_H
#define VECTOR_H

#include <string.h>


class Vector
{
public:
    enum Dimension{x, y, z, dimension};
    enum SortOrder{largest, mid, smallest};

    Vector();
    Vector(const Vector& other);
    Vector(const double X, const double Y, const double Z);

    Vector& operator=(const Vector& other);
    Vector operator+(const Vector& other) const;
    Vector& operator+=(const Vector& other);
    Vector operator-(const Vector& other) const;
    Vector& operator-=(const Vector& other);
    Vector operator*(const double factor) const;
    Vector& operator*=(const double factor);
    Vector operator/(const double factor) const;
    Vector& operator/=(const double factor);
    Vector operator*(const Vector& other) const;
    Vector& operator*=(const Vector& other);

    double length() const;
    double lengthSquared() const;
    double dot(const Vector& other) const;
    Vector cross(const Vector& other) const;
    Vector unit() const;
    void getSortOrder(int order[dimension]) const;

    inline bool operator==(const Vector& other) const
    {
        return (memcmp(val, other.val, sizeof(double[dimension])) == 0);
    }

    inline double X() const
    {
        return val[x];
    }

    inline double Y() const
    {
        return val[y];
    }

    inline double Z() const
    {
        return val[z];
    }

    inline void setX(double X)
    {
        val[x] = X;
    }

    inline void setY(double Y)
    {
        val[y] = Y;
    }

    inline void setZ(double Z)
    {
        val[z] = Z;
    }

    inline void clear()
    {
        memset(val, 0, sizeof(double[dimension]));
    }

    inline double& operator[](int i)
    {
        return val[i];
    }

    inline const double& operator[](const int i) const
    {
        return val[i];
    }

private:
    double val[dimension];
};

inline Vector operator*(const double left, const Vector& right)
{
    return right * left;
}

#endif // VECTOR_H
