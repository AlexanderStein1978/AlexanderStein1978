//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//

#include "vector.h"
#include <cmath>
#include <array>


Vector::Vector()
{
    memset(val, 0, sizeof(double[3]));
}

Vector::Vector(const Vector &other)
{
    memcpy(val, other.val, sizeof(double[3]));
}

Vector::Vector(const double X, const double Y, const double Z)
{
    val[x] = X;
    val[y] = Y;
    val[z] = Z;
}

double Vector::length() const
{
    return sqrt(lengthSquared());
}

double Vector::lengthSquared() const
{
    double result(0.0);
    for (int i=0; i < dimension; ++i) result += val[i] * val[i];
    return result;
}

Vector Vector::operator *(const double right) const
{
    Vector result;
    for (int i=0; i < dimension; ++i) result.val[i] = val[i] * right;
    return result;
}

Vector Vector::operator *(const Vector& right) const
{
    Vector result;
    for (int i=0; i < dimension; ++i) result.val[i] = val[i] * right.val[i];
    return result;
}

Vector& Vector::operator *=(const double factor)
{
    for (int i=0; i < dimension; ++i) val[i] *= factor;
    return *this;
}

Vector& Vector::operator *=(const Vector& other)
{
    for (int i=0; i < dimension; ++i) val[i] *= other.val[i];
    return *this;
}

Vector Vector::operator +(const Vector& other) const
{
    Vector result;
    for (int i=0; i < dimension; ++i) result.val[i] = val[i] + other.val[i];
    return result;
}

Vector& Vector::operator +=(const Vector& other)
{
    for (int i=0; i < dimension; ++i) val[i] += other.val[i];
    return *this;
}

Vector Vector::operator -(const Vector& other) const
{
    Vector result;
    for (int i=0; i < dimension; ++i) result.val[i] = val[i] - other.val[i];
    return result;
}

Vector& Vector::operator -=(const Vector& other)
{
    for (int i=0; i < dimension; ++i) val[i] -= other.val[i];
    return *this;
}

Vector Vector::operator /(const double right) const
{
    Vector result;
    for (int i=0; i < dimension; ++i) result.val[i] = val[i] / right;
    return result;
}

Vector& Vector::operator /=(const double right)
{
    for (int i=0; i < dimension; ++i) val[i] /= right;
    return *this;
}

Vector& Vector::operator =(const Vector& right)
{
    memcpy(val, right.val, sizeof(double[dimension]));
    return *this;
}

double Vector::dot(const Vector &other) const
{
    double result(0.0);
    for (int i=0; i < dimension; ++i) result += val[i] * other.val[i];
    return result;
}

Vector Vector::cross(const Vector& other) const
{
    return Vector(val[y] * other.val[z] - val[z] * other.val[y], val[z] * other.val[x] - val[x] * other.val[z], val[x] * other.val[y] - val[y] * other.val[x]);
}

Vector Vector::unit() const
{
    const double l = length();
    if (l == 0.0) return Vector();
    return *this / l;
}

void Vector::getSortOrder(int order[3]) const
{
    for (int i=0; i < dimension; ++i) order[i] = i;
    if (std::abs(val[0]) < std::abs(val[1])) std::swap(order[0], order[1]);
    if (std::abs(val[order[1]]) < std::abs(val[2])) std::swap(order[1], order[2]);
    if (std::abs(val[order[0]]) < std::abs(val[order[1]])) std::swap(order[0], order[1]);
}
