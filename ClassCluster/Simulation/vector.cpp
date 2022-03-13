#include "vector.h"
#include <cmath>

Vector::Vector() : x(0.0), y(0.0), z(0.0)
{
}

Vector::Vector(const Vector &other) : x(other.x), y(other.y), z(other.z)
{
}

Vector::Vector(const double X, const double Y, const double Z) : x(X), y(Y), z(Z)
{
}

double Vector::length() const
{
    return sqrt(lengthSquared());
}

double Vector::lengthSquared() const
{
    return x*x + y*y + z*z;
}

Vector Vector::operator *(const double right) const
{
    return Vector(x*right, y*right, z*right);
}

Vector Vector::operator *(const Vector& right) const
{
    return Vector(x*right.x, y*right.y, z*right.z);
}

Vector& Vector::operator *=(const double factor)
{
    x *= factor;
    y *= factor;
    z *= factor;
    return *this;
}

Vector& Vector::operator *=(const Vector& other)
{
    x *= other.x;
    y *= other.y;
    z *= other.z;
    return *this;
}

Vector Vector::operator +(const Vector& other) const
{
    return Vector(x + other.x, y + other.y, z + other.z);
}

Vector& Vector::operator +=(const Vector& other)
{
    x += other.x;
    y += other.y;
    z += other.z;
    return *this;
}

Vector Vector::operator -(const Vector& other) const
{
    return Vector(x - other.x, y - other.y, z - other.z);
}

Vector& Vector::operator -=(const Vector& other)
{
    x -= other.x;
    y -= other.y;
    z -= other.z;
    return *this;
}

Vector Vector::operator /(const double right) const
{
    return Vector(x / right, y / right, z / right);
}

Vector& Vector::operator /=(const double right)
{
    x /= right;
    y /= right;
    z /= right;
    return *this;
}

Vector& Vector::operator =(const Vector& right)
{
    x = right.x;
    y = right.y;
    z = right.z;
    return *this;
}

double Vector::dot(const Vector &other) const
{
    return other.x * x + other.y * y + other.z * z;
}
