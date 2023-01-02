#ifndef VECTOR_H
#define VECTOR_H


class Vector
{
public:
    Vector();
    Vector(const Vector& other);
    Vector(const double X, const double Y, const double Z);

    Vector& operator=(const Vector& other);
    Vector operator +(const Vector& other) const;
    Vector& operator+=(const Vector& other);
    Vector operator -(const Vector& other) const;
    Vector& operator-=(const Vector& other);
    Vector operator *(const double factor) const;
    Vector& operator*=(const double factor);
    Vector operator /(const double factor) const;
    Vector& operator/=(const double factor);
    Vector operator*(const Vector& other) const;
    Vector& operator*=(const Vector& other);

    double length() const;
    double lengthSquared() const;
    double dot(const Vector& other) const;
    Vector cross(const Vector& other) const;

    inline bool operator==(const Vector& other) const
    {
        return x == other.x && y == other.y && z == other.z;
    }

    inline double X() const
    {
        return x;
    }

    inline double Y() const
    {
        return y;
    }

    inline double Z() const
    {
        return z;
    }

    inline void setX(double X)
    {
        x = X;
    }

    inline void setY(double Y)
    {
        y = Y;
    }

    inline void setZ(double Z)
    {
        z = Z;
    }

    inline void clear()
    {
        x = y = z = 0.0;
    }

    inline Vector unit()
    {
        return *this / length();
    }

private:
    double x, y, z;
};

inline Vector operator*(const double left, const Vector& right)
{
    return right * left;
}

#endif // VECTOR_H
