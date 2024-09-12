#pragma once

#include <gmpxx.h>

#include "mesh.h"

class BarycentricVector {
public:
    BarycentricVector& operator+=(const BarycentricVector& v);

    friend BarycentricVector operator+(BarycentricVector a, const BarycentricVector& b);

    BarycentricVector& operator-=(const BarycentricVector& v);

    friend BarycentricVector operator-(BarycentricVector a, const BarycentricVector& b);

    BarycentricVector& operator*=(const mpq_class& s);

    friend BarycentricVector operator*(BarycentricVector v, const mpq_class& s);

    friend BarycentricVector operator*(const mpq_class& s, BarycentricVector v);

    BarycentricVector& operator/=(const mpq_class& s);

    friend BarycentricVector operator/(BarycentricVector v, const mpq_class& s);

    friend bool operator==(const BarycentricVector& a, const BarycentricVector& b);

    friend bool operator!=(const BarycentricVector& a, const BarycentricVector& b);

    friend std::ostream& operator<<(std::ostream& os, const BarycentricVector& v);

    std::vector<Vertex> vertices;
    std::vector<mpq_class> coordinates;
};
