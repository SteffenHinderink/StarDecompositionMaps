#include "bary.h"

BarycentricVector& BarycentricVector::operator+=(const BarycentricVector& v) {
    for (int i = 0; i < v.vertices.size(); i++) {
        bool added = false;
        for (int j = 0; j < vertices.size(); j++) {
            if (vertices[j] == v.vertices[i]) {
                coordinates[j] += v.coordinates[i];
                added = true;
            }
        }
        if (!added) {
            vertices.push_back(v.vertices[i]);
            coordinates.push_back(v.coordinates[i]);
        }
    }
    return *this;
}

BarycentricVector operator+(BarycentricVector a, const BarycentricVector& b) {
    a += b;
    return a;
}

BarycentricVector& BarycentricVector::operator-=(const BarycentricVector& v) {
    for (int i = 0; i < v.vertices.size(); i++) {
        bool subtracted = false;
        for (int j = 0; j < vertices.size(); j++) {
            if (vertices[j] == v.vertices[i]) {
                coordinates[j] -= v.coordinates[i];
                subtracted = true;
            }
        }
        if (!subtracted) {
            vertices.push_back(v.vertices[i]);
            coordinates.push_back(-v.coordinates[i]);
        }
    }
    for (int i = 0; i < vertices.size(); i++) {
        if (coordinates[i] == 0) {
            vertices.erase(vertices.begin() + i);
            coordinates.erase(coordinates.begin() + i);
            i--;
        }
    }
    return *this;
}

BarycentricVector operator-(BarycentricVector a, const BarycentricVector& b) {
    a -= b;
    return a;
}

BarycentricVector& BarycentricVector::operator*=(const mpq_class& s) {
    for (int i = 0; i < vertices.size(); i++) {
        coordinates[i] *= s;
    }
    return *this;
}

BarycentricVector operator*(BarycentricVector v, const mpq_class& s) {
    v *= s;
    return v;
}

BarycentricVector operator*(const mpq_class& s, BarycentricVector v) {
    v *= s;
    return v;
}

BarycentricVector& BarycentricVector::operator/=(const mpq_class& s) {
    for (int i = 0; i < vertices.size(); i++) {
        coordinates[i] /= s;
    }
    return *this;
}

BarycentricVector operator/(BarycentricVector v, const mpq_class& s) {
    v /= s;
    return v;
}

bool operator==(const BarycentricVector& a, const BarycentricVector& b) {
    if (a.vertices.size() != b.vertices.size()) {
        return false;
    }
    for (int i = 0; i < a.vertices.size(); i++) {
        bool found = false;
        for (int j = 0; j < b.vertices.size(); j++) {
            if (a.vertices[i] == b.vertices[j]) {
                if (a.coordinates[i] != b.coordinates[j]) {
                    return false;
                }
                found = true;
            }
        }
        if (!found) {
            return false;
        }
    }
    return true;
}

bool operator!=(const BarycentricVector& a, const BarycentricVector& b) {
    return !(a == b);
}

std::ostream& operator<<(std::ostream& os, const BarycentricVector& v) {
    os << "(";
    for (int i = 0; i < v.vertices.size(); i++) {
        os << v.coordinates[i] << " * <" << v.vertices[i].idx() << ">";
        if (i < v.vertices.size() - 1) {
            os << ", ";
        }
    }
    os << ")";
    return os;
}
