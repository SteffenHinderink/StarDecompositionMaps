#pragma once

#include <Eigen/Core>
#include <gmpxx.h>

template <>
struct Eigen::NumTraits<mpq_class> : Eigen::GenericNumTraits<mpq_class> {
    using Real = mpq_class;
    using NonInteger = mpq_class;
    using Nested = mpq_class;

    static inline Real epsilon() { return 0; }
    static inline Real dummy_precision() { return 0; }
    static inline int digits10() { return 0; }

    enum {
        IsInteger = 0,
        IsSigned = 1,
        IsComplex = 0,
        RequireInitialization = 1,
        ReadCost = 6,
        AddCost = 150,
        MulCost = 100
    };
};

using Vector2q = Eigen::Matrix<mpq_class, 2, 1>;
using Matrix2q = Eigen::Matrix<mpq_class, 2, 2>;
using Vector3q = Eigen::Matrix<mpq_class, 3, 1>;
using Matrix3q = Eigen::Matrix<mpq_class, 3, 3>;
