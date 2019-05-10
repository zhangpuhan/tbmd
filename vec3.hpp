//
// Created by puhan on 2/4/19.
//

#include <cmath>
#include <complex>
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <math.h>


#define VEC3_NDIMS 3

// Thin wrapper that adds bounds checking to std::vector
template<class T>
class Vec : public std::vector<T> {
public:
    // using std::vector<T>::vector;
    // Be explicit to avoid annoying Xcode errors
    Vec<T>(): std::vector<T>() {}
    Vec<T>(int size): std::vector<T>(size) {}
    Vec<T>(int size, T const& val): std::vector<T>(size, val) {}

    //Vec<T>(std::initializer_list<T> l): std::vector<T>(l) {}

    template <class InputIterator>
    Vec<T>(InputIterator first, InputIterator last): std::vector<T>(first, last) {}

    T& operator[](int i) {
        return std::vector<T>::at(i);
    }

    const T& operator[](int i) const {
        return std::vector<T>::at(i);
    }
};



template<typename T>
class Vec3 {
public:
    T x=0, y=0, z=0;

    Vec3(void) { };

    constexpr Vec3(T x, T y, T z): x(x), y(y), z(z) {}

    constexpr Vec3<T> operator-() const {
        return {-x, -y, -z};
    }

    constexpr Vec3<T> operator+(Vec3<T> that) const {
        return {x+that.x, y+that.y, z+that.z};
    }

    constexpr Vec3<T> operator-(Vec3<T> that) const {
        return {x-that.x, y-that.y, z-that.z};
    }

    constexpr Vec3<T> operator*(T a) const {
        return {x*a, y*a, z*a};
    }

    constexpr friend Vec3<T> operator *(T x, Vec3<T> y) {
        return y*x;
    }

    constexpr Vec3<T> operator/(T a) const {
        return {x/a, y/a, z/a};
    }

    constexpr T dot(Vec3<T> const& that) const {
        return x*that.x + y*that.y + z*that.z;
    }

    constexpr Vec3<T> cross(Vec3<T> const& that) const {
        return {y*that.z-z*that.y, z*that.x-x*that.z, x*that.y-y*that.x};
    }

    constexpr T norm2() const {
        return dot(*this);
    }

    constexpr T norm() const {
        return sqrt(norm2());
    }

    constexpr Vec3<T> normalized() const {
        return *this / norm();
    }

    void operator+=(Vec3<T> that) {
        x += that.x;
        y += that.y;
        z += that.z;
    }

    void operator-=(Vec3<T> that) {
        x -= that.x;
        y -= that.y;
        z -= that.z;
    }

    void operator*=(T a) {
        x *= a;
        y *= a;
        z *= a;
    }

    void operator/=(T a) {
        x /= a;
        y /= a;
        z /= a;
    }

    template <typename S>
    constexpr operator Vec3<S>() const {
        return Vec3<S>(x, y, z);
    }

    friend std::ostream& operator<< (std::ostream& os, Vec3<T> const& v) {
        return os << "<x=" << v.x << ", y=" << v.y << ", z=" << v.z << ">";
    }

    T & operator()(unsigned int k) {
        switch(k % 3) {
            case 0:
                return x;
                break;
            case 1:
                return y;
                break;
            case 2:
                return z;
                break;
            default:
                return x;
                break;
        }
    }

};

template <typename S>
constexpr Vec3<S> real(Vec3<std::complex<S>> v) {
    return {std::real(v.x), std::real(v.y), std::real(v.z)};
}

template <typename S>
constexpr Vec3<S> imag(Vec3<std::complex<S>> v) {
    return {std::imag(v.x), std::imag(v.y), std::imag(v.z)};
}

typedef Vec3<float> fvec3;
typedef Vec3<double> vec3;
