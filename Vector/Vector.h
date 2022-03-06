
#ifndef DEF_VECTOR
#define DEF_VECTOR
#pragma once
#include <math.h>
#include <vector>
#include <random>


static std::default_random_engine engine;
static std::uniform_real_distribution<double> uniform(0,1);



class Vector {
public:
    Vector(double x=0, double y=0, double z=0){
        coord[0] = x;
        coord[1] = y;
        coord[2] = z;
    }

    const double& operator[](int i) const { return coord[i];}
    double& operator[](int i) { return coord[i]; }

    double getNorm2() {
        return coord[0]*coord[0] + coord[1]*coord[1] + coord[2] * coord[2];
    }
    void normalize() {
        double norm = sqrt(getNorm2());
        coord[0] /= norm;
        coord[1] /= norm;
        coord[2] /= norm;

    }
    Vector getNormalized() {
        Vector result(*this);
        result.normalize();
        return result;
    }

    Vector& operator+=(const Vector& b){
        coord[0] += b[0];
        coord[1] += b[1];
        coord[2] += b[2];
        return *this;
    }

private:
    double coord[3];

};

Vector operator+(const Vector& a, const Vector &b);
Vector operator-(const Vector& a, const Vector &b);
Vector operator*(double a, const Vector &b);
Vector operator*(const Vector &b, double a);
Vector operator*(const Vector &a, const Vector &b);
Vector operator/(const Vector& a, double b);
double dot(const Vector&a, const Vector& b);
double sqr(const double&a);
Vector cross(const Vector&a, const Vector& b);

Vector random_cos(const Vector &N);




#endif