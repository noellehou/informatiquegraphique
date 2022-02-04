
#pragma once
#include <math.h>
#include <vector>
#include <random>
#define M_PI 3.1415926535897932

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


class Ray{
public:
    Ray(const Vector& o, const Vector& d): origin(o), direction(d) {};
    Vector origin, direction;
};

class Sphere {
public:
    Sphere(const Vector &origin, double rayon, const Vector &couleur, bool mirror = false, bool transparent = false): O(origin), R(rayon), albedo(couleur), miroir(mirror), transparency(transparent) {};

    bool intersection(const Ray& d, Vector& P, Vector& N, double &t) const {

    // resout a*t^2 + b*t + c = 0

    double a = 1;
    double b = 2 * dot(d.direction, d.origin - O);
    double c = (d.origin - O).getNorm2() - R*R;

    double delta = b*b - 4*a*c;
    if (delta <0) return false;
    double t1 = (-b - sqrt(delta))/(2*a);
    double t2 = (-b+sqrt(delta))/(2*a);

    if (t2<0) return false;

    if (t1>0)
        t = t1;
    else 
        t = t2;

    P = d.origin + t*d.direction;
    N = (P-O).getNormalized();
    return true;

}
    Vector O;
    double R;
    Vector albedo;
    bool miroir;
    bool transparency;
};

class Scene {
public:
    Scene() {};

    void addSphere(const Sphere& s) { spheres.push_back(s); }
    bool intersection(const Ray& d, Vector& P, Vector& N, int &sphere_id, double &min_t) const {

        bool has_inter = false;
        min_t = 1E99;

        for (int i=0; i<spheres.size();i++){
            Vector localP, localN;
            double t;
            bool local_has_inter = spheres[i].intersection(d, localP, localN, t);
            if (local_has_inter){
                has_inter = true;
                if (t < min_t) {
                    min_t = t;
                    P = localP;
                    N = localN;
                    sphere_id = i;
                }
            }
        }
        return has_inter;
    }
    std::vector<Sphere> spheres;
    Sphere *lumiere;
    double intensite_lumiere;
};
