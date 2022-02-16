
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


class Ray{
public:
    Ray(const Vector& o, const Vector& d): origin(o), direction(d) {};
    Vector origin, direction;
};

class Object {
public:
    Object() {}
    virtual bool intersection(const Ray& d, Vector& P, Vector& N, double &t) const = 0;

    Vector albedo;
    bool miroir;
    bool transparency;
};

class Sphere : public Object{
public:
    Sphere(const Vector &origin, double rayon, const Vector &couleur, bool mirror = false, bool transparent = false): O(origin), R(rayon) {
        albedo = couleur;
        miroir = mirror;
        transparency = transparent;
    };

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

};

class Triangle : public Object {
public:
    Triangle(const Vector& A, const Vector &B, const Vector& C, const Vector &couleur, bool mirror = false, bool transparent = false): A(A), B(B), C(C) {
        albedo = couleur;
        miroir = mirror;
        transparency = transparent;

    };

    bool intersection(const Ray& d, Vector& P, Vector& N, double &t) const{
        
        N = cross(B-A,C-A).getNormalized();
        t = dot(C - d.origin,N)/ dot(d.direction,N);
        if (t<0) return false;
        P = d.origin + t*d.direction;
        Vector u = B-A;
        Vector v = C-A;
        Vector w = P-A;
        double m11 = u.getNorm2();
        double m12 = dot(u,v);
        double m22 = v.getNorm2();
        double detm = m11 * m22 - m12 * m12;

        double b11 = dot(w,u);
        double b21 = dot(w,v);
        double detb = b11 * m22 - b21 * m12;
        double beta = detb/detm; //coord barycentrique w.r.t à B

        double g12 = b11;
        double g22 = b21;
        double detg = m11 * g22 - m12*g12;
        double gamma = detg/detm; //coord barycentrique w.r.t à A

        double alpha = 1 - beta - gamma;
        if (alpha <0 || alpha > 1) return false;
        if (beta < 0 || beta > 1 ) return false;
        if (gamma <0 || gamma > 1) return false;

        return true;
    }
    Vector A, B, C;
};

class Scene  {
public:
    Scene() {};

    void addSphere(const Sphere& s) { objects.push_back(&s); }
    void addTriangle(const Triangle& s) { objects.push_back(&s); }

    bool intersection(const Ray& d, Vector& P, Vector& N, int &sphere_id, double &min_t) const {

        bool has_inter = false;
        min_t = 1E99;

        for (int i=0; i<objects.size();i++){
            Vector localP, localN;
            double t;
            bool local_has_inter = objects[i]->intersection(d, localP, localN, t);
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
    std::vector<const Object*> objects;
    Sphere *lumiere;
    double intensite_lumiere;
};
