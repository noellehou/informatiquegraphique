#ifndef DEF_SPHERE
#define DEF_SPHERE

#include "../Vector/Vector.h"
#include "../Object/Object.h"
#include "../Ray/Ray.h"



#include <cmath>

class Sphere : public Object{
public:
    Sphere(const Vector &origin, double rayon, const Vector &couleur, bool mirror = false, bool transparent = false): O(origin), R(rayon) {
        albedo = couleur;
        miroir = mirror;
        transparency = transparent;
    };

    bool intersection(const Ray& d, Vector& P, Vector& N, double &t) const {
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
};
    Vector O;
    double R;

};

#endif