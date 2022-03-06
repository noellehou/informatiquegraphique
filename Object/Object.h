#ifndef DEF_OBJECT
#define DEF_OBJECT
#include "../Vector/Vector.h"
#include "../Ray/Ray.h"


class Object {
public:
    Object() {}
    virtual bool intersection(const Ray& d, Vector& P, Vector& N, double &t) const = 0;

    Vector albedo;
    bool miroir;
    bool transparency;
};

#endif