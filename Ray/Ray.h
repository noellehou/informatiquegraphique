#ifndef DEF_RAY
#define DEF_RAY


class Ray{
public:
    Ray(const Vector& o, const Vector& d): origin(o), direction(d) {};
    Vector origin, direction;
};

#endif