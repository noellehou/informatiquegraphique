#include <iostream>
#include <vector>
#include "Vector.h"

class Geometry : Object {
public :
    Geometry(const char* obj, double scaling, const Vector& offset, const Vector &couleur, bool mirror = false, bool transp = false);

    std::vector<int> faceGroup;
    std::vector<int> faces;
    std::vector<int> normalIds;
    std::vector<int> uvIds;
    std::vector<Vector> vertices;
    std::vector<Vector> normals;
    std::vector<Vector< uvs; //Vector en 3D mais on n'utilise que 2 composantes

    bool Intersection(const Ray& d, Vector& P, Vector& N, double &t) const;
};