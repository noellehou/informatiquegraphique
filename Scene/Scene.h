
#ifndef DEF_SCENE
#define DEF_SCENE

#include "../Sphere/Sphere.h"
#include "../Vector/Vector.h"
#include "../triangles/triangles.h"
#include "../Ray/Ray.h"


class Scene  {
public:
    Scene() {};

    void addSphere(const Sphere& s) { objects.push_back((Object*)&s); }
    //void addTriangle(const Triangle& s) { objects.push_back((Object*)&s); }
	void addMesh(TriangleMesh* m) {objects.push_back(m); }

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

#endif