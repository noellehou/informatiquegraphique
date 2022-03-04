#include <iostream>
#include <vector>
#include "Vector.h"

class BBox {
public:
    BBox() {};
	BBox(const Vector& bmin, const Vector& bmax) : bmin(bmin), bmax(bmax) {};
	bool intersection(const Ray& d) const {
		double t_1_x = (bmin[0] - d.origin[0]) / d.direction[0];
		double t_2_x = (bmin[0] - d.origin[0]) / d.direction[0];
		double t_min_x = std::min(t_1_x, t_2_x);
		double t_max_x = std::max(t_1_x, t_2_x);

		double t_1_y = (bmin[1] - d.origin[1]) / d.direction[1];
		double t_2_y = (bmin[1] - d.origin[1]) / d.direction[1];
		double t_min_y = std::min(t_1_y, t_2_y);
		double t_max_y = std::max(t_1_y, t_2_y);

		double t_1_z = (bmin[2] - d.origin[2]) / d.direction[2];
		double t_2_z = (bmin[2] - d.origin[2]) / d.direction[2];
		double t_min_z = std::min(t_1_z, t_2_z);
		double t_max_z = std::max(t_1_z, t_2_z);

		if (std::min(std::min(t_max_x, t_max_y),t_max_z) - std::max(std::max(t_min_x,t_min_y),t_min_z) > 0) return true;
		return false;
	}
	Vector bmin, bmax;
};



class Object {
public:
    Object() {}
    virtual bool intersection(const Ray& d, Vector& P, Vector& N, double &t) const = 0;

    Vector albedo;
    bool miroir;
    bool transparency;
};


class Geometry : Object {
public :
    Geometry(const char* obj, double scaling, const Vector& offset, const Vector &couleur, bool mirror = false, bool transp = false);

    std::vector<int> faceGroup;
    std::vector<int> faces;
    std::vector<int> normalIds;
    std::vector<int> uvIds;
    std::vector<Vector> vertices;
    std::vector<Vector> normals;
    std::vector<Vector> uvs; //Vector en 3D mais on n'utilise que 2 composantes

    bool Intersection(const Ray& d, Vector& P, Vector& N, double &t) const;

private:
    BBox bb;

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
    const Vector &A, &B, &C;
};






class Scene  {
public:
    Scene() {};

    void addSphere(const Sphere& s) { objects.push_back((Object*)&s); }
    void addTriangle(const Triangle& s) { objects.push_back((Object*)&s); }
    void addGeometry(const Geometry& s) { objects.push_back((Object*)&s); }

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
