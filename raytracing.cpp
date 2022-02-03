#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#include "Vector.h"

#include <math.h>


//#define M_PI 3.1415926535897932
Vector operator+(const Vector& a, const Vector &b) {
    return Vector(a[0]+b[0], a[1]+b[1], a[2]+b[2]);
}
Vector operator-(const Vector& a, const Vector &b) {
    return Vector(a[0]-b[0], a[1]-b[1], a[2]-b[2]);
}
Vector operator*(double a, const Vector &b) {
    return Vector(a*b[0], a*b[1], a*b[2]);
}
Vector operator*(const Vector &a, const Vector &b) {
    return Vector(a[0]*b[0], a[1]*b[1], a[2]*b[2]);
}
Vector operator*(const Vector &b, double a) {
    return Vector(a*b[0], a*b[1], a*b[2]);
}
Vector operator/(const Vector& a, double b) {
    return Vector(a[0]/b, a[1]/b, a[2]/b);
}

double dot(const Vector&a, const Vector& b) {
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}



int main() {
	int W = 512;
	int H = 512;
	double fov = 60*M_PI/180;

    Sphere s1(Vector(0,0,-55), 20, Vector(1,0,0));
    Sphere s2(Vector(0,-2000-20,0), 2000, Vector(1,1,1)); //sol
    Sphere s3(Vector(0,2000+100,0), 2000, Vector(1,1,1)); //plafond
    Sphere s4(Vector(-2000-50,0,0), 2000, Vector(0,1,0)); //mur gauche
    Sphere s5(Vector(2000+50,0,0), 2000, Vector(0,0,1)); //mur droit
    Sphere s6(Vector(0,0,-2000-100), 2000, Vector(0,1,1)); //mur fond

    Scene s;
    s.addSphere(s1);
    s.addSphere(s2);
    s.addSphere(s3);
    s.addSphere(s4);
    s.addSphere(s5);
    s.addSphere(s6);

    Vector position_lumiere(15, 70, -30);
    double intensite_lumiere = 100000000;

	std::vector<unsigned char> image(W*H * 3, 0);
	for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {

            Vector direction(j-W/2, i-H/2, -W/(2*tan(fov/2)));
            direction.normalize();

            Ray r(Vector(0,0,0), direction);
            Vector P, N;
            int sphere_id;
            double t;
            bool has_inter = s.intersection(r,P,N, sphere_id, t);

            Vector intensite_pixel(0,0,0);
            if (has_inter) {

                Ray ray_light(P + 0.01*N, (position_lumiere - P).getNormalized());
                Vector P_light, N_light;
                int sphere_id_light;
                double t_light;
                bool has_inter_light = s.intersection(ray_light, P_light, N_light, sphere_id_light, t_light);
                double d_light2 = (position_lumiere - P).getNorm2();
                if (has_inter_light && t_light*t_light < d_light2) {
                    intensite_pixel = Vector(0,0,0);
                } 
                else{
                    intensite_pixel = s.spheres[sphere_id].albedo * intensite_lumiere * std::max(0., dot((position_lumiere-P).getNormalized(),N)) / d_light2;
                }
            }

			image[((H-i-1)*W + j) * 3 + 0] = std::min(255., std::max(0.,pow(intensite_pixel[0],1/2.2))); //rouge
			image[((H-i-1)*W + j) * 3 + 1] = std::min(255., std::max(0.,pow(intensite_pixel[1],1/2.2)));  //vert
			image[((H-i-1)*W + j) * 3 + 2] = std::min(255., std::max(0.,pow(intensite_pixel[2],1/2.2)));  //bleu
		}
	}
	stbi_write_png("output.png", W, H, 3, &image[0], 0);

	return 0;
}
// 21:10 vidéo 2
