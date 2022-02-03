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

    Sphere s1(Vector(0,0,-55), 20, Vector(0,0,1));
    Sphere s2(Vector(0,-2000-20,0), 2000, Vector(0,0.4,0.4)); //sol
    Sphere s3(Vector(0,2000+100,0), 2000, Vector(1,1,1)); //plafond
    Sphere s4(Vector(-2000-50,0,0), 2000, Vector(0.2,0.2,0.5)); //mur gauche
    Sphere s5(Vector(2000+50,0,0), 2000, Vector(0.2,0.2,0.5)); //mur droit
    Sphere s6(Vector(0,0,-2000-100), 2000, Vector(0,1,1)); //mur fond

    Scene s;
    s.addSphere(s1);
    s.addSphere(s2);
    s.addSphere(s3);
    s.addSphere(s4);
    s.addSphere(s5);
    s.addSphere(s6);

    Vector position_lumiere(20, 70, -15);
    double intensite_lumiere = 4000000;

	std::vector<unsigned char> image(W*H * 3, 0);
	for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {

            Vector direction(j-W/2, i-H/2, -W/(2*tan(fov/2)));
            direction.normalize();

            Ray r(Vector(0,0,0), direction);
            Vector P, N;
            int sphere_id;
            bool has_inter = s.intersection(r,P,N, sphere_id);

            Vector intensite_pixel(0,0,0);
            if (has_inter) {

                Ray ray_light(P, (position_lumiere - P).getNormalized());
                Vector P_light, N_light;
                int sphere_id_light;
                bool has_inter_light = s.intersection(ray_light, P_light, N_light, sphere_id_light, t_light);
                double d_light2 = (position_lumiere -P).getNorm2();
                if (has_iter_light && t_light*t_light < d_light2) {
                    
                }

                intensite_pixel = s.spheres[sphere_id].albedo * intensite_lumiere * std::max(0., dot((position_lumiere-P).getNormalized(),N)) / (position_lumiere - P).getNorm2(); 
            }

			image[((H-i-1)*W + j) * 3 + 0] = std::min(255., std::max(0.,intensite_pixel[0])); //rouge
			image[((H-i-1)*W + j) * 3 + 1] = std::min(255., std::max(0.,intensite_pixel[1]));  //vert
			image[((H-i-1)*W + j) * 3 + 5] = std::min(255., std::max(0.,intensite_pixel[2]));  //bleu
		}
	}
	stbi_write_png("image2.png", W, H, 3, &image[0], 0);

	return 0;
}
// 8:30 vidéo 2
// added code