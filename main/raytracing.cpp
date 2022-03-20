#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "../lib/stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "../lib/stb_image.h" //permet de lire des fichiers images


#include "../Vector/Vector.h"
#include "../Ray/Ray.h"
#include "../Sphere/Sphere.h"
#include "../triangles/triangles.h"
#include "../Scene/Scene.h"


#include <algorithm>
#include <stdio.h>
#include <math.h>

#include <cmath>
#include <thread>
//#include<omp.h>
//#define OMP_NUM_THREADS = 16

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

double sqr(const double&a) {
    return pow(a,2);
}

Vector cross(const Vector&a, const Vector& b) {
    return Vector(a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2], a[0]*b[1] - a[1]*b[0]);
}

Vector random_cos(const Vector &N) {

double r1 = uniform(engine);
double r2 = uniform(engine);
Vector direction_aleatoire_repere_local(cos(2*M_PI*r1)*sqrt(1-r2),sin(2*M_PI*r1)*sqrt(1-r2),sqrt(r2));
Vector aleatoire(uniform(engine)-0.5,uniform(engine)-0.5,uniform(engine)-0.5);
Vector tangent1 = cross(N,aleatoire); tangent1.normalize();
Vector tangent2 = cross(tangent1,N);

return direction_aleatoire_repere_local[2]*N + direction_aleatoire_repere_local[0] * tangent1 + direction_aleatoire_repere_local[1] * tangent2;

}


Vector getColor(Ray &r, const Scene &s, int nbrebonds,bool show_lights=true) {

    if (nbrebonds==0) return Vector(0,0,0);

    Vector P, N;
    int sphere_id;
    double t;
    bool has_inter = s.intersection(r,P,N, sphere_id, t);

    Vector intensite_pixel(0,0,0);
    if (has_inter) {

        /*if (sphere_id ==0){
            intensite_pixel = show_lights?(s.lumiere->albedo * s.intensite_lumiere):Vector(0.,0.,0.);
        } else*/
        if (s.objects[sphere_id]->miroir) {
            Vector direction_miroir = r.direction - 2*dot(N, r.direction)*N;
            Ray rayon_miroir(P + 0.001*N, direction_miroir);
            intensite_pixel = getColor(rayon_miroir, s, nbrebonds - 1);
        
        } else{
            if (s.objects[sphere_id]->transparency) {
                double n1 = 1;
                double n2 = 1.3; //verre
                Vector normale_pour_transparence(N);
                if (dot(r.direction, N) > 0 ) { //on sort de la sphere
                    n1 = 1.3;
                    n2 = 1;
                    normale_pour_transparence = 0 - N;
                }
                double radical = 1-sqr(n1/n2)*(1-sqr(dot(normale_pour_transparence,r.direction)));
                if (radical >0 ) {
                    Vector direction_refracte = (n1/n2) * (r.direction - dot(r.direction, normale_pour_transparence)*normale_pour_transparence) - normale_pour_transparence * sqrt(radical);
                    Ray rayon_refracte(P - 0.001*normale_pour_transparence, direction_refracte);
                    intensite_pixel = getColor(rayon_refracte, s, nbrebonds - 1);
                }
            } else {

        //Contribution de l'éclairage direct
            
            Ray ray_light(P + 0.01*N, (s.position_lumiere - P).getNormalized());
            Vector P_light, N_light;
            int sphere_id_light;
            double t_light;
            bool has_inter_light = s.intersection(ray_light, P_light, N_light, sphere_id_light, t_light);
            double d_light2 = (s.position_lumiere - P).getNorm2();
            if (has_inter_light && t_light*t_light < d_light2) {
                intensite_pixel = Vector(0,0,0);
            } 
            else{
                intensite_pixel = s.objects[sphere_id]->albedo / M_PI * s.intensite_lumiere * std::max(0., dot((s.position_lumiere-P).getNormalized(),N)) / d_light2;
            } 
             /*
            Vector axeOP = (P - s.lumiere->O).getNormalized();
            Vector dir_aleatoire = random_cos(axeOP);
            Vector point_aleatoire = dir_aleatoire * s.lumiere->R + s.lumiere->O; 
            Vector wi = (point_aleatoire - P).getNormalized();
            double d_light2 = (point_aleatoire - P).getNorm2();
            Vector Np = dir_aleatoire;

            Ray ray_light(P + 0.01*N, wi);
            Vector P_light, N_light;
            int sphere_id_light;
            double t_light;
            bool has_inter_light = s.intersection(ray_light, P_light, N_light, sphere_id_light, t_light);

            if (has_inter_light && t_light*t_light < d_light2*0.99) {
                intensite_pixel = Vector(0,0,0);
            } else {
                intensite_pixel = ( s.intensite_lumiere / (4*M_PI*d_light2) * std::max(0.,dot(N,wi)) * dot(Np, 0-wi) / dot(axeOP, dir_aleatoire) ) * s.objects[sphere_id]->albedo ;
                // Vector BRDF = s.objects[sphere_id]->albedo / M_PI;
                // double J = 1. * dot(Np, 0-wi) / d_light2;
                // double proba = dot(axeOP, dir_aleatoire) / (M_PI * s.lumiere->R*s.lumiere->R);
                // intensite_pixel = ( s.intensite_lumiere * std::max(0.,dot(N,wi)) * dot(Np, 0-wi) / dot(axeOP, dir_aleatoire) ) * s.objects[sphere_id]->albedo ;
            }
         

        //Ajout de la contribution de l'éclairage indirect
        
            Vector direction_aleatoire = random_cos(N);
            Ray rayon_aleatoire (P + 0.001*N, direction_aleatoire);


            intensite_pixel += getColor(rayon_aleatoire, s, nbrebonds - 1) * s.objects[sphere_id]->albedo;
        
*/
            }
           
        }       
    }
    return intensite_pixel;

}

int main() {
    
	int W = 100;
	int H = 100;
    const int nrays = 10;
	double fov = 60*M_PI/180;

    

    Sphere slum(Vector(15, 70, -30), 15, Vector(1.,1.,1.));
 
    //Sphere s1(Vector(-15,0,-35), 20, Vector(1.,1.,1.),true,false);
    //Sphere s2(Vector(15,0,20), 20, Vector(1,1,1),false,true);
    //Sphere s3(Vector(15,0,-75), 10, Vector(1,1,1),true);
    const char* nom_image = "test_mesh.png";

    Sphere sol(Vector(0,-2000-20,0), 2000, Vector(1,1,1)); //sol
    Sphere plafond(Vector(0,2000+100,0), 2000, Vector(1,1,1)); //plafond
    Sphere murgauche(Vector(-2000-50,0,0), 2000, Vector(1,1,1)); //mur gauche
    Sphere murdroit(Vector(2000+50,0,0), 2000, Vector(1,1,1)); //mur droit
    Sphere murfond(Vector(0,0,-2000-100), 2000, Vector(1,1,1)); //mur fond
    TriangleMesh mesh_test(Vector(255,255,255), false,false);
    mesh_test.readOBJ("/Users/noelle/Documents/2 - CENTRALE LYON/MOS/MOS 2.2 - Infographie/informatiquegraphique/mesh/dog.obj");
    //Triangle tri(Vector(-10,-10,-20),Vector(10,-10,-20),Vector(0,10,-20),Vector(1,0,0));

    Scene s; 
    //s.addSphere(slum);
    //s.addSphere(s1);
    //s.addSphere(s2);
    //s.addSphere(s3);
    //s.addTriangle(tri);
    s.addMesh(mesh_test);
    s.addSphere(sol);
    s.addSphere(plafond);
    s.addSphere(murgauche);
    s.addSphere(murdroit);
    s.addSphere(murfond);
    s.lumiere = &slum;
    s.intensite_lumiere = 1000000000; 
    Vector position_camera(0.,0.,50.);
    double focus_distance = 55;
    double aperture = 0.5;


	std::vector<unsigned char> image(W*H * 3, 0);

#pragma omp parallel for 
	for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {

            Vector color(0.,0.,0.);
            for (int k=0; k < nrays; k++) {
                // méthode de Box Muller pour l'anti-aliasing
                double r1 = uniform(engine);
                double r2 = uniform(engine);
                double R = sqrt(-2*log(r1));
                double dx = R*cos(2*M_PI*r2);
                double dy = R*sin(2*M_PI*r2);

                double dx_aperture = (uniform(engine) - 0.5) * aperture;
                double dy_aperture = (uniform(engine) - 0.5) * aperture;

                Vector direction(j-W/2 + 0.5 +dx, i-H/2 +0.5 +dy, -W/(2*tan(fov/2)));
                direction.normalize();

                Vector destination = position_camera + focus_distance * direction;
                Vector new_origin = position_camera + Vector(dx_aperture, dy_aperture,0);

                Ray r(new_origin, (destination-new_origin).getNormalized());

                color += getColor(r, s, 5) / nrays;
        }
			image[((H-i-1)*W + j) * 3 + 0] = std::min(255., std::max(0.,pow(color[0],1/2.2))); //rouge
			image[((H-i-1)*W + j) * 3 + 1] = std::min(255., std::max(0.,pow(color[1],1/2.2)));  //vert
			image[((H-i-1)*W + j) * 3 + 2] = std::min(255., std::max(0.,pow(color[2],1/2.2)));  //bleu
		}
	}

	stbi_write_png("test_chien.png", W, H, 3, &image[0], 0);

	return 0;
}
