/*VERSION DU 07/02

---Fonctionnalités---
x Affichage d'une sphère
x Affichage d'une scène de sphères avec des couleurs différentes
x Eclairage de la sphère par une source
x Ombre portée des sphères
X Surface diffuse
X Surface Speculaire 
x Eclairage indirect avec Monte Carlo
x Chronomètre
x Anti-aliasing
x Multithread
x Ombres douces / soruces de lumière diffuses
x Profondeur de champ/ Diaphragme camera


-----Bugs-----
x Surfaces transparentes

*/


#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "lib/stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "lib/stb_image.h"

#include <list>
#include <iostream>

#include <random>

#include <chrono>
using namespace std::chrono;

std::default_random_engine engine; 
std::uniform_real_distribution<double> uniform(0,1);

//Pour lecture fichiers obj
#include <string>
#include <iostream>
#include <stdio.h>
#include <algorithm>


//DEFINITION DES CLASSES


//VECTOR (X,Y,Z)
class Vector{
	public: 
	explicit Vector(double x=0 , double y=0, double z=0){
		data[0] = x;
		data[1] = y; 
		data[2] = z;
		
	}
	double norm2() const {
		return data[0]*data[0] + data[1]*data[1] + data[2]*data[2];
	}
	double norm() const {
		return sqrt(norm2());
	}
	void normalize(){
		double n = norm();

		data[0] /= n;
		data[1] /= n;
		data[2] /=n;
	}



	double operator[](int i) const {return data[i];}; //lecture seule
	double& operator[](int i){return data[i];}; //modification
	double data[3];
};

Vector operator+(const Vector& a, const Vector& b){
	return Vector(a[0]+b[0], a[1]+b[1], a[2]+b[2]);
}

Vector operator-(const Vector&a, const Vector& b){
	return Vector(a[0]-b[0], a[1]-b[1], a[2]-b[2]);
}

Vector operator*(const double a, const Vector& b) {
    return Vector(a*b[0], a*b[1], a*b[2]);
}
Vector operator*(const Vector& a, const double b) {
    return Vector(a[0]*b, a[1]*b, a[2]*b);
}
Vector operator/(const Vector& a, const double b) {
    return Vector(a[0] / b, a[1] / b, a[2] / b);
}
double dot(const Vector& a, const Vector& b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}
Vector cross(const Vector& a, const Vector& b) {
    return Vector(a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]);
}

Vector operator*(const Vector& a, const Vector& b){
	return Vector(a[0]*b[0],a[1]*b[1],a[2]*b[2]);
}




//RAYON (ORIGINE,DIRECTION)
class Ray {
public:
	Ray(const Vector& origin, const Vector&  direction) : C(origin),u(direction){

};

Vector C,u;

};


// OBJET
class Objet {
public: 
	Objet() {};
	//méthode virtuelle de l'objet englobant
	virtual bool intersect(const Ray& r , Vector& P, Vector& N, double& t) const =0;


	Vector albedo;
	bool objet_miroir;
	bool objet_transparent;

};
 

//SPHERE (CENTRE, RAYON, ALBEDO, MIROIR, TRANSPARENCE) 
class Sphere : public Objet {
public:
	Sphere( const Vector& O, const double R, const Vector& couleur,bool miroir = false, bool transp = false) : O(O),R(R){
		albedo = couleur;
		objet_miroir = miroir; 
		objet_transparent = transp;

};

	
	//Routine d'intersection qui teste si un rayon touche la sphère
	bool intersect(const Ray& r , Vector& P, Vector& N, double& t) const{
		//solve a*t^2 +b*t +c = 0
		double a =1;
		double b= 2*dot(r.u , r.C - O );
		double c= (r.C - O).norm2()- R*R;


		double delta = b*b - 4*a*c;
		if (delta < 0){
			return false ;
		}
		double t1 = (-b - sqrt(delta))/(2*a);
		double t2 = (-b + sqrt(delta))/(2*a);

		if (t2 <0){
			return false;
		}


		if (t1>0){
			t = t1;
		}
		else{
			t=t2;
		}

		P= r.C + t*r.u;
		N=(P - O);
		N.normalize();

		return true;

	}
		double R=R;
		Vector O=O;
		
};

class Triangle : public Objet {
public:
	Triangle(const Vector &A, const Vector &B, const Vector &C, const Vector &couleur, bool miroir =false, bool transp = false) : A(A), B(B), C(C) {
		albedo = couleur;
		objet_miroir = miroir;
		objet_transparent = transp;
	};

	bool intersect(const Ray& r , Vector& P, Vector& N, double& t) const {

		N = cross(B-A,C-A);
		N.normalize();
		t = dot(C-r.C,N) / dot(r.u,N);
		if (t < 0) return false; 

		P = r.C + t*r.u;
		Vector a = B-A;
		Vector b = C-A;
		Vector c = P-A;

		double m11 = a.norm2();
		double m12 = dot(a,b);
		double m22 = b.norm2();
		double detm = m11*m22 - m12*m12;

		double b11 = dot(c,a);
		double b21 = dot(c,b);
		double detb = b11*m22 - b21*m12;
		double beta = detb/detm; //coordonnée barycentrique B

		double g12 = b11;
		double g22 = b21;
		double detg = m11*g22 - m12*g12;
		double gamma = detg/detm; //coordonnée barycentrique C

		double alpha = 1 - beta - gamma;
		if(alpha<0 || alpha>1) return false;
		if(beta<0 || beta>1) return false;
		if(gamma < 0 || gamma>1) return false;
		if(alpha+beta+gamma>1) return false; //incorrect --corrigé plus tard


		return true;
	}

	Vector A,B,C;

};


// BOITE ENGLOBANTE

class BoiteEnglobante{
public:
BoiteEnglobante(){};
BoiteEnglobante(const Vector& mins_b, const Vector & maxs_b): mins_b(mins_b),maxs_b(maxs_b){};


bool intersection_boite(const Ray& r)const{

	//Plan x
	double tx_1 = (mins_b[0]-r.C[0])/r.u[0];
	double tx_2 = (maxs_b[0]-r.C[0])/r.u[0];
	double tx_min = std::min(tx_1,tx_2);
	double tx_max = std::max(tx_1,tx_2);

	//Plan y
	double ty_1 = (mins_b[1]-r.C[1])/r.u[1];
	double ty_2 = (maxs_b[1]-r.C[1])/r.u[1];
	double ty_min = std::min(ty_1,ty_2);
	double ty_max = std::max(ty_1,ty_2);

	//Plan z
	double tz_1 = (mins_b[2]-r.C[2])/r.u[2];
	double tz_2 = (maxs_b[2]-r.C[2])/r.u[2];
	double tz_min = std::min(tz_1,tz_2);
	double tz_max = std::max(tz_1,tz_2);

	//test d'intersection des intervalles tx,ty,tz
	
	double min_des_maxs= std::min(std::min(tx_max,ty_max),tz_max);
	double max_des_mins = std::max(std::min(tx_min,ty_min),tz_min);


	if(min_des_maxs-max_des_mins>0){ //Si l'intersection des intervalles est positives la boite est touchée
		return true;
	}
	return false;
	
} 

Vector mins_b;
Vector maxs_b;

};

// OBJ FILE READER

// INDICE TRIANGLES
class TriangleIndices {
public:
	TriangleIndices(int vtxi = -1, int vtxj = -1, int vtxk = -1, int ni = -1, int nj = -1, int nk = -1, int uvi = -1, int uvj = -1, int uvk = -1, int group = -1, bool added = false) : vtxi(vtxi), vtxj(vtxj), vtxk(vtxk), uvi(uvi), uvj(uvj), uvk(uvk), ni(ni), nj(nj), nk(nk), group(group) {
	};
	int vtxi, vtxj, vtxk; // indices within the vertex coordinates array
	int uvi, uvj, uvk;  // indices within the uv coordinates array
	int ni, nj, nk;  // indices within the normals array
	int group;       // face group
};


// MESH
class  TriangleMesh :Objet {
public:
  ~TriangleMesh() {}
	TriangleMesh(const Vector &couleur, bool miroir =false, bool transp = false) {
		
		objet_miroir = miroir;
		objet_transparent = transp;
		albedo = couleur;
	};
	
	void readOBJ(const char* obj) {

		char matfile[255];
		char grp[255];

		FILE* f;
		f = fopen(obj, "r");
		int curGroup = -1;
		while (!feof(f)) {
			char line[255];
			if (!fgets(line, 255, f)) break;

			std::string linetrim(line);
			linetrim.erase(linetrim.find_last_not_of(" \r\t") + 1);
			strcpy(line, linetrim.c_str());

			if (line[0] == 'u' && line[1] == 's') {
				sscanf(line, "usemtl %[^\n]\n", grp);
				curGroup++;
			}

			if (line[0] == 'v' && line[1] == ' ') {
				Vector vec;

				Vector col;
				if (sscanf(line, "v %lf %lf %lf %lf %lf %lf\n", &vec[0], &vec[1], &vec[2], &col[0], &col[1], &col[2]) == 6) {
					col[0] = std::min(1., std::max(0., col[0]));
					col[1] = std::min(1., std::max(0., col[1]));
					col[2] = std::min(1., std::max(0., col[2]));

					vertices.push_back(vec);
					vertexcolors.push_back(col);

				} else {
					sscanf(line, "v %lf %lf %lf\n", &vec[0], &vec[1], &vec[2]);
					vertices.push_back(vec);
				}
			}
			if (line[0] == 'v' && line[1] == 'n') {
				Vector vec;
				sscanf(line, "vn %lf %lf %lf\n", &vec[0], &vec[1], &vec[2]);
				normals.push_back(vec);
			}
			if (line[0] == 'v' && line[1] == 't') {
				Vector vec;
				sscanf(line, "vt %lf %lf\n", &vec[0], &vec[1]);
				uvs.push_back(vec);
			}
			if (line[0] == 'f') {
				TriangleIndices t;
				int i0, i1, i2, i3;
				int j0, j1, j2, j3;
				int k0, k1, k2, k3;
				int nn;
				t.group = curGroup;

				char* consumedline = line + 1;
				int offset;

				nn = sscanf(consumedline, "%u/%u/%u %u/%u/%u %u/%u/%u%n", &i0, &j0, &k0, &i1, &j1, &k1, &i2, &j2, &k2, &offset);
				if (nn == 9) {
					if (i0 < 0) t.vtxi = vertices.size() + i0; else	t.vtxi = i0 - 1;
					if (i1 < 0) t.vtxj = vertices.size() + i1; else	t.vtxj = i1 - 1;
					if (i2 < 0) t.vtxk = vertices.size() + i2; else	t.vtxk = i2 - 1;
					if (j0 < 0) t.uvi = uvs.size() + j0; else	t.uvi = j0 - 1;
					if (j1 < 0) t.uvj = uvs.size() + j1; else	t.uvj = j1 - 1;
					if (j2 < 0) t.uvk = uvs.size() + j2; else	t.uvk = j2 - 1;
					if (k0 < 0) t.ni = normals.size() + k0; else	t.ni = k0 - 1;
					if (k1 < 0) t.nj = normals.size() + k1; else	t.nj = k1 - 1;
					if (k2 < 0) t.nk = normals.size() + k2; else	t.nk = k2 - 1;
					indices.push_back(t);
				} else {
					nn = sscanf(consumedline, "%u/%u %u/%u %u/%u%n", &i0, &j0, &i1, &j1, &i2, &j2, &offset);
					if (nn == 6) {
						if (i0 < 0) t.vtxi = vertices.size() + i0; else	t.vtxi = i0 - 1;
						if (i1 < 0) t.vtxj = vertices.size() + i1; else	t.vtxj = i1 - 1;
						if (i2 < 0) t.vtxk = vertices.size() + i2; else	t.vtxk = i2 - 1;
						if (j0 < 0) t.uvi = uvs.size() + j0; else	t.uvi = j0 - 1;
						if (j1 < 0) t.uvj = uvs.size() + j1; else	t.uvj = j1 - 1;
						if (j2 < 0) t.uvk = uvs.size() + j2; else	t.uvk = j2 - 1;
						indices.push_back(t);
					} else {
						nn = sscanf(consumedline, "%u %u %u%n", &i0, &i1, &i2, &offset);
						if (nn == 3) {
							if (i0 < 0) t.vtxi = vertices.size() + i0; else	t.vtxi = i0 - 1;
							if (i1 < 0) t.vtxj = vertices.size() + i1; else	t.vtxj = i1 - 1;
							if (i2 < 0) t.vtxk = vertices.size() + i2; else	t.vtxk = i2 - 1;
							indices.push_back(t);
						} else {
							nn = sscanf(consumedline, "%u//%u %u//%u %u//%u%n", &i0, &k0, &i1, &k1, &i2, &k2, &offset);
							if (i0 < 0) t.vtxi = vertices.size() + i0; else	t.vtxi = i0 - 1;
							if (i1 < 0) t.vtxj = vertices.size() + i1; else	t.vtxj = i1 - 1;
							if (i2 < 0) t.vtxk = vertices.size() + i2; else	t.vtxk = i2 - 1;
							if (k0 < 0) t.ni = normals.size() + k0; else	t.ni = k0 - 1;
							if (k1 < 0) t.nj = normals.size() + k1; else	t.nj = k1 - 1;
							if (k2 < 0) t.nk = normals.size() + k2; else	t.nk = k2 - 1;
							indices.push_back(t);
						}
					}
				}

				consumedline = consumedline + offset;

				while (true) {
					if (consumedline[0] == '\n') break;
					if (consumedline[0] == '\0') break;
					nn = sscanf(consumedline, "%u/%u/%u%n", &i3, &j3, &k3, &offset);
					TriangleIndices t2;
					t2.group = curGroup;
					if (nn == 3) {
						if (i0 < 0) t2.vtxi = vertices.size() + i0; else	t2.vtxi = i0 - 1;
						if (i2 < 0) t2.vtxj = vertices.size() + i2; else	t2.vtxj = i2 - 1;
						if (i3 < 0) t2.vtxk = vertices.size() + i3; else	t2.vtxk = i3 - 1;
						if (j0 < 0) t2.uvi = uvs.size() + j0; else	t2.uvi = j0 - 1;
						if (j2 < 0) t2.uvj = uvs.size() + j2; else	t2.uvj = j2 - 1;
						if (j3 < 0) t2.uvk = uvs.size() + j3; else	t2.uvk = j3 - 1;
						if (k0 < 0) t2.ni = normals.size() + k0; else	t2.ni = k0 - 1;
						if (k2 < 0) t2.nj = normals.size() + k2; else	t2.nj = k2 - 1;
						if (k3 < 0) t2.nk = normals.size() + k3; else	t2.nk = k3 - 1;
						indices.push_back(t2);
						consumedline = consumedline + offset;
						i2 = i3;
						j2 = j3;
						k2 = k3;
					} else {
						nn = sscanf(consumedline, "%u/%u%n", &i3, &j3, &offset);
						if (nn == 2) {
							if (i0 < 0) t2.vtxi = vertices.size() + i0; else	t2.vtxi = i0 - 1;
							if (i2 < 0) t2.vtxj = vertices.size() + i2; else	t2.vtxj = i2 - 1;
							if (i3 < 0) t2.vtxk = vertices.size() + i3; else	t2.vtxk = i3 - 1;
							if (j0 < 0) t2.uvi = uvs.size() + j0; else	t2.uvi = j0 - 1;
							if (j2 < 0) t2.uvj = uvs.size() + j2; else	t2.uvj = j2 - 1;
							if (j3 < 0) t2.uvk = uvs.size() + j3; else	t2.uvk = j3 - 1;
							consumedline = consumedline + offset;
							i2 = i3;
							j2 = j3;
							indices.push_back(t2);
						} else {
							nn = sscanf(consumedline, "%u//%u%n", &i3, &k3, &offset);
							if (nn == 2) {
								if (i0 < 0) t2.vtxi = vertices.size() + i0; else	t2.vtxi = i0 - 1;
								if (i2 < 0) t2.vtxj = vertices.size() + i2; else	t2.vtxj = i2 - 1;
								if (i3 < 0) t2.vtxk = vertices.size() + i3; else	t2.vtxk = i3 - 1;
								if (k0 < 0) t2.ni = normals.size() + k0; else	t2.ni = k0 - 1;
								if (k2 < 0) t2.nj = normals.size() + k2; else	t2.nj = k2 - 1;
								if (k3 < 0) t2.nk = normals.size() + k3; else	t2.nk = k3 - 1;								
								consumedline = consumedline + offset;
								i2 = i3;
								k2 = k3;
								indices.push_back(t2);
							} else {
								nn = sscanf(consumedline, "%u%n", &i3, &offset);
								if (nn == 1) {
									if (i0 < 0) t2.vtxi = vertices.size() + i0; else	t2.vtxi = i0 - 1;
									if (i2 < 0) t2.vtxj = vertices.size() + i2; else	t2.vtxj = i2 - 1;
									if (i3 < 0) t2.vtxk = vertices.size() + i3; else	t2.vtxk = i3 - 1;
									consumedline = consumedline + offset;
									i2 = i3;
									indices.push_back(t2);
								} else {
									consumedline = consumedline + 1;
								}
							}
						}
					}
				}

			}

		}
		fclose(f);

		//Calcul de la boite engloblante au chargement de l'objet
		Bbox.maxs_b = vertices[0];
		Bbox.mins_b = vertices[0];
		for (int i =1; i<vertices.size();i++){
			for (int j=0;j < 3;j++){
				Bbox.mins_b[j]= std::min(Bbox.mins_b[j],vertices[i][j]);
				Bbox.maxs_b[j]= std::max(Bbox.maxs_b[j],vertices[i][j]);
			}
		}

	}


	bool intersect(const Ray& r,Vector& P, Vector& N, double &t) const
	{
		//test d'intersection avec la boite
		if (!Bbox.intersection_boite(r)){ //boite non touchée
			return false;
		}

		//boite touchée on teste les triangles
		double t_global_tri = 1E99;
		bool inter_triangle = false;
		for (int i=0;i<indices.size()/3;i++){
			//coords triangle

			int t1 = indices[i*3].vtxi;
			int t2 = indices[i*3].vtxj;
			int t3 = indices[i*3].vtxk;

			Triangle triangle(vertices[t1],vertices[t2],vertices[t3],albedo,objet_miroir,objet_transparent);
			Vector P_triangle;
			Vector N_triangle;
			double t_triangle;
			if(triangle.intersect(r,P_triangle,N_triangle,t_triangle)){
				inter_triangle=true;
				if (t_triangle<t_global_tri){
					t = t_triangle;
					P= P_triangle;
					N=N_triangle;
				}
			}

		}

		return inter_triangle;

			

	}

	std::vector<TriangleIndices> indices;
	std::vector<Vector> vertices;
	std::vector<Vector> normals;
	std::vector<Vector> uvs;
	std::vector<Vector> vertexcolors;

	BoiteEnglobante Bbox;
	
};



// SCENE(OBJET)
class Scene{
	public:
	Scene(){};

	// AJOUT D'UNE SPHERE DANS LA SCENE
	void addSphere(const Sphere& s) {objets.push_back((Objet*)&s);}
	void addTriangle(const Triangle& obj) {objets.push_back((Objet*)&obj);}
	void addMesh(const TriangleMesh& mesh) {objets.push_back((Objet*)&mesh);}

	// ROUTINE D'INTERSECTION ENTRE UN RAYON ET TOUTES LES SPHERES DE LA SCENE
	bool intersect_scene(const Ray& r , Vector& P, Vector& N, int &rang_obj, double &t_min){

		bool inter_scene = false;
		t_min=1E99;

		for (int i=0;i<objets.size();i++){
			Vector P_obj, N_obj;
			double t;
			bool inter_obj = objets[i]->intersect(r,P_obj,N_obj,t);
			if(inter_obj){
				inter_scene=true;
				if(t<t_min){
					t_min=t;
					P = P_obj;
					N= N_obj;
					rang_obj=i;
				}
			}
		}

		return inter_scene;

	}


	// RECUPERE LA COULEUR RENVOYEE PAR UN RAYON
	Vector getColor(Ray &r, int rebond){

	Vector intensite_pixel(0,0,0); //On initialise le pixel en noir
	
	if(rebond<=0){ //Nombre max de rebond dépassé
		return intensite_pixel;
	}

	else{ //Nombre de rebond max non dépassé
		
		Vector P,N;
		int rang_obj;
		double t_min; 
		bool inter = intersect_scene(r,P,N,rang_obj,t_min);

		if (inter){ //le rayon intersecte avec un spherecou un triangle de la scène


			Vector l = (lumiere_position - P);
			l.normalize();

			// LA SPHERE EST MIROIR
			if(objets[rang_obj]->objet_miroir){
				//std::cout << "M";
				Vector vector_reflect = r.u - 2*dot(r.u,N)*N;
				Ray rayon_reflect = Ray(P+N*0.001,vector_reflect);
				return getColor(rayon_reflect,rebond-1);
			}

			else{

				// LA SPHERE EST TRANSPARENTE
				if(objets[rang_obj]->objet_transparent){
					//std::cout << "T";
					double n1 = 1;
					double n2 = 1.3;
					Vector N_transp = N;
					Vector dir = r.u;
					if(dot(r.u,N)>=0){ // Rayon sortant de la sphère, inversion de n1,n2. 
						n1 = 1.3;
						n2 = 1;
						N_transp = Vector(0,0,0)-N;
					}

					double racine = 1-(n1/n2)*(n1/n2)*(1-(dot(N_transp,dir)*(dot(N_transp,dir))));
					if(racine > 0){
						Vector vector_ref = (n1/n2)*(dir - dot(dir,N_transp)*N_transp) - N_transp*sqrt(racine);
						Ray rayon_ref = Ray(P ,vector_ref);
						return getColor(rayon_ref,rebond-1);
						
					}
					
					return Vector(0,0,0);
				}


				else{ //LA SPHERE EST QUELCONQUE
				
				// ECLAIRAGE DIRECT PONCTUEL
				/*
				//On envoie un rayon vers la source de lumière pour voir si on a une ombre (objet de la scène rencontré)
				Vector P_lum,N_lum;
				double t_min_lum;
				int rang_sphere_lum;
				if(intersect_scene(Ray(P+0.001*N,l),P_lum,N_lum,rang_sphere_lum,t_min_lum))
				{
					 
					if (t_min_lum < sqrt((lumiere_position - P).norm2())){ //objet devant la lumière donc ombre
						intensite_pixel = Vector(0,0,0);
						
					}

					else{ //objet derrière = pas d'ombre
						intensite_pixel = objets[rang_obj]->albedo * lumiere_intensite *dot(l,N) / (4 * M_PI * (lumiere_position - P).norm2());
						
					}
				}

				else{ //pas d'objet rencontré = pas d'ombre

					//intensite_pixel = intensite_lumiere * std::max(0.,dot(l,N)) / (4 * M_PI * (lumiere_position - P).norm2());
					intensite_pixel = objets[rang_obj]->albedo * lumiere_intensite * dot(l,N) / (4 * M_PI * (lumiere_position - P).norm2());
				
				}
				
				*/
				// ECLAIRAGE DIRECT DIFFUS
				
				double random_1_direct = uniform(engine);
				double random_2_direct = uniform(engine);

				Vector central_vector = (P-lumiere_position);
				central_vector.normalize();

				//vecteurs aléatoires
				Vector direction_random_local_direct(cos(2*M_PI*random_1_direct)*sqrt(1-random_2_direct),sin(2*M_PI*random_1_direct)*sqrt(1-random_2_direct),sqrt(random_2_direct))  ;
				Vector vector_random_direct(uniform(engine)-0.5,uniform(engine)-0.5,uniform(engine)-0.5);

				//changement de repère
				Vector tangente1_direct = cross(central_vector,vector_random_direct); 
				tangente1_direct.normalize();
				Vector tangente2_direct= cross(tangente1_direct,central_vector);

				Vector direction_random_direct = direction_random_local_direct[2]*central_vector + direction_random_local_direct[0]*tangente1_direct + direction_random_local_direct[1]*tangente2_direct;
				
				Vector point_random_lum = direction_random_direct * lumiere_rayon + lumiere_position;

				Vector w_i = (point_random_lum - P);
				w_i.normalize();

				//Visibilite

				//On envoie un rayon vers la source de lumière pour voir si on a une ombre (objet de la scène rencontré)
				Vector P_lum,N_lum;
				double t_min_lum;
				int rang_sphere_lum;
				int visibilite=1;
				if(intersect_scene(Ray(P+0.001*N,w_i),P_lum,N_lum,rang_sphere_lum,t_min_lum))
				{
					 
					if (t_min_lum < sqrt((lumiere_position - P).norm2())){ //objet devant la lumière donc ombre
						
						visibilite=0;
					}

				}


				intensite_pixel = objets[rang_obj]->albedo *(lumiere_intensite / (4*M_PI*(lumiere_position - P).norm2())*std::max(0.,dot(N,w_i))*std::max(0.,dot(direction_random_direct,(Vector(0,0,0)-w_i))) / dot(central_vector,direction_random_direct)*visibilite);
				
				// ECLAIRAGE INDIRECT
				
				double random_1 = uniform(engine);
				double random_2 = uniform(engine);

				//vecteurs aléatoires
				Vector direction_random_local(cos(2*M_PI*random_1)*sqrt(1-random_2),sin(2*M_PI*random_1)*sqrt(1-random_2),sqrt(random_2))  ;
				Vector vector_random(uniform(engine)-0.5,uniform(engine)-0.5,uniform(engine)-0.5);

				//changement de repère
				Vector tangente1_indirect = cross(N,vector_random); 
				tangente1_indirect.normalize();
				Vector tangente2_indirect= cross(tangente1_indirect,N);

				Vector direction_random_indirect = direction_random_local[2]*N + direction_random_local[0]*tangente1_indirect + direction_random_local[1]*tangente2_indirect;
				Ray rayon_random(P + 0.001*N, direction_random_indirect)	;			

				return intensite_pixel+getColor(rayon_random,rebond-1)*objets[rang_obj]->albedo;
				
				return intensite_pixel;

			}
			}
		
		}

		else{ //Le rayon n'intersecte avec aucune sphère
			//std::cout<<"test";
			return intensite_pixel;
			
		}

	}
};

	std::vector<Objet*> objets;
	Vector lumiere_position;
	double lumiere_intensite;
	double lumiere_rayon;

};





 

// MAIN 
int main() {
	auto start = high_resolution_clock::now();

	// RESOLUTION DE L'IMAGE
	int W = 512;
	int H = 512;

	//Monte Carlo
	const int Nb_rayons=10;

	//NOM DE L'IMAGE
	const char* nom_image = "test_mesh.png";

	// POSITION DE LA CAMERA
	Vector C(0,0,55);

	// CREATION DE LA SCENE

	Sphere s1(Vector(0,0,30),5,Vector(1,0,0),false,false); //albedo vert (couleur de la sphère)
	Sphere s2(Vector(6,10,20),5,Vector(0,1,0),false,false); //albedo rouge (couleur de la sphère)
	Sphere s9(Vector(-7,-6,40),5,Vector(0,0,1),false,false); //albedo rouge (couleur de la sphère)
	Sphere s3(Vector(0,0,-1000),1000-80,Vector(1,0,1),false,false); //mur du fond 
	Sphere s5(Vector(0,-1000,0),1000-10,Vector(1,1,0),false,false); //sol  
	Sphere s6(Vector(1000,0,0),1000-40,Vector(0.7,0.5,0.2),false,false); //mur de droite
	Sphere s7(Vector(-1000,0,0),1000-50,Vector(0,1,1),false,false); //mur de gauche
	Sphere s8(Vector(0,1000,0),1000-50,Vector(1,0,0),false,false); //plafond


	//TRIANGLE:
	Vector T1(-5,-5,30);
	Vector T2(5,-5,30);
	Vector T3(0,5,30);
	Triangle triangle_test(T1,T2,T3,Vector(1,0,0));


	//MAILLAGE


	TriangleMesh mesh_test(Vector(255,255,255),false,false);
	mesh_test.readOBJ("/Users/noelle/Documents/2 - CENTRALE LYON/MOS/MOS 2.2 - Infographie/informatiquegraphique/mesh/dog.obj");


	Scene scene;

	//scene.addSphere(s1);
	//scene.addSphere(s2);
	scene.addSphere(s3);
	//scene.addSphere(s4);
	scene.addSphere(s5);
	scene.addSphere(s6);
	scene.addSphere(s7);
	scene.addSphere(s8);
    //scene.addSphere(s9);

	//scene.addTriangle(triangle_test);

	scene.addMesh(mesh_test);

	// NOMBRE DE REBONDS MAX POUR LES SPHERES MIROIRS
	int nb_rebond_max = 10; //scene.spheres.size();


	// ANGLE DE CHAMPS
	double fov=60*3.14/180;
	double tanfov2 = tan(fov/2);

	
	// LUMIERE POSITION ET INTENSITE ET RAYON POUR LES SOURCES DIFFUSES
	scene.lumiere_position = Vector(-20,40,80);
	scene.lumiere_intensite = 8000000000;
	scene.lumiere_rayon = 7;



	// CALCUL DE L'IMAGE
	std::vector<unsigned char> image(W * H * 3, 0);
	

#pragma omp parallel for
	for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {

			

			//RECUPERATION DE LA COULEUR DU RAYON
			
			Vector pixel(0.,0.,0.);
			for (int k=0; k<Nb_rayons;k++){
				
			//CREATION DU RAYON POUR CHAQUE PIXEL


				//Box-Muller pour antialiasing
				double random_Muller_1 = uniform(engine);
				double random_Muller_2 = uniform(engine);

				double dx = sqrt(-2*log(random_Muller_1))*cos(2*M_PI*random_Muller_2);
				double dy = sqrt(-2*log(random_Muller_1))*sin(2*M_PI*random_Muller_2);


				//MODELE D'OUVERTURE DIAPHRAGME

				double ouverture_diaph=5;

				double distance_map=25;

				double x_diaph = (uniform(engine)-0.5)*2*ouverture_diaph;
				double y_diaph = (uniform(engine)-0.5)*2*ouverture_diaph;


				Vector u(j-W/2+0.5+dx, H-i-H/2+0.5+dy,-W/(2*tanfov2));
				u.normalize();

				Vector plan_net = C + distance_map*u;

				Vector new_C = C+Vector(x_diaph,y_diaph,0);

				Vector new_u = plan_net-new_C;
				new_u.normalize();

				//SANS DIAPH
				Ray r(C,u);

				//AVEC DIAPH
				//Ray r(new_C,new_u);
				
				pixel = pixel + scene.getColor(r,nb_rebond_max)/Nb_rayons;
			}
			//Vector intensite_pixel = scene.getColor(r,nb_rebond_max);

			//STOCKAGE DANS L'IMAGE
			//Correction gamma avec la puissance 1/2.2
			image[(i * W + j) * 3 + 0] = std::min(255.,std::max(0., std::pow(pixel[0],1/2.2)));
			image[(i * W + j) * 3 + 1] = std::min(255.,std::max(0.,std::pow(pixel[1],1/2.2)));
			image[(i * W + j) * 3 + 2] = std::min(255.,std::max(0.,std::pow(pixel[2],1/2.2)));


			
			
			
		}
	}
	stbi_write_png(nom_image, W, H, 3, &image[0], 0);
	auto stop = high_resolution_clock::now();
    auto duration = duration_cast<seconds>(stop - start);
    std::cout << duration.count() << " s"<< std::endl ;

	return 0;
}


//Pour utilisation multithreads compiler avec : g++ main.cpp -lm -fopenmp