#ifndef DEF_TRIANGLES
#define DEF_TRIANGLES
#pragma once
#include <string.h>
#include <string>
#include <iostream>
#include <stdio.h>
#include <algorithm>
#include <vector>
#include <list>

#include "triangle_indices.h"
#include "../Object/Object.h"
#include "../Ray/Ray.h"
#include "../Vector/Vector.h"




class BBox {
public:
	BBox() {};
	BBox(const Vector& bmin, const Vector& bmax) : bmin(bmin), bmax(bmax) {};
	bool intersection_boite(const Ray& d) const {
		double t_1_x = (bmin[0] - d.origin[0]) / d.direction[0];
		double t_2_x = (bmax[0] - d.origin[0]) / d.direction[0];
		double t_min_x = std::min(t_1_x, t_2_x);
		double t_max_x = std::max(t_1_x, t_2_x);

		double t_1_y = (bmin[1] - d.origin[1]) / d.direction[1];
		double t_2_y = (bmax[1] - d.origin[1]) / d.direction[1];
		double t_min_y = std::min(t_1_y, t_2_y);
		double t_max_y = std::max(t_1_y, t_2_y);

		double t_1_z = (bmin[2] - d.origin[2]) / d.direction[2];
		double t_2_z = (bmax[2] - d.origin[2]) / d.direction[2];
		double t_min_z = std::min(t_1_z, t_2_z);
		double t_max_z = std::max(t_1_z, t_2_z);

		if (std::min(std::min(t_max_x, t_max_y),t_max_z) - std::max(std::max(t_min_x,t_min_y),t_min_z) > 0) return true;
		return false;
	}
	Vector bmin;
	Vector bmax;
};

class Triangle : public Object {
 public:
     Triangle(const Vector& A, const Vector &B, const Vector& C, const Vector &couleur, bool mirror = false, bool transparent = false): A(A), B(B), C(C) {
         albedo = couleur;
         miroir = mirror;
         transparency = transparent;

     };

     bool intersection(const Ray& d, Vector& P, Vector& N, double &t) const{

         N = 0-cross(B-A,C-A).getNormalized();
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

		N.normalize();
         return true;
     }
     const Vector &A, &B, &C;
 };

class TriangleMesh :  Object
{
public:
  ~TriangleMesh(){}
	TriangleMesh(const Vector &couleur, bool mirror = false, bool transparent = false){
	albedo = couleur;
    miroir = mirror; 
    transparency = transparent;

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

	bb.bmax = vertices[0];
	bb.bmin = vertices[0];
	for (int i=1; i<vertices.size();i++){
		for (int j=0;j<3;j++){
			bb.bmin[j] = std::min(bb.bmin[j], vertices[i][j]);
			bb.bmax[j] = std::max(bb.bmax[j], vertices[i][j]);
		}
	}
	}

	bool intersection(const Ray& d, Vector& P, Vector& N, double &t) const {
		if (!bb.intersection_boite(d)) return false;

    	t = 1E99;
    	bool has_inter = false;
    	for (int i=0; i < indices.size()/3; i++){
        	int i0 = indices[i * 3].vtxi;
        	int i1 = indices[i*3 + 1].vtxj;
        	int i2 = indices[i*3 + 2].vtxk;
        	Triangle tri(vertices[i0], vertices[i1], vertices[i2], albedo, miroir, transparency);
        	Vector localP;
			Vector localN;
        	double localt;
        if (tri.intersection(d, localP, localN, localt)) {
            has_inter = true;
            if (localt < t){
                 t = localt;
                 P = localP;
                 N = localN;
            }
        }
    }
    return has_inter;
	}


	std::vector<TriangleIndices> indices;
	std::vector<Vector> vertices;
	std::vector<Vector> normals;
	std::vector<Vector> uvs;
	std::vector<Vector> vertexcolors;
	
	

	BBox bb;
};


#endif