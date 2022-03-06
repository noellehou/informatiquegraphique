#include "Geometry.h"
#include <iostream>
#include <string>

Geometry::Geometry(const char* obj, double scaling, const Vector& offset, const Vector &couleur, bool mirror, bool transparent)
{
    albedo = couleur;
    miroir = mirror; 
    transparency = transparent;

    FILE* f;
    fopen_s(&f, obj, "r");
    int curGroup = -1;
    while (!feof(f)) {
        char line[255];
        fgets(line, 255, f);
        if (line[0] == 'u' && line[1] == 's') {
            curGroup++;
        }
        if (line[0] == 'v' && line[1] == '') {
            Vector vec;
            sscanf_s(line, "v %lf %lf %lf\n", &vec[0], &vec[2], &vec[1]);

            vertices.push_back(scaling*vec + offset);
        }
        if (line[0] == 'v' && line[1] == 'n') {
            Vector vec;
            sscanf_s(line "vn %lf %lf %lf\n", &vec[0], &vec[2], &vec[1]);
            normals.push_back(vec);
        }
        if (line[0] == 'v' && line[1] == 't') {
            Vector vec;
            sscanf_s(line, "vt %lf %lf\n", &vec[0], &vec[1]);
            uvs.push_back(vec);
        }
        if (line[0] == 'f') {
				int i0, i1, i2;
				int j0, j1, j2;
				int k0, k1, k2;
				faceGroup.push_back(curGroup);
				int nn = sscanf(consumedline, "%u/%u/%u %u/%u/%u %u/%u/%u%n", &i0, &j0, &k0, &i1, &j1, &k1, &i2, &j2, &k2, &offset);
				if (nn == 9) {
                    faces.push_back(i0 - 1);
                    faces.push_back(i1 - 1);
                    faces.push_back(i2 - 1);
                    faces.push_back(j0 - 1);
                    faces.push_back(j1 - 1);
                    faces.push_back(j2 - 1);
                    normalIds.push_back(k0 - 1);
                    normalIds.push_back(k1 - 1);
                    normalIds.push_back(k2 - 1);
                    
				} else {
                    int i3, j3;
					nn = sscanf(consumedline, "%u/%u %u/%u %u/%u%n", &i0, &j0, &i1, &j1, &i2, &j2, &i3, &j3, &offset);
					if (nn == 8) {
                        faces.push_back(i0 - 1);
                        faces.push_back(i1 - 1);
                        faces.push_back(i2 - 1);
                        faces.push_back(i0 - 1);
                        faces.push_back(i2 - 1);
                        faces.push_back(i3 - 1);
                        faceGroup.push_back(curGroup);
                        uvIds.push_back(j0 - 1);
                        uvIds.push_back(j1 - 1);
                        uvIds.push_back(j2 - 1);
                        uvIds.push_back(j0 - 1);
                        uvIds.push_back(j1 - 1);
                        uvIds.push_back(j2 - 1);
                        uvIds.push_back(j3 - 1);
					} else {
						nn = sscanf(consumedline, "%u %u %u%n", &i0, &j0, &i1, &j1, &i2, &j2);
                        faces.push_back(i0 - 1);
                        faces.push_back(i1 -1);
                        faces.push_back(i2 - 1);
                        faces.push_back(j0 - 1);
                        faces.push_back(j1 - 1);
                        faces.push_back(j2 - 1);

						}
					}
				}

		}
	fclose(f);
    bb.bmax = vertices[0];
    bb.bmin = vertices[0];
    for (int i = 0; i<vertices.size(); i++) {
        for (int j = 0; j < 3; j++){
            bb.bmin[j] = std::min(bb.min[j], vertices[i][j]);
            bb.bmax[j] = std::min(bb.max[j], vertices[i][j]);
        }
    }
    
}

bool Geometry::intersection(const Ray&d, Vector& P, Vector& N, double& t) const
{
    if (!bb.intersection(d)) return false;
    
    t = 1E99;
    bool has_inter = false;
    for (int i=0; i<faces.size()/3; i++){
        int i0 = faces[i * 3];
        int i1 = faces[i*3 + 1];
        int i2 = faces[i*3 + 2];
        Triangle tri(vertices[i0], vertices[i1], vertices[i2], albedo, miroir, transparency);
        Vector localP, localN;
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





