#include "triangles.h"



// TriangleMesh::TriangleMesh(const char* obj, double scaling, const Vector& offset, const Vector& couleur, bool mirror, bool transparent) {
//     albedo = couleur;
//     miroir = mirror; 
//     transparency = transparent;

// 	char matfile[255];
// 	char grp[255];

// 	FILE* f;
// 	f = fopen(obj, "r");
// 	int curGroup = -1;
// 	while (!feof(f)) {
// 		char line[255];
// 		if (!fgets(line, 255, f)) break;

// 		std::string linetrim(line);
// 		linetrim.erase(linetrim.find_last_not_of(" \r\t") + 1);
// 		strcpy(line, linetrim.c_str());

// 		if (line[0] == 'u' && line[1] == 's') {
// 			sscanf(line, "usemtl %[^\n]\n", grp);
// 			curGroup++;
// 		}

// 		if (line[0] == 'v' && line[1] == ' ') {
// 			Vector vec;

// 			Vector col;
// 			if (sscanf(line, "v %lf %lf %lf %lf %lf %lf\n", &vec[0], &vec[1], &vec[2], &col[0], &col[1], &col[2]) == 6) {
// 				col[0] = std::min(1., std::max(0., col[0]));
// 				col[1] = std::min(1., std::max(0., col[1]));
// 				col[2] = std::min(1., std::max(0., col[2]));

// 				vertices.push_back(vec);
// 				vertexcolors.push_back(col);

// 			} else {
// 				sscanf(line, "v %lf %lf %lf\n", &vec[0], &vec[1], &vec[2]);
// 				vertices.push_back(vec);
// 			}
// 		}
// 		if (line[0] == 'v' && line[1] == 'n') {
// 			Vector vec;
// 			sscanf(line, "vn %lf %lf %lf\n", &vec[0], &vec[1], &vec[2]);
// 			normals.push_back(vec);
// 		}
// 		if (line[0] == 'v' && line[1] == 't') {
// 			Vector vec;
// 			sscanf(line, "vt %lf %lf\n", &vec[0], &vec[1]);
// 			uvs.push_back(vec);
// 		}
// 		if (line[0] == 'f') {
// 			TriangleIndices t;
// 			int i0, i1, i2, i3;
// 			int j0, j1, j2, j3;
// 			int k0, k1, k2, k3;
// 			int nn;
// 			t.group = curGroup;

// 			char* consumedline = line + 1;
// 			int offset;

// 			nn = sscanf(consumedline, "%u/%u/%u %u/%u/%u %u/%u/%u%n", &i0, &j0, &k0, &i1, &j1, &k1, &i2, &j2, &k2, &offset);
// 			if (nn == 9) {
// 				if (i0 < 0) t.vtxi = vertices.size() + i0; else	t.vtxi = i0 - 1;
// 				if (i1 < 0) t.vtxj = vertices.size() + i1; else	t.vtxj = i1 - 1;
// 				if (i2 < 0) t.vtxk = vertices.size() + i2; else	t.vtxk = i2 - 1;
// 				if (j0 < 0) t.uvi = uvs.size() + j0; else	t.uvi = j0 - 1;
// 				if (j1 < 0) t.uvj = uvs.size() + j1; else	t.uvj = j1 - 1;
// 				if (j2 < 0) t.uvk = uvs.size() + j2; else	t.uvk = j2 - 1;
// 				if (k0 < 0) t.ni = normals.size() + k0; else	t.ni = k0 - 1;
// 				if (k1 < 0) t.nj = normals.size() + k1; else	t.nj = k1 - 1;
// 				if (k2 < 0) t.nk = normals.size() + k2; else	t.nk = k2 - 1;
// 				indices.push_back(t);
// 			} else {
// 				nn = sscanf(consumedline, "%u/%u %u/%u %u/%u%n", &i0, &j0, &i1, &j1, &i2, &j2, &offset);
// 				if (nn == 6) {
// 					if (i0 < 0) t.vtxi = vertices.size() + i0; else	t.vtxi = i0 - 1;
// 					if (i1 < 0) t.vtxj = vertices.size() + i1; else	t.vtxj = i1 - 1;
// 					if (i2 < 0) t.vtxk = vertices.size() + i2; else	t.vtxk = i2 - 1;
// 					if (j0 < 0) t.uvi = uvs.size() + j0; else	t.uvi = j0 - 1;
// 					if (j1 < 0) t.uvj = uvs.size() + j1; else	t.uvj = j1 - 1;
// 					if (j2 < 0) t.uvk = uvs.size() + j2; else	t.uvk = j2 - 1;
// 					indices.push_back(t);
// 				} else {
// 					nn = sscanf(consumedline, "%u %u %u%n", &i0, &i1, &i2, &offset);
// 					if (nn == 3) {
// 						if (i0 < 0) t.vtxi = vertices.size() + i0; else	t.vtxi = i0 - 1;
// 						if (i1 < 0) t.vtxj = vertices.size() + i1; else	t.vtxj = i1 - 1;
// 						if (i2 < 0) t.vtxk = vertices.size() + i2; else	t.vtxk = i2 - 1;
// 						indices.push_back(t);
// 					} else {
// 						nn = sscanf(consumedline, "%u//%u %u//%u %u//%u%n", &i0, &k0, &i1, &k1, &i2, &k2, &offset);
// 						if (i0 < 0) t.vtxi = vertices.size() + i0; else	t.vtxi = i0 - 1;
// 						if (i1 < 0) t.vtxj = vertices.size() + i1; else	t.vtxj = i1 - 1;
// 						if (i2 < 0) t.vtxk = vertices.size() + i2; else	t.vtxk = i2 - 1;
// 						if (k0 < 0) t.ni = normals.size() + k0; else	t.ni = k0 - 1;
// 						if (k1 < 0) t.nj = normals.size() + k1; else	t.nj = k1 - 1;
// 						if (k2 < 0) t.nk = normals.size() + k2; else	t.nk = k2 - 1;
// 						indices.push_back(t);
// 					}
// 				}
// 			}

// 			consumedline = consumedline + offset;

// 			while (true) {
// 				if (consumedline[0] == '\n') break;
// 				if (consumedline[0] == '\0') break;
// 				nn = sscanf(consumedline, "%u/%u/%u%n", &i3, &j3, &k3, &offset);
// 				TriangleIndices t2;
// 				t2.group = curGroup;
// 				if (nn == 3) {
// 					if (i0 < 0) t2.vtxi = vertices.size() + i0; else	t2.vtxi = i0 - 1;
// 					if (i2 < 0) t2.vtxj = vertices.size() + i2; else	t2.vtxj = i2 - 1;
// 					if (i3 < 0) t2.vtxk = vertices.size() + i3; else	t2.vtxk = i3 - 1;
// 					if (j0 < 0) t2.uvi = uvs.size() + j0; else	t2.uvi = j0 - 1;
// 					if (j2 < 0) t2.uvj = uvs.size() + j2; else	t2.uvj = j2 - 1;
// 					if (j3 < 0) t2.uvk = uvs.size() + j3; else	t2.uvk = j3 - 1;
// 					if (k0 < 0) t2.ni = normals.size() + k0; else	t2.ni = k0 - 1;
// 					if (k2 < 0) t2.nj = normals.size() + k2; else	t2.nj = k2 - 1;
// 					if (k3 < 0) t2.nk = normals.size() + k3; else	t2.nk = k3 - 1;
// 					indices.push_back(t2);
// 					consumedline = consumedline + offset;
// 					i2 = i3;
// 					j2 = j3;
// 					k2 = k3;
// 				} else {
// 					nn = sscanf(consumedline, "%u/%u%n", &i3, &j3, &offset);
// 					if (nn == 2) {
// 						if (i0 < 0) t2.vtxi = vertices.size() + i0; else	t2.vtxi = i0 - 1;
// 						if (i2 < 0) t2.vtxj = vertices.size() + i2; else	t2.vtxj = i2 - 1;
// 						if (i3 < 0) t2.vtxk = vertices.size() + i3; else	t2.vtxk = i3 - 1;
// 						if (j0 < 0) t2.uvi = uvs.size() + j0; else	t2.uvi = j0 - 1;
// 						if (j2 < 0) t2.uvj = uvs.size() + j2; else	t2.uvj = j2 - 1;
// 						if (j3 < 0) t2.uvk = uvs.size() + j3; else	t2.uvk = j3 - 1;
// 						consumedline = consumedline + offset;
// 						i2 = i3;
// 						j2 = j3;
// 						indices.push_back(t2);
// 					} else {
// 						nn = sscanf(consumedline, "%u//%u%n", &i3, &k3, &offset);
// 						if (nn == 2) {
// 							if (i0 < 0) t2.vtxi = vertices.size() + i0; else	t2.vtxi = i0 - 1;
// 							if (i2 < 0) t2.vtxj = vertices.size() + i2; else	t2.vtxj = i2 - 1;
// 							if (i3 < 0) t2.vtxk = vertices.size() + i3; else	t2.vtxk = i3 - 1;
// 							if (k0 < 0) t2.ni = normals.size() + k0; else	t2.ni = k0 - 1;
// 							if (k2 < 0) t2.nj = normals.size() + k2; else	t2.nj = k2 - 1;
// 							if (k3 < 0) t2.nk = normals.size() + k3; else	t2.nk = k3 - 1;								
// 							consumedline = consumedline + offset;
// 							i2 = i3;
// 							k2 = k3;
// 							indices.push_back(t2);
// 						} else {
// 							nn = sscanf(consumedline, "%u%n", &i3, &offset);
// 							if (nn == 1) {
// 								if (i0 < 0) t2.vtxi = vertices.size() + i0; else	t2.vtxi = i0 - 1;
// 								if (i2 < 0) t2.vtxj = vertices.size() + i2; else	t2.vtxj = i2 - 1;
// 								if (i3 < 0) t2.vtxk = vertices.size() + i3; else	t2.vtxk = i3 - 1;
// 								consumedline = consumedline + offset;
// 								i2 = i3;
// 								indices.push_back(t2);
// 							} else {
// 								consumedline = consumedline + 1;
// 							}
// 						}
// 					}
// 				}
// 			}

// 		}

// 	}
// 	fclose(f);

//     bb.bmax = vertices[0];
//     bb.bmin = vertices[0];
//     for (int i = 1; i<vertices.size(); i++) {
//         for (int j = 0; j < 3; j++){
//             bb.bmin[j] = std::min(bb.bmin[j], vertices[i][j]);
//             bb.bmax[j] = std::max(bb.bmax[j], vertices[i][j]);
//         }
//     }
    
// }

// bool TriangleMesh::intersection(const Ray& d, Vector& P, Vector& N, double &t) const
// {
// 	if (!bb.intersection(d)) return false;

//     t = 1E99;
//     bool has_inter = false;
//     for (int i=0; i < faces.size()/3; i++){
//         int i0 = faces[i * 3];
//         int i1 = faces[i*3 + 1];
//         int i2 = faces[i*3 + 2];
//         Triangle tri(vertices[i0], vertices[i1], vertices[i2], albedo, miroir, transparency);
//         Vector localP, localN;
//         double localt;
//         if (tri.intersection(d, localP, localN, localt)) {
//             has_inter = true;
//             if (localt < t){
//                  t = localt;
//                  P = localP;
//                  N = localN;
//             }
//         }
//     }
//     return has_inter;
// }