#include "Geometry.h"

Geometry::Geometry(const char* obj, double scaling, const Vector& offset, const Vector &couleur, bool mirror, bool transp)
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
            int i0, j1, i2;
            int j0, j1, j2;
            int k0, k1, k2;
            faceGroup.push_back(curGroup);
            int nn = sscanf_s(line, "f %u")

        }
    }


}




////




