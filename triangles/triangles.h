#ifndef DEF_TRIANGLES
#define DEF_TRIANGLES

#include <string.h>
#include <string>
#include <iostream>
#include <stdio.h>
#include <algorithm>
#include <vector>
#include <list>

#include "triangle_indices.h"
#include "../Object/Object.h"



class BBox
{
public:
	Vector m, M;
	BBox();
	BBox(const Vector&, const Vector&);
	~BBox();
};


class Node
{
public:
	BBox box;

	Node* left;
	Node* right;

	unsigned int i1;
	unsigned int i2;
	
	Node();
	~Node();
};



class TriangleMesh : public Object
{
public:
	TriangleMesh();
	~TriangleMesh();

	Node root_box;
	
	bool intersect_tri(const TriangleIndices&, const Ray&, Vector&, Vector&, double&);
	bool intersect_BBox(const Ray&, const BBox&);

	virtual bool intersect(const Ray&);
    virtual bool intersect(const Ray&, Vector&, Vector&);
    virtual bool intersect(const Ray&, Vector&, Vector&, double&);

	// double center_of_triangle(unsigned int, unsigned int);

	void readOBJ(const char* obj);
	void create_BBox();
	BBox create_bounding_box(int, int);

	std::vector<TriangleIndices> indices;
	std::vector<Vector> vertices;
	std::vector<Vector> normals;
	std::vector<Vector> uvs;
	std::vector<Vector> vertexcolors;
	
};




#endif