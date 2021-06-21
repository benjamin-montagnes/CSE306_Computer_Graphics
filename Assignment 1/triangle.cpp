#pragma once
#include <string>
#include <iostream>
#include <stdio.h>
#include <string.h>
#include <algorithm>
#include <vector>
#include "vector.cpp"
#include "ray.cpp"
 
class TriangleIndices {
public:
	TriangleIndices(int vtxi = -1, int vtxj = -1, int vtxk = -1, int ni = -1, int nj = -1, int nk = -1, int uvi = -1, int uvj = -1, int uvk = -1, int group = -1, bool added = false) : vtxi(vtxi), vtxj(vtxj), vtxk(vtxk), uvi(uvi), uvj(uvj), uvk(uvk), ni(ni), nj(nj), nk(nk), group(group) {
	};
	int vtxi, vtxj, vtxk; // indices within the vertex coordinates array
	int uvi, uvj, uvk;  // indices within the uv coordinates array
	int ni, nj, nk;  // indices within the normals array
	int group;       // face group
};

class TriangleMesh : public Geometry {
public:
    ~TriangleMesh() {}
	TriangleMesh(int index) {
        this->albedo = Vector(1.,1.,1.);
        this->mirror = false;
        this->transparent = false;
        this->index = index;
        this->root = new Node();
    }
    TriangleMesh(int index, const char* obj) {
        this->albedo = Vector(1.,1.,1.);
        this->mirror = false;
        this->transparent = false;
        this->index = index;
        this->readOBJ(obj);
        this->root = new Node();
    }
	
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
	}

	BoundingBox compute_bound(int& starting_triangle, int& ending_triangle){
        BoundingBox bound;
        double bmax_x,bmax_y,bmax_z = -std::numeric_limits<double>::max();
        double bmin_x,bmin_y,bmin_z = std::numeric_limits<double>::max();
        for (int i = starting_triangle; i < ending_triangle; i++){
            TriangleIndices T = indices[i];
            std::vector<Vector> points;
            points.push_back(vertices[T.vtxi]);
            points.push_back(vertices[T.vtxj]);
            points.push_back(vertices[T.vtxk]);
            for (int i = 0; i < 3; i++){
                if (points[i][0] > bmax_x){ bmax_x = points[i][0]; }
                if (points[i][0] < bmin_x){ bmin_x = points[i][0]; }
                if (points[i][1] > bmax_y){ bmax_y = points[i][1]; }
                if (points[i][1] < bmin_y){ bmin_y = points[i][1]; }
                if (points[i][2] > bmax_z){ bmax_z = points[i][2]; }
                if (points[i][2] < bmin_z){ bmin_z = points[i][2]; }
            }
        }
        bound.min = Vector(bmin_x,bmin_y,bmin_z);
        bound.max = Vector(bmax_x,bmax_y,bmax_z);
        return bound;
    }

    void bvh(Node* node, int& starting_triangle, int& ending_triangle){
        BoundingBox bound = this->compute_bound(starting_triangle , ending_triangle) ;
        node->bound = bound;
        node->starting_triangle = starting_triangle; 
        node->ending_triangle = ending_triangle;
        Vector diag = node->bound.max - node->bound.min;
        Vector middle_diag = node->bound.min + diag*0.5; 
		int index = 0;
		double max = diag[0];
        if (diag[1] > max){ 
            max = diag[1];
            index = 1;
        }
        if (diag[2] > max){ 
            max = diag[2];
            index = 2; 
        }
        int longest_axis = index;
        int pivot_index = starting_triangle;
        for (int i=starting_triangle ; i<ending_triangle ; i++) {
			TriangleIndices T = indices[i];
			Vector barycenter = (vertices[T.vtxi] + vertices[T.vtxj] + vertices[T.vtxk]) * (0.33333);
            if (barycenter[longest_axis] < middle_diag[longest_axis]) { 
                std::swap(indices[i], indices[pivot_index]); 
                pivot_index++;
            }
        }
        if (pivot_index<=starting_triangle || pivot_index>=ending_triangle-1 || ending_triangle-starting_triangle<5) {
            node->is_leaf = true;
            return ;
        } else { node->is_leaf = false; }

        node->child_left = new Node();
        node->child_right = new Node();
        bvh(node->child_left, starting_triangle, pivot_index); 
        bvh(node->child_right, pivot_index, ending_triangle);
    }

   virtual Intersection intersect(Ray& ray){ 
        double unused_var;
        if(!root->bound.intersect(ray,unused_var)){ return Intersection(false); }
        std::list<Node*> nodes_to_visit; 
        nodes_to_visit.push_front(root);
        double best_inter_distance = std::numeric_limits<double>::max(); 
        Intersection inter;
        Node* curNode;
        while (!nodes_to_visit.empty()){
            curNode = nodes_to_visit.back();
            nodes_to_visit.pop_back();
            if (!curNode->is_leaf){
                double inter_distance;
                if (curNode->child_left->bound.intersect(ray, inter_distance)) {
                    if(inter_distance < best_inter_distance){ nodes_to_visit.push_back(curNode->child_left); }
                }
                if (curNode->child_right->bound.intersect(ray, inter_distance)){
                    if(inter_distance < best_inter_distance){ nodes_to_visit.push_back(curNode->child_right); }
                } 
            } else {
                double min_distance = std::numeric_limits<double>::max();
                Intersection closest_I = Intersection(false);
                for (int i = curNode->starting_triangle; i < curNode->ending_triangle; i++){
                    TriangleIndices T = indices[i];
                    Vector A = vertices[T.vtxi];
                    Vector e1 = vertices[T.vtxj]-A;
                    Vector e2 = vertices[T.vtxk]-A;
                    Vector N = cross(e1,e2);
                    double uN = 1. / dot(ray.u,N);
                    double beta = dot(e2,cross(A-ray.O,ray.u)) * uN;
                    double gamma = - dot(e1,cross(A-ray.O,ray.u)) * uN;
                    double alpha = 1. - beta - gamma;
                    if (alpha > 0. && gamma > 0. && beta > 0.){
                        double t = dot(A-ray.O,N) * uN;
                        if (t < min_distance && t > 0.){ 
                            min_distance = t;
                            closest_I = Intersection(ray.O + t*ray.u, N.normalise(), true, t, this->albedo, index); 
                        }
                    }
                }
                if (min_distance < best_inter_distance){
                    best_inter_distance = min_distance;
                    inter = closest_I;
                }
            }
        }
        return inter;
    }

	std::vector<TriangleIndices> indices;
	std::vector<Vector> vertices;
	std::vector<Vector> normals;
	std::vector<Vector> uvs;
	std::vector<Vector> vertexcolors;
    BoundingBox bound;
    Node* root;

	virtual bool checkvis(Ray& ray, double& d){
        Intersection inter = this->intersect(ray);
        if (inter.intersect){ if (inter.t < d){ return false; } }
        return true;
    }
};