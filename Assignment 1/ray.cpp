#pragma once
#include <iostream>
#include "vector.cpp"

class Ray {
public:
    Ray(Vector origin, Vector direction, double n1){
        O = origin;
        u = normalize(direction);
        n = n1;
    }
    Vector O, u;
    double n;
};

struct Intersection {
    bool intersect;
    Vector P, N, albedo;
    double t;
    int index;
    Intersection(){ intersect = false; }
    Intersection(bool i){ intersect = i; }
    Intersection(Vector pos, Vector norm, bool intersect, double t, Vector albedo, int index){
        this->intersect = intersect;
        this->P = pos;
        this->N = norm.normalise();
        this->t = t;
        this->albedo = albedo;
        this->index = index;
    }
};

class Geometry {
public:
    Vector albedo; 
    int index;
    bool mirror, transparent;
    double n;
    virtual Intersection intersect(Ray& ray) = 0;
    virtual bool checkvis(Ray& ray, double& d)=0;
};

class BoundingBox {
    public:
    Vector min, max;

    BoundingBox(){}

    bool intersect(Ray& ray, double& distance){
        double u_min1 = (min[0] - ray.O[0]) / ray.u[0];
        double u_min2 = (min[1] - ray.O[1]) / ray.u[1];
        double u_min3 = (min[2] - ray.O[2]) / ray.u[2];
        double u_max1 = (max[0] - ray.O[0]) / ray.u[0];
        double u_max2 = (max[1] - ray.O[1]) / ray.u[1];
        double u_max3 = (max[2] - ray.O[2]) / ray.u[2];
        
        double smallest = std::max(std::min(u_min1,u_max1),std::max(std::min(u_min2,u_max2),std::min(u_min3,u_max3)));
        double biggest = std::min(std::max(u_min1,u_max1),std::min(std::max(u_min2,u_max2),std::max(u_min3,u_max3)));
        if (biggest > smallest && smallest > 0){ 
            distance = smallest;
            return true; 
        }
        return false;
    }
};

class Node {
public:
    Node* child_left;
    Node* child_right;
    BoundingBox bound;
    int starting_triangle, ending_triangle;
    bool is_leaf;
    
    Node(){
        this->child_left = NULL;
        this->child_right = NULL;
    }
};