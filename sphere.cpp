#pragma once
#include <iostream>
#include "vector.cpp"
#include "ray.cpp"

class Sphere : public Geometry{
public:
    Vector C;
    double R;
    bool interior;

    Sphere(int id, Vector center, double radius, Vector albedo){
        this->C = center;
        this->R = radius;
        this->albedo = albedo;
        this->index = index;
        this->mirror = false;
        this->transparent = false;
        this->interior = false;
        this->n = 1.;
    }
    Sphere(int index, Vector center, double radius, Vector albedo, bool mirror, bool transparent, bool interior, double n){
        this->C = center;
        this->R = radius;
        this->albedo = albedo;
        this->index = index;
        this->mirror = mirror;
        this->transparent = transparent;
        this->interior = interior;
        this->n = n;
    }
    // see if a ray intersects with this sphere
    virtual Intersection intersect(Ray& ray){
        double delta = pow(dot(ray.u, ray.O-C),2) - pow((ray.O-C).norm(),2) + pow(R,2);
        double t;
        if (delta < 0.){ return Intersection(false); }
        else {
            double t_1 = dot(ray.u, C-ray.O) - sqrt(delta);
            double t_2 = dot(ray.u, C-ray.O) + sqrt(delta);
            if (t_2 < 0.){ return Intersection(false); }
            else{
                if (t_1 >= 0.){ t = t_1; }
                else { t = t_2; }
            }
        }
        Vector P = ray.O + t*ray.u;
        if (interior){ return Intersection(P,-(P-C)/(P-C).norm(),true,t,albedo,index); }
        else { return Intersection(P,(P-C).normalise(),true,t,albedo,index); }
    }

    virtual bool checkvis(Ray& ray, double& d){
        double delta = pow(dot(ray.u, ray.O-C),2) - pow((ray.O-C).norm(),2) + pow(R,2);
        double t;
        if (delta < 0.){ return true; }
        else {
            double t_1 = dot(ray.u, C-ray.O) - sqrt(delta);
            double t_2 = dot(ray.u, C-ray.O) + sqrt(delta);
            if (t_2 < 0.){ return true; } 
            else{ 
                if (t_1 >= 0.){ t = t_1; }
                else { t = t_2; }
                if (t < d){ return false; }
                else { return true; }
            }
        }
    }
};