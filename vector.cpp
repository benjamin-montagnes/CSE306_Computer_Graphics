#pragma once
#include <iostream>
#include <cmath>

class Vector { 
public :
    explicit Vector(double x = 0. , double y = 0. , double z = 0.){ 
        coords[0] = x;
        coords[1] = y;
        coords[2] = z;
    };
    Vector& operator+=(const Vector& b) { 
        coords[0] += b[0]; 
        coords[1] += b[1]; 
        coords[2] += b[2];
        return *this ;
    }
    Vector& normalise(){
        double c = 1 / this->norm();
        coords[0] *= c; 
        coords[1] *= c; 
        coords[2] *= c; 
        return *this;
    }
    Vector& operator/=(const double& c) { 
        double inv = 1/c;
        coords[0] *= inv; 
        coords[1] *= inv; 
        coords[2] *= inv;
        return *this ;   
    }
    const double& operator [ ] (int i) const { return coords [ i ] ; } 
    double& operator []( int i) { return coords[i]; }
    double norm(){ return sqrt(pow(coords[0],2) + pow(coords[1],2) + pow(coords[2],2)) ; }
private:
    double coords [  3 ] ;
};

Vector operator+(const Vector& a, const Vector& b) {
    return Vector(a[0] + b[0], a[1] + b[1], a[2] + b[2]); 
}
Vector operator+(const Vector& a, const double b) {
    return Vector(a[0] + b, a[1] + b, a[2] + b);
}
Vector operator+(const double b, const Vector& a) {
    return Vector(a[0] + b, a[1] + b, a[2] + b);
}
Vector operator-(const Vector& a, const Vector& b) {
    return Vector(a[0] - b[0], a[1] - b[1], a[2] - b[2]); 
}
Vector operator-(const Vector& a) {
    return Vector(-a[0], -a[1], -a[2]); 
}
Vector operator-(const Vector& a, const double b) {
    return Vector(a[0] - b, a[1] - b, a[2] - b);
}
Vector operator*(const Vector& a, const double b) {
    return Vector(a[0]*b, a[1]*b, a[2]*b); 
}
Vector operator*(const double b, const Vector& a) {
    return Vector(a[0]*b, a[1]*b, a[2]*b);
}
Vector operator*(const Vector& a, const Vector& b) {
    return Vector(a[0]*b[0], a[1]*b[1], a[2]*b[2]); 
}
Vector operator/(const Vector& a, const double b) {
    return Vector(a[0]/b, a[1]/b, a[2]/b); 
}
double dot(const Vector& a, const Vector& b) { 
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}
double square(const Vector& a) {
    return dot(a, a);
}
Vector cross(const Vector& a, const Vector& b) {
    return Vector(a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]);
}
double norm(Vector& a) {
    return sqrt(pow(a[0],2) + pow(a[1],2) + pow(a[2],2));
}
Vector normalize(Vector &a){
    return a/a.norm();
}