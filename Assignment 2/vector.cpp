#pragma once
#include <iostream>
#include <math.h>
#include <vector>
#include <limits>
#include <random>



class Vector{
    public:
        explicit Vector(double x = 0., double y = 0., double z = 0.){
            coords[0] = x;
            coords[1] = y;
            coords[2] = z;
        };
        Vector& operator+=(const Vector& b){
            coords[0] += b[0];
            coords[1] += b[1];
            coords[2] += b[2];
            return *this;
        }
        Vector& operator*=(const double t){
            coords[0] *= t;
            coords[1] *= t;
            coords[2] *= t;
            return *this;
        }
        const double &operator[](int i) const { return coords[i]; }
        double &operator[](int i) { return coords[i]; }
    private:
        double coords[3];
};

bool operator==(const Vector &a, const Vector &b){
    bool temp = ((a[0] == b[0]) && (a[1] == b[1]) && (a[2] == b[2])) ? true : false;
    return temp;
}

bool operator!=(const Vector &a, const Vector &b){
    bool temp = ((a[0] != b[0]) || (a[1] != b[1]) || (a[2] != b[2])) ? true : false;
    return temp;
}

Vector operator+(const Vector &a, const Vector &b){
    return Vector(a[0] + b[0], a[1] + b[1], a[2] + b[2]);
}

Vector operator-(const Vector &a, const Vector &b){
    return Vector(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
}

Vector operator-(const Vector &b){
    return Vector(- b[0], - b[1], - b[2]);
}

double dot(const Vector &a, const Vector &b){
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

double norm(const Vector &a){
    return sqrt(pow(a[0], 2) + pow(a[1], 2) + pow(a[2], 2));
}

double norm_square(const Vector &a){
    return dot(a,a);
}

Vector cross(const Vector &a, const Vector &b){
    return Vector(a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]);
}

Vector normalization(const Vector &a){
    if (norm_square(a) == 1.) return a;
    double nrm = norm(a);
    return Vector(a[0] / nrm, a[1] / nrm, a[2] / nrm);
}

Vector operator*(const double t, const Vector &a){
    return Vector(t * a[0], t * a[1], t * a[2]);
}

Vector operator*(const Vector &a, const double t){
    return Vector(t * a[0], t * a[1], t * a[2]);
}

Vector operator*(const Vector &a, const Vector &b){
    return Vector(a[0] * b[0], a[1] * b[1], a[2] * b[2]);
}

Vector operator/(const Vector &a, const double t){
    return Vector(a[0]/t, a[1]/t, a[2]/t);
}

Vector operator/(const double t, const Vector &a){
    return Vector(a[0] / t, a[1] / t, a[2] / t);
}

Vector perpen(const Vector &a){
    return Vector(a[1],-a[0],a[2]);
}

double distance_square(const Vector &a, const Vector &b){
    return pow(a[0] - b[0], 2) + pow(a[1] - b[1], 2) + pow(a[2] - b[2], 2);
}

#define PI 3.14159265
#define rtd double(180 / PI)
#define dtr double(PI / 180)

double inf = std::numeric_limits<double>::infinity();

static std::default_random_engine engine_scene;
static std::uniform_real_distribution<double> uniform_scene(0, 1);

std::vector<Vector> random_points(int n){
    double u,v;
    std::vector<Vector> res;
    for (int i = 0; i < n; i++){
        u = uniform_scene(engine_scene);
        v = uniform_scene(engine_scene);
        res.push_back(Vector(u,v,0));
    }
    return res;
}
