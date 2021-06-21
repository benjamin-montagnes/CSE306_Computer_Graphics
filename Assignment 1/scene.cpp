#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include "vector.cpp"
#include "ray.cpp"

static std::default_random_engine engine(10) ; // random seed = 10 
static std::uniform_real_distribution<double> uniform(0, 1);

void boxMuller(double stdev , double &x, double &y) {
    double r1 = uniform(engine);
    double r2 = uniform(engine);
    x = sqrt(-2. * log(r1))*cos(2.*M_PI*r2)*stdev;
    y = sqrt(-2. * log(r1))*sin(2.*M_PI*r2)*stdev;
}

Vector random_cos(const Vector& N){
    double r1 = uniform(engine);
    double r2 = uniform(engine);

    double x = cos(2*M_PI*r1) * sqrt(1-r2);
    double y = sin(2*M_PI*r1) * sqrt(1-r2);
    double z = sqrt(r2);

    Vector T1;
    double min_component = N[0];
    double min_index = 0;
    for (int i=1;i<3;i++){
        if (N[i] < min_component){ 
            min_component = N[i];
            min_index = i;
        }
    }
    if (min_index == 0){ T1 = Vector(0.,N[2],-N[1]); }
    else if (min_index == 1){ T1 = Vector(N[2],0.,-N[0]); }
    else { T1 = Vector(N[1],-N[0],0.); }

    T1 = normalize(T1);

    Vector T2 = cross(N,T1);
    
    return (x*T1 + y*T2 + z*N).normalise();
}

//An array/std::vector of Spheres
class Scene {
public:
    Vector cam;
    Vector light_source;
    std::vector<Geometry*> arrayG;
    double n;
    double light_intensity;

    Scene(std::vector<Geometry*> arrayG, Vector cam, Vector light_source, double n, double light_intensity){
        this->arrayG = arrayG;
        this->cam = cam;
        this->light_source = light_source;
        this->n = n;
        this->light_intensity = light_intensity;
    }
    // Computes the point of intersection between a Ray and the sphere, if any
    Intersection intersection_spheres(Ray& ray){
        double dist = std::numeric_limits<double>::max();
        Intersection best;
        best.intersect = false;
        Geometry* geo_index;
        for(std::vector<Geometry*>::iterator i = arrayG.begin(); i != arrayG.end(); ++i) { 
            geo_index = *i;
            Intersection inter = geo_index->intersect(ray);
            if (inter.intersect && inter.t < dist){
                dist = inter.t;
                best = inter;
            }
        }
        return best;
    }

    double visibility(Ray& ray, double& d){
        for (int i = 0; i < arrayG.size(); i++){
            bool vi = arrayG[i]->checkvis(ray, d);
            if(!vi){ return 0.; }
        }
        return 1.0;
    }

    Vector getColor(Ray& ray, const int& ray_depth){
        if (ray_depth < 0.){ return Vector(0.,0.,0.); } // terminates recursion at some point
        Intersection inter = this -> intersection_spheres(ray);
        if (inter.intersect){ 
            double e = 0.0001;
            if (arrayG[inter.index]->transparent){ 
                //handle refractive surfaces
                double n1 = ray.n;
                double n2;
                Vector N;
                if ( dot(ray.u,inter.N) > 0. ){
                    N = -inter.N;
                    n2 = this->n;
                } else {
                    N = inter.N;
                    n2 = arrayG[inter.index]->n;
                }
                double k0 = (n1-n2)*(n1-n2) / ((n1+n2)*(n1+n2));
                double u = (double) rand()/RAND_MAX ;
                if (u < k0 + (1 - k0) * pow(1 - std::abs(dot(N,ray.u)), 5.)){ // reflect
                    Ray reflected_ray = Ray(inter.P+e*inter.N, ray.u - ( 2 * dot(ray.u,inter.N) * inter.N ), ray.n);
                    return getColor(reflected_ray, ray_depth-1);
                }
                double ratio = n1 / n2;
                double dot_uN = dot(ray.u,N);
                double x = 1 - ratio*ratio * (1 - dot_uN*dot_uN);   
                if ( x < 0. ){
                    Ray int_reflec_ray = Ray(inter.P+e*inter.N, ray.u - ( 2 * dot(ray.u,inter.N) * inter.N ), n2);
                    return getColor(int_reflec_ray, ray_depth-1);
                } else {
                    Vector wT = ratio * (ray.u - dot(ray.u, N) * N);       
                    Vector wN = - N * sqrt(x); 
                    Vector w = wT + wN;
                    Ray refracted_ray = Ray(inter.P - e*N, w, n2); 
                    return getColor(refracted_ray, ray_depth-1);
                }
            } else if (arrayG[inter.index]->mirror){
                //handle mirror surfaces
                Ray reflected_ray = Ray(inter.P+e*inter.N, ray.u - ( 2 * dot(ray.u,inter.N) * inter.N ), ray.n);
                return getColor(reflected_ray, ray_depth-1);
            } else {
                // handle diffuse surfaces
                double distance_light = (light_source - inter.P).norm();
                Vector light_direction = (light_source - inter.P)/distance_light;
                Ray direct_ray = Ray(inter.P + e*inter.N, light_direction, ray.n);
                // add direct lighting
                double visibility = this->visibility(direct_ray, distance_light); // computes the visibility term by launching a ray towards the light source
                Vector Lo = light_intensity/(4.*M_PI*pow( distance_light , 2)) * inter.albedo/M_PI * visibility * std::max(dot(inter.N,light_direction),0.);
                
                // add indirect lighting
                Ray randomRay = Ray(inter.P + e*inter.N, random_cos(inter.N), ray.n); // randomly sample ray using random_cos
                Lo += inter.albedo  * getColor(randomRay, ray_depth-1);
                
                return Lo;
            }  
        } else { return Vector(0.,0.,0.); }
    }
};
