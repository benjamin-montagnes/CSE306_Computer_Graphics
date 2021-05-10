#include <iostream>
#include <chrono>
#include <list>
#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb/stb_image_write.h"
#define STB_IMAGE_IMPLEMENTATION
#include "stb/stb_image.h"

#include "vector.cpp"
#include "sphere.cpp"
#include "triangle.cpp"
#include "scene.cpp"


int main() {
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

    int W = 512;
    int H = 512;
    std::vector<unsigned char> image(W*H * 3, 0);

    // first define the scene , variables , 4 spheres on the side...
    std::vector<Geometry*> arrayS;
    arrayS.push_back(new Sphere(0, Vector(0.,1000.,0.),   940., Vector(1.,0.,0.)));
    arrayS.push_back(new Sphere(1, Vector(0.,-1000.,0.),  990., Vector(0.,0.,1.)));
    arrayS.push_back(new Sphere(2, Vector(0.,0., -1000.), 940., Vector(0.,1.,0.)));
    arrayS.push_back(new Sphere(3, Vector(0.,0.,1000.),   940., Vector(1.,0.,1.)));
    arrayS.push_back(new Sphere(4, Vector(-1000.,0.,0.),  940., Vector(0.,1.,1.)));
    arrayS.push_back(new Sphere(5, Vector(1000.,0.,0.),   940., Vector(1.,1.,0.)));

    //For the white sphere :
    // arrayS.push_back(new Sphere(6, Vector(0.,0.,0.), 10., Vector(1.,1.,1.)));

    //For the 3 spheres :
    // arrayS.push_back(new Sphere(6, Vector(-20.,0.,0.), 10.,  Vector(0.,0.,0.),  true,  false, false, 1.));//mirror
    // arrayS.push_back(new Sphere(7, Vector(0.,0.,0.),   10.,  Vector(0.,0.,0.),  false, true,  false, 1.5));//glass
    // arrayS.push_back(new Sphere(8,Vector(20.,0.,0.),   10.,  Vector(0.,0.,0.),  false, true,  false, 1.5));//hollow
    // arrayS.push_back(new Sphere(9,Vector(20.,0.,0.),   9.5,  Vector(0.,0.,0.),  false, true,  true,  1.5));//hollow

    //For the cat :
    TriangleMesh* triangleMesh = new TriangleMesh(6, "model/cat.obj");
    for (int i = 0; i < (triangleMesh->vertices).size(); i++){
        (triangleMesh->vertices)[i] = 0.6*(triangleMesh->vertices)[i] + Vector(0.,-10.,0.);
    }
    int s = 0;
    int e = (triangleMesh->indices).size();
    BoundingBox bound = triangleMesh->compute_bound(s, e);
    triangleMesh->bound = bound;
    triangleMesh->bvh(triangleMesh->root, s, e);
    arrayS.push_back(triangleMesh);

    //We set up the camera and the scene
    int NB_PATHS = 64;
    int max_path_length = 5;
    double gamma = 1./2.2;
    Vector camera = Vector(0.,0.,55.);
    Scene scene = Scene(arrayS, camera, Vector(-10.,20.,40.), 1., 2e10);

    // Parallelization
    #pragma omp parallel for schedule(dynamic,1) 
    // then scan all pixels
    for (int i = 0; i < H; i++) {
        for (int j = 0; j < W; j++) {
            Vector color;
            for (int k = 0; k < NB_PATHS; k++){
                //we do the antialias
                double x,y;
                boxMuller(0.4,x,y);
                Vector pixel = Vector(camera[0] + j + 0.5 - W*0.5, camera[1] + H*0.5 - i - 0.5, camera[2]-(0.5 * W / tan(M_PI / 6.)));
                pixel[0] += x;
                pixel[1] += y;
                Ray ray = Ray(camera, (pixel-camera).normalise(), 1.); // cast a ray from the camera center to pixel i , j
                color += scene.getColor(ray, max_path_length); //stores color channel
            }
            // we do the color correction with gamma
            image[(i*W + j) * 3 + 0] = std::min(255., pow(color[0]/NB_PATHS, gamma)); 
            image[(i*W + j) * 3 + 1] = std::min(255., pow(color[1]/NB_PATHS, gamma)); 
            image[(i*W + j) * 3 + 2] = std::min(255., pow(color[2]/NB_PATHS, gamma)); 
        }
    }
    // save image and return 0
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    stbi_write_png("image.png", W, H, 3, &image[0], 0);
    std::cout << "Success! Saved in 'image.png'" << std::endl;
    std::cout << "Total time : " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << " ms" << std::endl;
    return 0;
}