#include <iostream>
#include <chrono>
#include <map>

#include "vector.cpp"

class Polygon{
    public:
        std::vector<Vector> vertices;
        std::vector<std::pair<Vector,Vector>> edges;
    explicit Polygon(std::vector<Vector> V = {}){
        vertices = V;
        for (size_t i = 0; i < vertices.size(); i++){
            std::pair<Vector, Vector> edge(V[i], V[(i != 0) ? i - 1 : vertices.size() - 1]);
            edges.push_back(edge);
        }
    }
};

Vector intersect(Vector prevVertex, Vector curVertex, std::pair<Vector,Vector> clipEdge){
    Vector u = clipEdge.first;
    Vector v = clipEdge.second;
    Vector N = Vector(v[1]-u[1],u[0]-v[0],0);
    double t = dot(u-prevVertex,N)/dot(curVertex-prevVertex,N);
    if (t < 0 || t > 1) return Vector(inf,inf,inf);
    return prevVertex + t*(curVertex-prevVertex);
};

Vector intersect_voronoi(Vector &A, Vector &B, Vector &Pi, Vector &Pj, Vector M){
    double t = dot(M-A,Pi-Pj)/dot(B-A,Pi-Pj);
    return A + t*(B-A);
}

bool inside_voronoi(Vector X, Vector Pi, Vector Pj, Vector M){return dot(X-M,Pj-Pi)<0;}
bool inside_voronoi(Vector P, std::pair<Vector,Vector> clipEdge){
    Vector u = clipEdge.first;
    Vector v = clipEdge.second;
    Vector N = Vector(v[1] - u[1], u[0] - v[0],0);
    return dot(P-u,N) <= 0 ? true : false;
}

double max_dist(Polygon polygon,Vector P){
    double dist = 0;
    for (auto &vert: polygon.vertices) dist = std::max(dist, distance_square(vert, P));
    return dist;
}

Polygon clip_poly(Polygon subjectPolygon, Polygon clipPolygon){
    Polygon outPolygon;    
    for(auto &clipEdge : clipPolygon.edges){
        outPolygon = Polygon();
        for (size_t i = 0; i < subjectPolygon.vertices.size(); i++){
            Vector curVertex = subjectPolygon.vertices[i];
            Vector prevVertex = subjectPolygon.vertices[(i>0)?(i-1):subjectPolygon.vertices.size()-1];
            Vector intersection = intersect(prevVertex,curVertex,clipEdge);
            if (inside_voronoi(curVertex,clipEdge)){
                if (!inside_voronoi(prevVertex,clipEdge)) outPolygon.vertices.push_back(intersection);
                outPolygon.vertices.push_back(curVertex);
            }
            else if (inside_voronoi(prevVertex,clipEdge)) outPolygon.vertices.push_back(intersection);
        }
        subjectPolygon = outPolygon;
    }
    return outPolygon;
};

std::vector<Polygon> create_voronoi(std::vector<Vector> points, std::vector<double> weights = {}, Polygon space = Polygon({Vector(0, 0, 0), Vector(0, 1, 0), Vector(1, 1, 0), Vector(1, 0, 0)}));
std::vector<Polygon> create_voronoi(std::vector<Vector> points, std::vector<double> weights, Polygon space){
    std::vector<Polygon> subjectPolygon;
    Polygon outPolygon;
    if (weights.size() == 0){
        for (int i = 0; i < points.size(); i++) weights.push_back(0.1);
    }
    for (int i = 0; i < points.size();i++){
        auto curPoint = points[i]; 
        std::map<double,Vector> otherPoints;
        std::vector<Vector> otherPoints_vector;
        for(auto &tempPoint: points){
            if (tempPoint != curPoint){
                otherPoints[distance_square(curPoint,tempPoint)] = tempPoint;
                otherPoints_vector.push_back(tempPoint);
            }
        }
        outPolygon = space;
        for (int j = 0; j < otherPoints_vector.size();j++){
            Vector tempPoint = otherPoints_vector[j];
            size_t n = outPolygon.vertices.size();
            Polygon tempPolygon = Polygon();
            for(size_t t = 0; t < n; t++){
                Vector curVertex = outPolygon.vertices[t];
                Vector prevVertex = outPolygon.vertices[(t>0)?(t-1):n-1];
                Vector M_temp = (curPoint + tempPoint)/2;
                Vector M = M_temp + (weights[i] - weights[j])/(2*norm_square(curPoint - tempPoint)) * (tempPoint - curPoint);
                Vector intersection = intersect_voronoi(prevVertex,curVertex,curPoint,tempPoint,M);
                if(inside_voronoi(curVertex,curPoint,tempPoint,M)){
                    if (!inside_voronoi(prevVertex, curPoint,tempPoint, M)) tempPolygon.vertices.push_back(intersection);
                    tempPolygon.vertices.push_back(curVertex);
                }
                else if (inside_voronoi(prevVertex, curPoint,tempPoint, M)) tempPolygon.vertices.push_back(intersection);
            }
            outPolygon = tempPolygon;
        }
        subjectPolygon.push_back(outPolygon);
    }
    return subjectPolygon;
}

std::vector<Vector> tesselation(std::vector<Vector> points, int iter = 20){
    if (iter == 0) return points;
    std::vector<Polygon> vor = create_voronoi(points);
    std::vector<Vector> new_points;
    for (auto &cell: vor){
        int n = cell.vertices.size();
        double A = 0;   
        Vector C = Vector(); 
        for (size_t i = 0; i < n;i++){
            Vector cur = cell.vertices[i];
            auto next = cell.vertices[(i == n-1) ? 0 : (i+1)];
            A += (cur[0]*next[1] - next[0]*cur[1]);
        }
        A /= 2;
        for (size_t i = 0; i < n;i++){
            auto cur = cell.vertices[i];
            auto next = cell.vertices[(i == n-1) ? 0 : (i+1)];
            C[0] += (cur[0]+next[0])*(cur[0]*next[1]-next[0]*cur[1]);
            C[1] += (cur[1]+next[1])*(cur[0]*next[1]-next[0]*cur[1]);
        }
        C[0] /= 6*A;
        C[1] /= 6*A;
        new_points.push_back(C);
    }
    return tesselation(new_points,iter-1);
}

void save_svg(std::string filename, const std::vector<Polygon> &polygons = {}, const std::vector<Vector> &vectors = {}, std::string fillcol = "none");
void save_svg(std::string filename, const std::vector<Polygon> &polygons, const std::vector<Vector> &vectors, std::string fillcol){
    FILE *f = fopen(filename.c_str(), "w+");
    fprintf(f, "<svg xmlns = \"http://www.w3.org/2000/svg\" width = \"1000\" height = \"1000\">\n");
    if (!polygons.empty()){
        for (size_t i = 0; i < polygons.size(); i++){
            fprintf(f, "<g>\n");
            fprintf(f, "<polygon points = \"");
            for (size_t j = 0; j < polygons[i].vertices.size(); j++) fprintf(f, "%3.3f, %3.3f ", (polygons[i].vertices[j][0] * 1000), (1000 - polygons[i].vertices[j][1] * 1000));
            fprintf(f, "\"\nfill = \"%s\" stroke = \"black\"/>\n", fillcol.c_str());
            fprintf(f, "</g>\n");
        }
    }
    if (!vectors.empty()){
        for (auto &vect : vectors){
            fprintf(f, "<g>\n");
            fprintf(f, "<circle cx = \"%3.3f\" cy = \"%3.3f\" r=\"5\"/>", vect[0]*1000, 1000-vect[1]*1000);
            fprintf(f, "</g>\n");
        }
    }
    fprintf(f, "</svg>\n");
    fclose(f);
}
void save_svg(std::string filename, const std::vector<Vector> &vectors = {}, const std::vector<Polygon> &polygons = {}, std::string fillcol = "none");
void save_svg(std::string filename, const std::vector<Vector> &vectors, const std::vector<Polygon> &polygons, std::string fillcol){
    FILE *f = fopen(filename.c_str(), "w+");
    fprintf(f, "<svg xmlns = \"http://www.w3.org/2000/svg\" width = \"1000\" height = \"1000\">\n");
    for (size_t i = 0; i < polygons.size(); i++){
        fprintf(f, "<g>\n");
        fprintf(f, "<polygon points = \"");
        for (size_t j = 0; j < polygons[i].vertices.size(); j++) fprintf(f, "%3.3f, %3.3f ", (polygons[i].vertices[j][0] * 1000), (1000 - polygons[i].vertices[j][1] * 1000));
        fprintf(f, "\"\nfill = \"%s\" stroke = \"black\"/>\n", fillcol.c_str());
        fprintf(f, "</g>\n");
    }
    for (auto &vect : vectors){
        fprintf(f, "<g>\n");
        fprintf(f, "<circle cx = \"%3.3f\" cy = \"%3.3f\" r=\"5\"/>", vect[0] * 1000, 1000 - vect[1] * 1000);
        fprintf(f, "</g>\n");
    }
    fprintf(f, "</svg>\n");
    fclose(f);
}

int main(int argc, char **argv){
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    
    // For Polygon clipping (lab 6)

    // Polygon subjectPolygon = Polygon({Vector(0.4, 0.7, 0), Vector(0.15, 0.8, 0), Vector(0.2, 0.65, 0), Vector(0.1, 0.5, 0), Vector(0.45, 0.3, 0), Vector(0.5, 0.2, 0), Vector(0.55, 0.3, 0), Vector(0.9, 0.5, 0), Vector(0.8, 0.65, 0), Vector(0.85, 0.8, 0), Vector(0.6, 0.7, 0),});
    // Polygon clipPolygon = Polygon({Vector(0.6, 0.55, 0), Vector(0.6, 0.85, 0), Vector(0.9, 0.85, 0), Vector(0.9, 0.55, 0)});
    // Polygon clippedPolygon = clip_poly(subjectPolygon,clipPolygon);
    // save_svg("image1.svg", {subjectPolygon, clipPolygon});
    // save_svg("image2.svg", {clippedPolygon});
    
    // For the Voronoi diagrams (lab 6)
    std::vector<Vector> p3 = random_points(100);
    // std::vector<Polygon> vor = create_voronoi(p3);
    // save_svg("image.svg",p3,vor);

    // For the Centroidal Voronoi Tessellation (lab 6)
    // auto new_points = tesselation(p3);
    // save_svg("image.svg",new_points,create_voronoi(new_points));

    // For Power diagrams (lab 7)
    std::vector<double> weight = {};
    for (int i = 0; i < p3.size();i++) {
        if (p3[i][0] < 0.2 || p3[i][0] < 0.7 || p3[i][1] > 0.7) weight.push_back(0.2);
        else weight.push_back(0.8);
    } 
    save_svg("image.svg", p3, create_voronoi(p3,weight));

    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cout << "Success! Saved in 'image.svg'" << std::endl;
    std::cout << "Total time : " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << " ms" << std::endl;
    return 0;
}
