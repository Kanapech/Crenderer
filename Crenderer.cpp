#include <cmath>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <vector>
#include "tgaimage.h"
#include "Objects.h"
#include <sstream>
using namespace std;

const TGAColor white = TGAColor(255, 255, 255, 255);
const TGAColor red   = TGAColor(255, 0,   0,   255);
const TGAColor blue  = TGAColor(0, 0,   255,   255);
const TGAColor green  = TGAColor(0, 255,   0,   255);

const int width = 1000;
const int height = 1000;

void line(int x0, int y0, int x1, int y1, TGAImage &image, TGAColor color) {
    bool steep = false;
    if (std::abs(x0-x1)<std::abs(y0-y1)) {
        std::swap(x0, y0);
        std::swap(x1, y1);
        steep = true;
    }
    if (x0>x1) {
        std::swap(x0, x1);
        std::swap(y0, y1);
    }
    int dx = x1-x0;
    int dy = y1-y0;
    float derror = std::abs(dy/float(dx));
    float error = 0;
    int y = y0;
    for (int x=x0; x<=x1; x++) {
        if (steep) {
            image.set(y, x, color);
        } else {
            image.set(x, y, color);
        }
        error += derror;

        if (error>.5) {
            y += (y1>y0?1:-1);
            error -= 1.;
        }
    }
}

Point3D cross(Point3D v1,Point3D v2) {
    return Point3D{v1.y*v2.z - v1.z*v2.y, v1.z*v2.x - v1.x*v2.z, v1.x*v2.y - v1.y*v2.x};
}

Point3D barycentric(Point3D A, Point3D B, Point3D C, Point3D P) {
    Point3D s[2];

    s[0].x = C.x-A.x;
    s[0].y = B.x-A.x;
    s[0].z = A.x-P.x;

    s[1].x = C.y-A.y;
    s[1].y = B.y-A.y;
    s[1].z = A.y-P.y;

    Point3D u = cross(s[0], s[1]);
    if (std::abs(u.z)>1e-2)
        return Point3D{1.f-(u.x+u.y)/u.z, u.y/u.z, u.x/u.z};
    return Point3D{-1,1,1};
}

void triangle(Point3D *pts, float *zbuffer, TGAImage &image, TGAColor color) {
    Point2D bboxmin{std::numeric_limits<float>::max(), std::numeric_limits<float>::max()};
    Point2D bboxmax{-std::numeric_limits<float>::max(), -std::numeric_limits<float>::max()};
    Point2D clamp{(float)image.get_width()-1, (float) image.get_height()-1};

    for (int i=0; i<3; i++) {
        bboxmin.x = std::max(0.f, std::min(bboxmin.x, pts[i].x));
        bboxmax.x = std::min(clamp.x, std::max(bboxmax.x, pts[i].x));

        bboxmin.y = std::max(0.f, std::min(bboxmin.y, pts[i].y));
        bboxmax.y = std::min(clamp.y, std::max(bboxmax.y, pts[i].y));
    }

    Point3D P;
    for (P.x=bboxmin.x; P.x<=bboxmax.x; P.x++) {
        for (P.y=bboxmin.y; P.y<=bboxmax.y; P.y++) {
            Point3D bc_screen  = barycentric(pts[0], pts[1], pts[2], P);
            if (bc_screen.x<0 || bc_screen.y<0 || bc_screen.z<0) continue;
            P.z = 0;
            P.z += pts[0].z*bc_screen.x;
            P.z += pts[1].z*bc_screen.y;
            P.z += pts[2].z*bc_screen.z;

            if (zbuffer[int(P.x+P.y*width)]<P.z) {
                zbuffer[int(P.x+P.y*width)] = P.z;
                image.set(P.x, P.y, color);
            }
        }
    }
}

Point3D world2screen(Point3D v) {
    return Point3D{(float) int((v.x+1.)*width/2.+.5), (float) int((v.y+1.)*height/2.+.5), v.z};
}

int main(int argc, char** argv) {
    string v, ligne;
    float vert[3][1258];
    int faces[3][2492];
    float x, y, z; 
    int vcount = 0;
    int fcount = 0;
    string a, b, c;
    TGAImage image(width, height, TGAImage::RGB);
    
    ifstream file("african_head.obj");
	if (!file){
        cout << "Fichier non trouvé" << endl;
        exit(EXIT_FAILURE);
    }
    while (getline(file, ligne)) {   
        std::istringstream iss(ligne.c_str());
        if(ligne.rfind("v ", 0) == 0){
            iss >> v >> x >> y >> z;
            vert[0][vcount] = x; vert[1][vcount] = y; vert[2][vcount] = z;
            //cout << vert[0][vcount] << vert[1][vcount] << vert[2][vcount] << endl;
            //line(vert[0][vcount], vert[1][vcount], vert[0][vcount], vert[1][vcount], image, red);
            vcount++;            
        }
        if(ligne.rfind("f ", 0) == 0){
            iss >> v >> a >> b >> c;
            a = a.substr(0, a.find("/"));
            b = b.substr(0, b.find("/"));
            c = c.substr(0, c.find("/"));
            faces[0][fcount] = stoi(a)-1 ; faces[1][fcount] = stoi(b)-1; faces[2][fcount] = stoi(c)-1;
            //cout << faces[0][fcount] << faces[1][fcount] << faces[2][fcount] << endl;
            fcount++;
        }
    }

    Point3D light_dir{0,0,-1};
    float *zbuffer = new float[width*height];
    for (int i=width*height; i--; zbuffer[i] = -std::numeric_limits<float>::max());

    for (int i=0; i<fcount; i++) {
        std::vector<int> face{faces[0][i], faces[1][i], faces[2][i]};
        Point3D pts[3];
        Point3D world_coords[3]; 
        for (int i=0; i<3; i++){
            Point3D v{vert[0][face[i]], vert[1][face[i]], vert[2][face[i]]};
            world_coords[i] = v;
            pts[i] = world2screen(v);
        } 
        Point3D n = (world_coords[2]-world_coords[0])^(world_coords[1]-world_coords[0]);
        n.normalize();
        float intensity = n*light_dir;
        if (intensity>0){
            triangle(pts, zbuffer, image, TGAColor(intensity*255, intensity*255, intensity*255, 255));
        }
    }

    image.flip_vertically(); // i want to have the origin at the left bottom corner of the image
    image.write_tga_file("output.tga");

    file.close(); 
    return EXIT_SUCCESS;
}