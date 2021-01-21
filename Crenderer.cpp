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

const int width = 600;
const int height = 600;

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

void triangle(Point2D t0, Point2D t1, Point2D t2, TGAImage &image, TGAColor color) {
    if (t0.y==t1.y && t0.y==t2.y) return;
    if (t0.y>t1.y) std::swap(t0, t1);
    if (t0.y>t2.y) std::swap(t0, t2);
    if (t1.y>t2.y) std::swap(t1, t2);
    int total_height = t2.y-t0.y;
    for (int i=0; i<total_height; i++) {
        bool second_half = i>t1.y-t0.y || t1.y==t0.y;
        int segment_height = second_half ? t2.y-t1.y : t1.y-t0.y;
        float alpha = (float)i/total_height;
        float beta  = (float)(i-(second_half ? t1.y-t0.y : 0))/segment_height;
        Point2D A = t0 + (t2-t0)*alpha;
        Point2D B = second_half ? t1 + (t2-t1)*beta : t0 + (t1-t0)*beta;
        if (A.x>B.x) std::swap(A, B);
        for (int j=A.x; j<=B.x; j++) {
            image.set(j, t0.y+i, color);
        }
    }
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
            //vert[0][vcount] = (x+1.)*width/2. ; vert[1][vcount] = (y+1.)*height/2.; vert[2][vcount] = (z+1.)*height/2.;
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

    /*for(int i = 0; i<fcount; i++){
        Point2D t0, t1, t2;
        TGAColor couleur(std::rand()%255,std::rand()%255,std::rand()%255,255);
        t0.x = vert[0][faces[0][i]];
        t0.y = vert[1][faces[0][i]];

        t1.x = vert[0][faces[1][i]];
        t1.y = vert[1][faces[1][i]];

        t2.x = vert[0][faces[2][i]];
        t2.y = vert[1][faces[2][i]];

        triangle(t0, t1, t2, image, couleur);
    }*/
    int tcount=0;
    Point3D light_dir; light_dir.x = 0; light_dir.y = 0; light_dir.z = -1;
    for (int i=0; i<fcount; i++) { 
        std::vector<int> face{faces[0][i], faces[1][i], faces[2][i]}; 
        Point2D screen_coords[3]; 
        Point3D world_coords[3]; 
        for (int j=0; j<3; j++) { 
            Point3D v; v.x = vert[0][face[j]]; v.y = vert[1][face[j]]; v.z = vert[2][face[j]];
            Point2D v2; v2.x = (v.x+1.)*width/2.; v2.y = (v.y+1.)*height/2.;
            screen_coords[j] = v2;
            world_coords[j]  = v; 
        } 
        Point3D n = (world_coords[2]-world_coords[0])^(world_coords[1]-world_coords[0]);
        n.normalize();
        float intensity = n*light_dir;
        if (intensity>0){
            tcount++;
            triangle(screen_coords[0], screen_coords[1], screen_coords[2], image, TGAColor(intensity*255, intensity*255, intensity*255, 255)); 
        }
    } 
    cout << tcount << endl;
    image.flip_vertically(); // i want to have the origin at the left bottom corner of the image
    image.write_tga_file("output.tga");

    file.close(); 
    return 0;
}