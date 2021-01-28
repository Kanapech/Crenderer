#include <cmath>
#include <iostream>
#include <limits>
#include <vector>
#include "tgaimage.h"
#include "Model.h"
#include "geometry.h"
using namespace std;

const TGAColor white = TGAColor(255, 255, 255, 255);
const TGAColor red   = TGAColor(255, 0,   0,   255);
const TGAColor blue  = TGAColor(0, 0,   255,   255);
const TGAColor green  = TGAColor(0, 255,   0,   255);

const int width = 1000;
const int height = 1000;

Model *model;
TGAImage image(width, height, TGAImage::RGB);

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

Vec3f cross(Vec3f v1,Vec3f v2) {
    return Vec3f{v1.y*v2.z - v1.z*v2.y, v1.z*v2.x - v1.x*v2.z, v1.x*v2.y - v1.y*v2.x};
}

Vec3f barycentric(Vec3f A, Vec3f B, Vec3f C, Vec3f P) {
    Vec3f s[2];
    for (int i=2; i--; ) {
        s[i][0] = C[i]-A[i];
        s[i][1] = B[i]-A[i];
        s[i][2] = A[i]-P[i];
    }
    Vec3f u = cross(s[0], s[1]);
    if (std::abs(u[2])>1e-2)
        return Vec3f(1.f-(u.x+u.y)/u.z, u.y/u.z, u.x/u.z);
    return Vec3f(-1,1,1);
}

void triangle(Model* model, Vec3f *pts, float *zbuffer, TGAImage &image, float intensity, Vec2i* uv) {
    Vec2f bboxmin(std::numeric_limits<float>::max(), std::numeric_limits<float>::max());
    Vec2f bboxmax(-std::numeric_limits<float>::max(), -std::numeric_limits<float>::max());
    Vec2f clamp(image.get_width()-1, image.get_height()-1);

    for (int i=0; i<3; i++) {
        bboxmin.x = std::max(0.f, std::min(bboxmin.x, pts[i].x));
        bboxmax.x = std::min(clamp.x, std::max(bboxmax.x, pts[i].x));

        bboxmin.y = std::max(0.f, std::min(bboxmin.y, pts[i].y));
        bboxmax.y = std::min(clamp.y, std::max(bboxmax.y, pts[i].y));
    }

    Vec3f P;
    for (P.x=bboxmin.x; P.x<=bboxmax.x; P.x++) {
        for (P.y=bboxmin.y; P.y<=bboxmax.y; P.y++) {
            Vec3f bary  = barycentric(pts[0], pts[1], pts[2], P);
            if (bary.x<0 || bary.y<0 || bary.z<0) continue;
            P.z = 0;
            P.z += pts[0].z*bary.x;
            P.z += pts[1].z*bary.y;
            P.z += pts[2].z*bary.z;

            float u=0;
            u += uv[0].x*bary.x;
            u += uv[1].x*bary.y;
            u += uv[2].x*bary.z;

            float v=0;
            v += uv[0].y*bary.x;
            v += uv[1].y*bary.y;
            v += uv[2].y*bary.z;
            Vec2i uvP(u, v);
            

            if (zbuffer[int(P.x+P.y*width)]<P.z) {
                zbuffer[int(P.x+P.y*width)] = P.z;
                TGAColor color = model->diffuse(uvP);
                image.set(P.x, P.y, TGAColor(color.r*intensity, color.g*intensity, color.b*intensity));
            }
        }
    }
}

Vec3f world2screen(Vec3f v) {
    return Vec3f(int((v.x+1.)*width/2.+.5), (float) int((v.y+1.)*height/2.+.5), v.z);
}

int main(int argc, char** argv) {

    model = new Model("african_head.obj");
    Vec3f light_dir{0,0,-1};
    float *zbuffer = new float[width*height];
    for (int i=width*height; i--; zbuffer[i] = -std::numeric_limits<float>::max());

    for (int i=0; i<model->nbfaces(); i++) {
        vector<int> face = model->face(i);

        Vec3f pts[3];
        Vec3f world_coords[3]; 
        for (int j=0; j<3; j++){
            Vec3f v = model->vert(face[j]);
            world_coords[j] = v;
            pts[j] = world2screen(v);
        }

        Vec3f n = (world_coords[2]-world_coords[0])^(world_coords[1]-world_coords[0]);
        n.normalize();
        float intensity = n*light_dir;
        if (intensity>0){
            Vec2i uv[3];
            for(int k=0; k<3; k++){
                uv[k] = model->uv(i, k);
            } 
            triangle(model, pts, zbuffer, image, intensity, uv);
            //triangle(pts, zbuffer, image, TGAColor(255*intensity, 255*intensity, 255*intensity));
        }
    }

    image.flip_vertically(); // i want to have the origin at the left bottom corner of the image
    image.write_tga_file("output.tga");
    return EXIT_SUCCESS;
}