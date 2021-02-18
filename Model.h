#ifndef __MODEL_H__
#define __MODEL_H__

#include <vector>
#include <iostream>
#include "tgaimage.h"
#include "geometry.h"
using namespace std;

class Model {
    private:
        vector<Vec3f> verts_;
        vector<vector<Vec3i>> faces_;
        vector<Vec3f> norms_;
        vector<Vec2f> uv_;
        TGAImage texture_;
        TGAImage normalmap_;
        TGAImage specularmap_;
        void load_texture(string filename, const char *suffix, TGAImage &img);
    public:
        Model(string filename);
        ~Model();
        int nbverts();
        int nbfaces();
        Vec3f vert(int i);
        Vec3f vert(int iface, int nthvert);
        Vec2i uv(int iface, int nvert);
        TGAColor diffuse(Vec2i uv);
        float specular(Vec2f uvf);
        vector<int> face(int idx);
        Vec3f norm(int i, int j);
        Vec3f normal(Vec2f uvf);
};

#endif //__MODEL_H__