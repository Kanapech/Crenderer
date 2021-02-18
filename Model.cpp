#include <cmath>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <vector>
#include "tgaimage.h"
#include "Objects.h"
#include "Model.h"
#include <sstream>
using namespace std;

Model::Model(string filename) : verts_(), faces_(), norms_(), uv_(), normalmap_(), specularmap_() {
    string v, ligne;
    ifstream file(filename);
	if (!file){
        cout << "Fichier non trouvé" << endl;
        exit(EXIT_FAILURE);
    }
    while (getline(file, ligne)) {   
        std::istringstream iss(ligne.c_str());
        char trash;
        if(ligne.rfind("v ", 0) == 0){
            Vec3f p;
            iss >> v >> p.x >> p.y >> p.z;
            verts_.push_back(p);         
        }
        if(ligne.rfind("f ", 0) == 0){
            vector<Vec3i> f;
            Vec3i tmp;
            iss >> trash;
            while (iss >> tmp[0] >> trash >> tmp[1] >> trash >> tmp[2]) {
                for (int i=0; i<3; i++) tmp[i]--;
                f.push_back(tmp);
            }
            faces_.push_back(f);
        }
        if(ligne.rfind("vt ", 0) == 0){
            Vec2f p;
            iss >> v >> p.x >> p.y >> trash;
            uv_.push_back(p);
        }
        if(ligne.rfind("vn ", 0) == 0){
            Vec3f p;
            iss >> v >> p.x >> p.y >> p.z;
            norms_.push_back(p);
        }
    }
    load_texture(filename, "_diffuse.tga", texture_);
    load_texture(filename, "_nm.tga", normalmap_);
    load_texture(filename, "_spec.tga", specularmap_);
}

Model::~Model() {
}

int Model::nbverts() {
    return (int)verts_.size();
}

int Model::nbfaces() {
    return (int)faces_.size();
}

vector<int> Model::face(int idx) {
    std::vector<int> face;
    for (int i=0; i<(int)faces_[idx].size(); i++) face.push_back(faces_[idx][i][0]);
    return face;
}

Vec3f Model::vert(int i) {
    return verts_[i];
}

Vec3f Model::vert(int iface, int nthvert) {
    return verts_[faces_[iface][nthvert][0]];
}

Vec3f Model::norm(int iface, int nvert) {
    int idx = faces_[iface][nvert][2];
    return norms_[idx].normalize();
}

Vec3f Model::normal(Vec2f uvf) {
    Vec2i uv(uvf[0]*normalmap_.get_width(), uvf[1]*normalmap_.get_height());
    TGAColor c = normalmap_.get(uv[0], uv[1]);
    Vec3f res;
    for (int i=0; i<3; i++)
        res[2-i] = (float)c[i]/255.f*2.f - 1.f;
    return res;
}

void Model::load_texture(string filename, const char *suffix, TGAImage &img) {
    string texfile(filename);
    size_t dot = texfile.find_last_of(".");
    if (dot!=string::npos) {
        texfile = texfile.substr(0,dot) + std::string(suffix);
        img.read_tga_file(texfile.c_str());
        img.flip_vertically();
    }
}

TGAColor Model::diffuse(Vec2i uv) {
    return texture_.get(uv.x, uv.y);
}

float Model::specular(Vec2f uvf) {
    Vec2i uv(uvf[0]*specularmap_.get_width(), uvf[1]*specularmap_.get_height());
    return specularmap_.get(uv[0], uv[1])[0]/1.f;
}

Vec2i Model::uv(int iface, int nvert) {
    int idx = faces_[iface][nvert][1];
    return Vec2i(uv_[idx].x*texture_.get_width(), uv_[idx].y*texture_.get_height());
}