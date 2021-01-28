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

Model::Model(string filename) : verts_(), faces_(), norms_(), uv_() {
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
        }
    }
    load_texture(filename, "_diffuse.tga",texture_);
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

void Model::load_texture(string filename, const char *suffix, TGAImage &img) {
    string texfile(filename);
    size_t dot = texfile.find_last_of(".");
    if (dot!=string::npos) {
        texfile = texfile.substr(0,dot) + std::string(suffix);
        std::cerr << "texture_ file " << texfile << " loading " << (img.read_tga_file(texfile.c_str()) ? "ok" : "failed") << std::endl;
        img.flip_vertically();
    }
}

TGAColor Model::diffuse(Vec2i uv) {
    return texture_.get(uv.x, uv.y);
}

Vec2i Model::uv(int iface, int nvert) {
    int idx = faces_[iface][nvert][1];
    return Vec2i(uv_[idx].x*texture_.get_width(), uv_[idx].y*texture_.get_height());
}