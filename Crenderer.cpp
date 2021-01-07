#include <cmath>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <vector>
#include "tgaimage.h"
using namespace std;

const TGAColor white = TGAColor(255, 255, 255, 255);
const TGAColor red   = TGAColor(0 , 0,   255,   255);

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
            image.set(y, x, TGAColor(255, 1));
        } else {
            image.set(x, y, TGAColor(255, 1));
        }
        error += derror;

        if (error>.5) {
            y += (y1>y0?1:-1);
            error -= 1.;
        }
    }
}

int main(int argc, char** argv) {
    string v, ligne;
    //float vert[2][1258];
    float x, y, z; 
    //int count = 0;
    TGAImage image(200, 200, TGAImage::RGB);
    
    ifstream file("african_head.obj");
	if (!file) 
		cout << "Fichier non trouvé";
	
    while (getline(file, ligne)) {
        if(ligne.rfind("v ", 0) == 0){
            file >> v >>  x >> y >> z;
            //vert[0][count] = x ; vert[1][count] = z;
            //cout << vert[0][count] << vert[1][count] << endl;
            //count++;
            
            line(x*100+100, y*100+100, x*100+100, y*100+100, image, white);     
        }
    }

    image.flip_vertically(); // i want to have the origin at the left bottom corner of the image
    image.write_tga_file("output.tga");

    file.close(); 
    return 0;
}