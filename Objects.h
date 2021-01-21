struct Point2D { 
    float x, y;
    inline Point2D operator-(Point2D p) {
        return Point2D{x-p.x, y-p.y};
    }
    inline Point2D operator+(Point2D p) {
        return Point2D{x+p.x, y+p.y};
    }
    inline Point2D operator *(float f){
        return Point2D{x*f, y*f};
    }          
};

struct Point3D { 
    float x, y, z;
    inline Point3D operator-(Point3D p) {
        return Point3D{x-p.x, y-p.y, z-p.z};
    }
    inline Point3D operator ^(Point3D p) {
        return Point3D{y*p.z-z*p.y, z*p.x-x*p.z, x*p.y-y*p.x}; 
    }
    inline Point3D operator *(float f){
        return Point3D{x*f, y*f, z*f};
    }  
    inline float operator *(Point3D p){
        return x*p.x + y*p.y+ z*p.z;
    }  
    float norm () const { return std::sqrt(x*x+y*y+z*z); }
	Point3D & normalize(int l=1) { *this = (*this)*(l/norm()); return *this; } 
};

struct Triangle
{
    Point2D points[3];
};
