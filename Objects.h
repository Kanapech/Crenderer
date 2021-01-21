struct Point2D { 
    int x, y;
    inline Point2D operator-(Point2D p) {
        Point2D r;
        r.x = x-p.x;
        r.y = y-p.y;
        return r;
    }
    inline Point2D operator+(Point2D p) {
        Point2D r;
        r.x = x+p.x;
        r.y = y+p.y;
        return r;
    }
    inline Point2D operator *(float f){
        Point2D r;
        r.x = x*f;
        r.y = y*f;
        return r;
    }          
};

struct Point3D { 
    float x, y, z;
    inline Point3D operator-(Point3D p) {
        Point3D r;
        r.x = x-p.x;
        r.y = y-p.y;
        r.z = z-p.z;
        return r;
    }
    inline Point3D operator ^(Point3D p) {
        Point3D r;
        r.x = y*p.z-z*p.y;
        r.y = z*p.x-x*p.z;
        r.z = x*p.y-y*p.x;
        return r; 
    }
    inline Point3D operator *(float f){
        Point3D r;
        r.x = x*f;
        r.y = y*f;
        r.z = z*f;
        return r;
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
