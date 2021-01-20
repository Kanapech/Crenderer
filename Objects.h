struct Point2D { 
    float x, y;
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
};

struct Triangle
{
    Point2D points[3];
};
