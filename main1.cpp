#include <iostream>
#include <fstream>
#include <cmath>
#include <float.h>
#include <list>

using namespace std;

class Vec {
public:

    double X;
    double Y;
    double Z;

    Vec(double x, double y, double z) : X(x), Y(y), Z(z) {}

    Vec operator +(const Vec& vec) const{
        return Vec(X+vec.X,Y+vec.Y,Z+vec.Z);
    }

    Vec operator -(const Vec& vec) const {
        return Vec(X-vec.X,Y-vec.Y,Z-vec.Z);
    }

    Vec operator *(double scalar) const{
        return Vec(scalar*X,scalar*Y,scalar*Z);
    }


    double dot(const Vec& vec) const {
        return X*vec.X+Y*vec.Y+Z*vec.Z;
    }

    double norm() const {
        return sqrt(dot(*this));
    }
	double norm2() const {
		return dot(*this);
	}
    
    Vec normalize() {
        return (*this)*(1/norm());
    }
    
};

class Ray {
public:
    Vec origin;
    Vec direction;

	Ray(Vec ori, Vec dir) : origin(ori), direction(dir) {}
	
	Vec get_point(double t) const{
		return origin + direction * t;
	}
	Vec reflect_by(const Vec& normal) const{
		return direction - normal * normal.dot(direction) * 2;
	}

};

class Color  {

public:
	int r, g, b;
	Color(int red, int green, int blue) : r(red), g(green), b(blue) {}
	
	Color scale_by(double scalar)const {
		return (scalar > 0) ? Color_trunc(scalar * r, scalar * g, scalar * b) : Color(0,0,0);
	}
	Color scale_by2(double scalar)const {
		if (scalar < 0)
			return  Color(0, 0, 0);
		scalar *= scalar;
		return Color_trunc(scalar * r, scalar * g, scalar * b);
	}
	Color mix_with(const Color& rhs) const {
		return Color(r * rhs.r, g * rhs.g, b * rhs.b);
	}
	Color operator *(double scalar)const {
		return Color(scalar * r, scalar * g, scalar * b);
	}
	Color operator +(const Color& rhs) const {
		return Color_trunc(r + rhs.r, g + rhs.g, b + rhs.b);
	}
	static int trunc(int val) {
		return (val > 255) ? 255 : val;
	}
	static Color Color_trunc(int red, int green, int blue) {
		return Color(trunc(red), trunc(green), trunc(blue));
	}
};


class Sphere  {
    Vec Center;
    double Radius;
    
public:
    Color color;
    
    Sphere(Vec center, double radius, Color color_) : Center(center), Radius(radius),color(color_){}
    
    Vec get_center() const{
        return Center;
    }

    Vec get_normal(const Vec& p) const {
		return ((p - Center)*(-1/Radius)).normalize();
    }

    bool intersect(const Ray& ray, double &t) const {
        Vec v = ray.origin - Center;

        const double b = 2 * v.dot(ray.direction);
        const double c = v.dot(v) - Radius*Radius;
        double delta = b*b - 4 * c;

		if (delta < 0) {
			t = FLT_MAX; 
			return false;
		}

        const double t1 = (-b - sqrt(delta))/2;
        const double t2 = (-b + sqrt(delta))/2;

		if (t2 < 1e-2) { 
			t = FLT_MAX; 
			return false;
		}

		t = (t1 >= 1e-2) ? t1 : t2; 
        
        return true;
    }
     Color Shading(const Ray& ray, const Sphere& object, double t, int depth) {          
            Vec intersect_point = ray.origin + ray.direction * t;
            Vec normal = object.get_normal(intersect_point);
    
				return (object.color).scale_by(normal.dot(ray.direction) * 0.5);
			
			
			}
	
};


class Lightsource {
public:
	Vec position;
	Color color;
	double intensity = 100;
	Lightsource(Vec position_, Color color_ = Color(255,255,255), double intensity_=100.0) : position(position_), color(color_), intensity(intensity_) {}
};

class Scene {
   
    list<Sphere*> objects;
	list<Lightsource> lightsources;

public:
    Scene() {}

    void add(Sphere *object) {
		objects.push_back(object);
    }
	
	void add(Lightsource light) {
		lightsources.push_back(light);
	}
    
    Color Shading(const Ray& ray, const Sphere& object, double t, int depth) {          
            Vec intersect_point = ray.origin + ray.direction * t;
            Vec normal = object.get_normal(intersect_point);
    
				return (object.color).scale_by(normal.dot(ray.direction) * 0.5);
			
			
			}
	
    Color trace(int x, int y) {
		Vec ray_origin = Vec(0, 0, -1000);
		Vec ray_direction = Vec(x, y, 1250).normalize();

		return trace_ray(Ray(ray_origin, ray_direction), 0, 50);
	}
	Color trace_ray(const Ray& ray, const Sphere* exclude_obj, int depth) {
		double min_t = FLT_MAX;
		const Sphere* nearest_obj = nullptr;

		double t = FLT_MAX;
		for (const Sphere* object : objects) {
				if ((*object).intersect(ray, t)) {
					if (min_t > t) {
						nearest_obj = object;
			
						min_t = t;
					}
				}
		}
		if (nearest_obj != nullptr) {
			return Shading(ray, *nearest_obj, min_t, depth);
		}
		return Color(128, 128, 128);
	}
	

};




int main() {
    
    int Height = 700;
    int Width = 700;
    
    Color pix_col(0,0,0);
    
    Scene scene = Scene();
    
	Sphere sphere1 (Vec(100, 140, 265), 75, Color(40, 255, 100));
	Sphere sphere2 (Vec(300, 400, 265), 120, Color(255, 0, 255));

	Lightsource light = Lightsource(Vec(0, 0, 0), Color(255,255,255));
    
    scene.add(&sphere1);
    scene.add(&sphere2);

	
	scene.add(light);

    ofstream my_Image ("image.ppm");

    if (my_Image.is_open ()) {
        my_Image << "P3\n" << Width << " " << Height << " 255\n";
        for (int i = 0; i < Height; i++) {
            for (int j = 0; j < Width; j++)  {

            pix_col = scene.trace(i, j);


            my_Image << (int)pix_col.r << ' ' << (int)pix_col.g << ' ' << (int)pix_col.b << "\n";
            }
        }
        my_Image.close();
    }
    else
      cout << "Could not open the file";
    
    return 0;
}
