#ifndef CCVT_POINT_H
#define CCVT_POINT_H

#include <assert.h>
#include <list>
#include "ccvt_motion.h"
#include <vector>

namespace ccvt {

	struct Point2 {

		typedef std::list<Point2> List;
		typedef std::vector<Point2> Vector;

		Point2() :x(0), y(0) {
		}

		Point2(const double x, const double y) :x(x), y(y) {
		}

		Point2(const double x, const double y, int father) :x(x), y(y), father(father) {
		}
		
		Point2(const double dx, const double dy, const double dradiu, const double dw, const double dangle_Original) {
			this->MC = MotionCircle(dx, dy, dradiu, dw, dangle_Original);
			this->MC.update();
			this->x = this->MC.getX();
			this->y = this->MC.getY();
		}

		Point2(const MotionCircle& MC) :MC(MC) {
			this->MC.update();
			x = this->MC.getX();
			y = this->MC.getY();
			//std::cout << this->MC.radiu << std::endl;
		}
		
		Point2(const Point2 & p) {
			x = p.x;
			y = p.y;
			MC = p.MC;
			father = p.father;
		}

		inline const double operator[](const int i) const {
			assert(i == 0 || i == 1);
			if (i == 0) {
				return x;
			}
			return y;
		}

		inline double operator[](const int i) {
			assert(i == 0 || i == 1);
			if (i == 0) {
				return x;
			}
			return y;
		}

		inline int getFather() {
			return father;
		}

		inline void changeFather(int fa) {
			father = fa;
		}

		inline void update() {
			//printf("%lf %lf||", x, y);
			//x = radiu * cos((angle_Original + w * t)*PI / 180);
			//y = radiu * sin((angle_Original + w * t)*PI / 180);
			MC.update();
			x = MC.getX();
			int xx = static_cast<int>(x * 100+0.5);
			x = xx / 100.0;
			y = MC.getY();
			int yy = static_cast<int>(y * 100 + 0.5);
			y = yy / 100.0;
			//printf("%lf %lf\n", x, y);
		}

		double x;
		double y;
		int father = -1;
		MotionCircle MC;
/*		double dx, dy;
		double radiu;
		double w;
		double angle_Original;
		int t = 0;
*/
		static const int D = 2;
	};


	struct Point3 {

		typedef std::list<Point3> List;
		typedef std::vector<Point3> Vector;

		Point3() :x(0), y(0), z(0) {
		}

		Point3(const double x, const double y, const double z) :x(x), y(y), z(z) {
		}

		Point3(const Point3& p) :x(p.x), y(p.y), z(p.z) {
		}

		inline const double operator[](const int i) const {
			assert(i == 0 || i == 1 || i == 2);
			if (i == 0) {
				return x;
			}
			if (i == 1) {
				return y;
			}
			return z;
		}

		inline double operator[](const int i) {
			assert(i == 0 || i == 1 || i == 2);
			if (i == 0) {
				return x;
			}
			if (i == 1) {
				return y;
			}
			return z;
		}

		double x;
		double y;
		double z;

		static const int D = 3;
	};

	template<int d>
	struct Point {
		typedef std::list<Point> List;
		typedef std::vector<Point> Vector;

		Point() {
			for (int i = 0; i < d; i++) {
				coordinates[i] = 0.0;
			}
		}

		Point(const double coordinates[d]) {
			for (int i = 0; i < d; i++) {
				this->coordinates[i] = coordinates[i];
			}
		}

		Point(const Point& p) {
			for (int i = 0; i < d; i++) {
				this->coordinates[i] = p.coordinates[i];
			}
		}

		inline const double operator[](const int i) const {
			assert(i >= 0 && i < d);
			return coordinates[i];
		}

		inline double operator[](const int i) {
			assert(i >= 0 && i < d);
			return coordinates[i];
		}

		static const int D = d;

	private:
		double coordinates[d];
	};
}

#endif // !1