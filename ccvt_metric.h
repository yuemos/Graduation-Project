#ifndef CCVT_METRIC_H
#define CCVT_METRIC_H

#include <math.h>
#include "ccvt_point.h"

namespace ccvt {
	struct MetricEuclidean2 {
		Point2 point;
		inline double distance(const Point2& p1, const Point2& p2) const {
			return sqrt(distance_square(p1, p2));
		}

		inline double distance_square(const Point2& p1, const Point2& p2) const {
			return (p1.x - p2.x) * (p1.x - p2.x) + (p1.y - p2.y) * (p1.y - p2.y);
		}

		//把所有的点的横纵坐标的平均数算出来了
		inline Point2 centroid(const Point2&, const Point2::Vector& points) const {
			Point2 centroid;//
			int pointSize = static_cast<int>(points.size());
			for (int j = 0; j < pointSize; j++) {
				centroid.x += points[j].x;
				centroid.y += points[j].y;
			}
			centroid.x /= pointSize;
			centroid.y /= pointSize;
			return centroid;
		}
	};

	struct MetricEuclidean3 {
		inline double distance(const Point3& p1, const Point3& p2) const{
			return sqrt(distance_square(p1, p2));
		}

		inline double distance_square(const Point3& p1, const Point3& p2) const {
			return (p1.x - p2.x) * (p1.x - p2.x) + (p1.y - p2.y) * (p1.y - p2.y) + (p1.z - p2.z) * (p1.z - p2.z);
		}

		inline Point3 centroid(const Point3&, const Point3::Vector& points) const {
			Point3 centroid;
			int pointSize = static_cast<int>(points.size());
			for (int j = 0; j < pointSize; j++) {
				centroid.x += points[j].x;
				centroid.y += points[j].y;
				centroid.z += points[j].z;
			}
			centroid.x /= pointSize;
			centroid.y /= pointSize;
			centroid.z /= pointSize;
			return centroid;

		}
	};

	template<class Point>
	struct MetricEuclidean {
		inline double distance(const Point& p1, const Point& p2) const {
			return sqrt(distance_square(p1, p2));
		}

		inline double distance_square(const Point& p1, const Point& p2)const {
			double distanceSquare = 0;
			for (int i = 0; i < Point::D; i++) {
				distanceSquare += (p1[i] - p2[i])*(p1[i] - p2[i]);
			}
			return distanceSquare;
		}

		inline Point centroid(const Point&, const typename Point::Vector points)  const {
			Point centroid;
			int pointSize = static_cast<int>(points.size());
			for (int j = 0; j < pointSize; j++) {
				for (int i = 0; i < Point::D; i++) {
					centroid[i] += points[j][i];
				}
			}
			for (int j = 0; j < pointSize; j++) {
				centroid[j] /= pointSize;
			}
			return centroid;
		}
	};

	struct MetricToroidalEuclidean2 {
		Point2 size;//这里的size到底代表了什么？

		MetricToroidalEuclidean2() :size(1, 1) {

		}

		MetricToroidalEuclidean2(const Point2& size) :size(size) {

		}

		inline double distance(const Point2& p1, const Point2& p2) const {
			//return distance_square(p1, p2);
			return sqrt(distance_square(p1, p2));
		}

		inline double distance_square(const Point2& p1, const Point2& p2) const {
			double dx = p1.x - p2.x;
			if (fabs(dx) > size.x / 2) {
				if (p1.x < size.x / 2) {
					dx = p1.x - (p2.x - size.x);
				}
				else {
					dx = p1.x - (p2.x + size.x);
				}
			}
			double dy = p1.y - p2.y;
			if (fabs(dy) > size.y / 2) {
				if (p1.y < size.x / 2) {
					dy = p1.y - (p2.y - size.y);
				}
				else {
					dy = p1.y - (p2.y + size.y);
				}
			}
			return dx * dx + dy * dy;
		}

		inline Point2 centroid(const Point2& center, const Point2::Vector& points) const {
			Point2 centroid;
			int pointSize = static_cast<int>(points.size());
			for (int j = 0; j < pointSize; j++) {
				double p = points[j].x;
				if (fabs(center.x - p) > size.x / 2) {
					if (center.x < size.x / 2) {
						p -= size.x;
					}
					else {
						p += size.x;
					}
				}
				centroid.x += p;
				p = points[j].y;
				if (fabs(center.y - p) > size.y / 2) {
					if (center.y < size.y) {
						p -= size.y;
					}
					else {
						p += size.y;
					}
				}
				centroid.y += p;
			}
			centroid.x /= pointSize;
			centroid.y /= pointSize;
			if (centroid.x < 0) {
				centroid.x += size.x;
			}
			if (centroid.x >= size.x) {
				centroid.x -= size.x;
			}
			if (centroid.y < 0) {
				centroid.y += size.y;
			}
			if (centroid.y >= size.y) {
				centroid.y -= size.y;
			}
			return centroid;
		}
	};

	struct MetricToroidalEuclidean3 {
		
		Point3 size;

		MetricToroidalEuclidean3() :size(1, 1, 1) {

		}

		MetricToroidalEuclidean3(const Point3& size) :size(size) {

		}

		inline double distance(const Point3& p1, const Point3& p2) const {
			return sqrt(distance_square(p1, p2));
		}

		inline double distance_square(const Point3& p1, const Point3& p2) const {
			double dx = p1.x - p2.x;
			if (fabs(dx) > size.x / 2) {
				if (p1.x < size.x / 2) {
					dx = p1.x - (p2.x - size.x);
				}
				else {
					dx = p1.x - (p2.x + size.x);
				}
			}
			double dy = p1.y - p2.y;
			if (fabs(dy) > size.y / 2) {
				if (p1.y < size.y / 2) {
					dy = p1.y - (p2.y - size.y);
				}
				else {
					dy = p1.y - (p2.y + size.y);
				}
			}
			double dz = p1.z - p2.z;
			if (fabs(dx) > size.z / 2) {
				if (p1.z < size.z / 2) {
					dz = p1.z - (p2.z - size.z);
				}
				else {
					dz = p1.z - (p2.z + size.z);
				}
			}
			return dx * dx + dy * dy + dz * dz;
		}

		inline Point3 centroid(const Point3& center, const Point3::Vector points) const {
			Point3 centroid;
			int pointSize = static_cast<int>(points.size());
			for (int j = 0; j < pointSize; j++) {
				double p = points[j].x;
				if (fabs(center.x - p) > size.x / 2) {
					if (center.x < size.x / 2) {
						p -= size.x;
					}
					else {
						p += size.x;
					}
				}
				centroid.x += p;
				p = points[j].y;
				if (fabs(center.y - p) > size.y / 2) {
					if (center.y < size.y / 2) {
						p -= size.y;
					}
					else {
						p += size.y;
					}
				}
				centroid.y += p;
				p = points[j].z;
				if (fabs(center.z - p) > size.z / 2) {
					if (center.z < size.z / 2) {
						p -= size.z;
					}
					else {
						p += size.z;
					}
				}
				centroid.z += p;
			}
			centroid.x /= pointSize;
			centroid.y /= pointSize;
			centroid.z /= pointSize;
			if (centroid.x < 0) {
				centroid.x += size.x;
			}
			if (centroid.x >= size.x) {
				centroid.x -= size.x;
			}
			if (centroid.y < 0) {
				centroid.y += size.y;
			}
			if (centroid.y >= size.y) {
				centroid.y -= size.y;
			}
			if (centroid.z < 0) {
				centroid.z += size.z;
			}
			if (centroid.z >= size.z) {
				centroid.z -= size.z;
			}
			return centroid;
		}
	};

	template<class Point>
	struct MetricToroidalEuclidean {
		Point size;

		MetricToroidalEuclidean() {
			for (int i = 0; i < Point::D; i++) {
				size[i] = 1;
			}
		}


		MetricToroidalEuclidean(const Point size) :size(size) {

		}

		inline double distance(const Point& p1, const Point& p2) const {
			return distance_square(p1, p2);
		}

		inline double distance_square(const Point& p1, const Point& p2) const {
			double distanceSquare = 0;
			for (int i = 0; i < Point::D; ++i) {
				double di = p1[i] - p2[i];
				if (fabs(di) > size[i] / 2) {
					if (p1[i] < size[i] / 2) {
						di = p1[i] - (p2[i] - size[i]);
					}
					else {
						di = p1[i] - (p2[i] + size[i]);
					}
				}
				distanceSquare += di * di;
			}
			return distanceSquare;
		}

		inline Point centroid(const Point& center, const typename Point::Vector& points) const {
			Point centroid;
			int pointsSize = static_cast<int>(points.size());
			for (int j = 0; j < pointsSize; ++j) {
				for (int i = 0; i < Point::D; ++i) {
					double pi = points[j][i];
					if (fabs(center[i] - pi) > size[i] / 2) {
						if (center[i] < size[i] / 2) {
							pi -= size[i];
						}
						else {
							pi += size[i];
						}
					}
					centroid[i] += pi;
				}
			}
			for (int i = 0; i < Point::D; ++i) {
				centroid[i] /= pointsSize;
				if (centroid[i] < 0) {
					centroid[i] += size[i];
				}
				if (centroid[i] >= size[i]) {
					centroid[i] -= size[i];
				}
			}
			return centroid;
		}
	};
}

#endif // !1

