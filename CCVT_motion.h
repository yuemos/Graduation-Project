#pragma once
#ifndef CCVT_MOTION_H
#define CCVT_MOTION_H

#include <cstdio>
#include <iostream>
#include <vector>
#include <math.h>

#define PI 3.141592653590

namespace ccvt {

	struct MotionCircle {

		MotionCircle() :x(0), y(0), w(0), t(0), radiu(0), angle_Original(0) {

		}

		MotionCircle(const double dx, const  double dy, const  double dradiu, const double dw, const double dangle_Original) {
			x = dx;
			y = dy;
			radiu = dradiu;
			w = dw;
			angle_Original = dangle_Original;
		}

		//以下两个函数原设计为获得点的当前位置，后来发现无用
		inline double getX() {
			return radiu * cos((angle_Original + w * t)*PI / 180);
		}

		inline double getY() {
			return radiu * sin((angle_Original + w * t)*PI / 180);
		}

		//初始化设置点的位置
		void setMotionCircle(double dx, double dy, double dradiu, double dw, double dangle_Original) {
			x = dx;
			y = dy;
			radiu = dradiu;
			w = dw;
			angle_Original = dangle_Original;
		}

		//更新点的位置
		inline void update() {
			this->t++;
		}

		//预测下一步的点的位置
		inline void forecast(double& point_x, double& point_y) {
			point_x = cos((angle_Original + w * (t + 1))*PI / 180);
			point_y = sin((angle_Original + w * (t + 1))*PI / 180);
		}

		double x, y;// 
		double radiu;
		double w;
		double angle_Original;
		double t = 0;
	};

}

#endif // !CCVT_MOTION_H



