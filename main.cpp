#include <fstream>
#include <sstream>
#include <cstring>
#include <time.h>
#include <iostream>
#include "ccvt_metric.h"
#include "ccvt_optimizer.h"
#include "ccvt_point.h"
#include "ccvt_site.h"
#include "CCVT_motion.h"

using namespace ccvt;

std::ofstream oFile("D:/workspace/CCVT/out2.txt");

// discrete space with constant density;
// the points form a regular grid
void constant_regular_density(Point2::List& points, const int numberOfPoints, const double torusSize) {
	double n = sqrt(static_cast<double>(numberOfPoints));
	for (int x = 0; x < n; ++x) {
		for (int y = 0; y < n; ++y) {
			double dx = x / n * torusSize;
			double dy = y / n * torusSize;
			points.push_back(Point2(dx, dy));
		}
	}
}

// discrete space with constant density;
// the points are randomly distributed
void constant_random_density(Point2::List& points, const int numberOfPoints, const double torusSize) {
	for (int i = 0; i < numberOfPoints; ++i) {
		double x = static_cast<double>(rand() % RAND_MAX) / RAND_MAX * torusSize;
		double y = static_cast<double>(rand() % RAND_MAX) / RAND_MAX * torusSize;
		points.push_back(Point2(x, y));
	}
}

// discrete space with the density function e^(-20x^2-20y^2)+0.2sin^2(PIx)sin^2(PIy);
// the points are generated via rejection sampling
void nonconstant_density(Point2::List& points, const int numberOfPoints, const double torusSize) {
	const double E = 2.718281828459;
	while (points.size() < static_cast<unsigned int>(numberOfPoints)) {
		double x = static_cast<double>(rand() % RAND_MAX) / RAND_MAX * 2 - 1;
		double y = static_cast<double>(rand() % RAND_MAX) / RAND_MAX * 2 - 1;
		double p = pow(E, -20.0 * x * x - 20.0 * y * y) + 0.2 * sin(PI * x) * sin(PI * x) * sin(PI * y) * sin(PI * y);
		double r = static_cast<double>(rand() % RAND_MAX) / RAND_MAX;
		if (p >= r) {
			points.push_back(Point2((x + 1) / 2 * torusSize, (y + 1) / 2 * torusSize));
		}
	}
}

//���ļ�����̶�����,�������ļ���һ����site,һ����points
void readFile(int &NUMBER_SITES,int &NUMBER_POINTS,Point2::List& points, Site<Point2>::List& sites) {
	std::ifstream iFile("D:/workspace/CCVT/in.txt");
	std::string s;
	std::stringstream ss;

	if (iFile) {
		getline(iFile, s);
		int n = s.find(" ");
		std::string a = s.substr(0, n);
		std::string b = s.substr(n + 1, s.length() - n - 1);
		ss << a, ss >> NUMBER_SITES;
		ss.clear();
		ss << b, ss >> NUMBER_POINTS;
		ss.clear();

		for (int i = 0; i < NUMBER_SITES; i++) {
			double x, y;
			getline(iFile, s);
			int n = s.find(" ");
			std::string a = s.substr(0, n);
			std::string b = s.substr(n + 1, s.length() - n - 1);
			ss << a, ss >> x;
			ss.clear();
			ss << b, ss >> y;
			ss.clear();
			
			sites.push_back(Site<Point2>(i, NUMBER_POINTS, Point2(x, y)));
		}
		for (int i = 0; i < NUMBER_SITES*NUMBER_POINTS; i++) {
			double dx, dy;
			getline(iFile, s);
			int n = s.find(" ");
			std::string a = s.substr(0, n);
			std::string b = s.substr(n + 1, s.length() - n - 1);
			ss << a, ss >> dx;
			ss.clear();
			ss << b, ss >> dy;
			ss.clear();
			double dradiu,dangle;
			getline(iFile, s);
			n = s.find(" ");
			a = s.substr(0, n);
			b = s.substr(n + 1, s.length() - n - 1);
			ss << a, ss >> dradiu;
			ss.clear();
			ss << b, ss >> dangle;
			ss.clear();
			double dw;
			getline(iFile, s);
			ss << s, ss >> dw;
			ss.clear();
			//std::cout << radiu << " " << angle << std::endl;
			double d = 0.0;
			double w = 1.0;
			points.push_back(Point2(dx, dy, dradiu, dw, double(dangle)));
		}
	}
//	printf("(%f,%f)\n", sites.front().location.x, sites.front().location.y);
//	printf("(%f,%f)\n", sites.front().location.x, sites.front().location.y);
//	printf("阿斯顿发斯蒂芬阿萨德安抚");
//	for (int i = 0; i < points.size(); i++) {
//		Point2 p = points.front();
//		printf("(%f,%f)\n", p.x, p.y);
//		points.pop_front();
//	}
	if (iFile.is_open()) {
		iFile.close();
	}
}

// export sites to an EPS image
bool save_eps(const char* filename, const Site<Point2>::Vector& sites, const double width, const double height, const double radius, Optimizer<Site<Point2>, Point2, MetricEuclidean2> optimizer) {
	std::ofstream stream(filename, std::ios::out);
	if (stream.bad()) {
		return false;
	}

	stream << "%!PS-Adobe EPSF-3.0\n";
	stream << "%%HiResBoundingBox: " << -width << " " << -height << " " << width << " " << height << "\n";
	stream << "%%BoundingBox: " << -static_cast<int>(width) << " " << -static_cast<int>(height) << " " << static_cast<int>(width) << " " << static_cast<int>(height) << "\n";
	stream << "\n";
	stream << "%% Sites: " << sites.size() << "\n";
	stream << "\n";
	stream << "/radius { " << radius << " } def\n";
	stream << "\n";
	stream << "/p { radius 0 360 arc closepath fill } def\n";
	stream << "\n";
	stream << "0 0 0 setrgbcolor\n";
	stream << "\n";
	for (unsigned int i = 0; i < sites.size(); ++i) {
		stream << sites[i].location.x << " " << sites[i].location.y << " p\n";
	}
	for (int i = 0; i < optimizer.entries_.size(); i++) {
		switch (i) {
		case 1:
			stream << "1 0 0 setrgbcolor\n";
			break;
		case 2:
			stream << "0 1 0 setrgbcolor\n";
			break;
		case 3:
			stream << "0 0 1 setrgbcolor\n";
			break;
		case 4:
			stream << "1 0 1 setrgbcolor\n";
			break;
		case 5:
			stream << "1 1 0 setrgbcolor\n";
			break;
		case 0:
			stream << "0 1 1 setrgbcolor\n";
			break;
		}
		for (int j = 0; j < optimizer.entries_[i].points.size(); j++) {
			stream << optimizer.entries_[i].points[j].x << " " << optimizer.entries_[i].points[j].y << " p\n";
		}
	}
	stream << "\n";
	stream << "showpage\n";

	stream.close();
	return true;
}

void writeFile(std::stringstream &ss, std::vector<Point2> points) {
	for (int i = 0; i < static_cast<int>(points.size()); i++) {
		ss << points[i].x << " " << points[i].y << std::endl;
	}
}

int main(int, char*[]) {
	int     NUMBER_SITES = 56;					//վ������ 
	int     NUMBER_POINTS = 24 * NUMBER_SITES;	//������ ������ 
	const double  TORUS_SIZE = 10;					//���Ĵ�С 
	const bool    CONSTANT_DENSITY = true;					//�ܶ��Ƿ�㶨 
	const bool    CENTROIDAL = false;					//�Ƿ����site
	const bool    RESULT_PRINT = false;					//�Ƿ��ӡ��� 
	const bool    RESULT_FILE = true;					//����ļ� 
	const char*   RESULT_FILENAME = "result.eps";
	const double  RESULT_RADIUS = 0.05;						//����뾶 

	typedef Optimizer<Site<Point2>, Point2, MetricEuclidean2> Optimizer;//�Ż����� 

	// intializing the underlying discrete space����ʼ��Ǳ�ڵ���ɢ�ռ� 
	Point2::List points;//List<Point2> points 
//	if (CONSTANT_DENSITY) {//�ܶȺ㶨 
//		constant_regular_density(points, NUMBER_POINTS, TORUS_SIZE);
//	}
//	else {
//		nonconstant_density(points, NUMBER_POINTS, TORUS_SIZE);
//	}

	// initializing the Voronoi sites with equal capacity����ʼ��������ȵ�Voronoiվ�� 
//	unsigned int overallCapacity = static_cast<int>(points.size());//���е������=����վ��������� 
	Site<Point2>::List sites;
//	for (int i = 0; i < NUMBER_SITES; ++i) {//����վ������ 
//		double x = static_cast<double>(rand() % RAND_MAX) / RAND_MAX * TORUS_SIZE;
//		double y = static_cast<double>(rand() % RAND_MAX) / RAND_MAX * TORUS_SIZE;
//		int capacity = overallCapacity / (NUMBER_SITES - i);
//		overallCapacity -= capacity;
//		sites.push_back(Site<Point2>(i, capacity, Point2(x, y)));
//	}

	readFile(NUMBER_SITES, NUMBER_POINTS, points, sites);

	clock_t start = clock();

	// initializing the CCVT
	clock_t startInitialization = clock();
	printf("initialization...");
	Optimizer optimizer;
	MetricEuclidean2 metric;
	optimizer.initialize(sites, points, metric);
	printf("done\n");
	clock_t endInitialization = clock();

	// optimization
	int iteration = 0;
	bool stable;
	do {
		printf("iteration %d...", ++iteration);
		stable = optimizer.optimize(CENTROIDAL);
		printf("swap %d times,", optimizer.getSwapNumber());
		printf("done\n");
	} while (!stable);

	clock_t end = clock();
	printf("The state has reached a stable\n");
	
	int num_swap = 0;
	int mm = 0;
	optimizer.initializeGird(100,10.0/100);//格子个数、每个格子边长
	do {
		optimizer.updatePointFather();
		iteration = 0;
		optimizer.clear2();
		optimizer.updatePointLocation();
		optimizer.checkOrder();
		clock_t sss = clock();
		int all_swap = 0;
		printf("The status changes for the %d time.", ++mm);
		do {
			//printf("iteration %d...", ++iteration);
			stable = optimizer.myOptimize();
			num_swap = optimizer.getSwapNumber();
			all_swap += num_swap;
			//printf("swap %d times,",num_swap);
			//printf("done\n");
		} while (!stable);
		//std::cout << "OK" << std::endl;
		num_swap = optimizer.getSwapNumber();
		clock_t eee = clock();
		printf("The state has reached a stable,the all num_swap is %d,the time is %.3f\n", all_swap,static_cast<double>(eee-sss) / CLOCKS_PER_SEC);
		if (mm == 20) {
			break;
		}
		//std::system("pause");
	} while (num_swap == 0);
	
	const Site<Point2>::Vector& result = optimizer.sites();

	std::stringstream ss;
	ss << NUMBER_SITES << " " << NUMBER_POINTS << std::endl;

	for (int i = 0; i < NUMBER_SITES; i++) {
		ss << optimizer.x_of_sites(i) << " " << optimizer.y_of_sites(i) << std::endl;
	}

	for (int i = 0; i < NUMBER_SITES; i++) {
		const std::vector<Point2> p = optimizer.entryPoint(i);
		//printf("%d\n", p.size());
		for (int j = 0; j < p.size(); j++) {
			ss << p[j].x << " " << p[j].y << std::endl;
		}
	}
	oFile << ss.str();
	ss.clear();

	// writing the Voronoi sites to console
	if (RESULT_PRINT) {
		printf("\nresult:\n");
		for (unsigned int i = 0; i < result.size(); ++i) {
			printf("site %d: %f, %f\n", result[i].id, result[i].location.x, result[i].location.y);
		}
	}

	printf("\ninitialization time: %.3f sec\n", static_cast<double>(endInitialization - startInitialization) / CLOCKS_PER_SEC);
	printf("computation time: %.3f sec\n", static_cast<double>(end - start) / CLOCKS_PER_SEC);

	// writing the Voronoi sites to EPS file
	if (RESULT_FILE) {
		if (save_eps(RESULT_FILENAME, result, TORUS_SIZE, TORUS_SIZE, RESULT_RADIUS,optimizer)) {
			printf("\nresult saved in '%s'\n", RESULT_FILENAME);
		}
		else {
			printf("\nresult could not be saved in '%s'\n", RESULT_FILENAME);
		}
	}

	printf("\n");

	return 0;
}