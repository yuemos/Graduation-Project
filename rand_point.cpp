#include <cstdio>
#include <cstdlib>
#include <time.h>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>

using namespace std;

int main(){
	srand((unsigned)time(NULL));
	int number_of_sites = 0,number_of_points=0;
	scanf("%d%d",&number_of_sites,&number_of_points);
	
	ofstream oFile("D:/workspace/CCVT/in.txt") ;
	stringstream ss;
	
	ss<<number_of_sites<<" "<<number_of_points<<endl; 
	
	for(int i = 0;i<number_of_sites;i++){
		double x,y;
		x=rand()%1000 / 100.0;
		y=rand()%1000 / 100.0;
		ss<<x<<" "<<y<<endl;
	}
	
	for(int i = 0;i<number_of_sites*number_of_points;i++){
		int angle;
		double x,y,w,radiu;
		x=rand()%1000 / 100.0;//生成圆心x 
		y=rand()%1000 / 100.0;//生成圆心y
		ss<<x<<" "<<y<<endl;
		radiu=rand()%100/10.0;//生成半径 
		angle=rand()%360;//生成初始角度 
		ss<<radiu<<" "<<angle<<endl;
		w=rand()%20/10.0;//生成半径 
		ss<<w<<" "<<endl; 
	}
	
	oFile<<ss.str();
	
	return 0;
}
