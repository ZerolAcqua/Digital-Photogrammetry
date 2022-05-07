#include <stdio.h>
#include <iostream>

#include <opencv2/core/core.hpp>   
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/imgproc/types_c.h>

#include "CCM.h"
#include "LSM.h"

#define MY_DEBUG

using namespace std;
using namespace cv;

int main()
{
	//创建矩阵img并读取图像
	Mat imgL, imgR;
	//灰度化的图像
	Mat grayImgL, grayImgR;
	Mat temp;
	//提取的特征点
	vector<Point2i>featPoints;
	//匹配点
	vector<pairedPoint2i>matchedPoints;
	Mat comparation;

#ifndef MY_DEBUG
	/*------------------------------------------第一对影像--------------------------------------------*/
	imgL = imread("data/航空影像立体像对1.jpg");
	imgR = imread("data/航空影像立体像对2.jpg");
	if (imgL.empty()|| imgR.empty())
	{
		return 0;
	}

	//化为灰度图像
	cvtColor(imgL, grayImgL, CV_BGR2GRAY);
	cvtColor(imgR, grayImgR, CV_BGR2GRAY);
	//提取特征点
	MoravecPoints(grayImgL, featPoints, 15, 1000, 5);

	//绘制在左像片上提取的特征点
	imgL.copyTo(temp);
	for (auto it = featPoints.begin(); it != featPoints.end(); it++)
	{
		circle(temp, *it, 5, Scalar(0, 0, 255), 1, FILLED);
	}
	imwrite("Moravec.jpg", temp);

	//进行相关系数匹配
	correlationCoefMatching(grayImgL, grayImgR, featPoints, matchedPoints, 10, 0.90);

	//最小二乘匹配
	vector<pairedPoint2f> refinedMatchedPoints;
	LeastSquaresMatching(grayImgL, grayImgR, matchedPoints, refinedMatchedPoints);

	//绘制匹配点
	showMatchedPoints(imgL, imgR, comparation, refinedMatchedPoints);
	imshow("match", comparation);
	imwrite("match.jpg", comparation);

	cout << endl;
	for (auto it = refinedMatchedPoints.begin(); it < refinedMatchedPoints.end(); it++)
	{
		cout << it->first << " " << it->second << endl;
}

	waitKey(0);


#endif // !MY_DEBUG

#ifdef MY_DEBUG
	/*------------------------------------------第二对影像--------------------------------------------*/
	imgL = imread("data/低空飞艇影像立体像对1.jpg");
	imgR = imread("data/低空飞艇影像立体像对2.jpg");
	if (imgL.empty() || imgR.empty())
	{
		return 0;
	}

	//初始四对同名点确定大概的透视变换参数
	vector<Point2f> src;
	vector<Point2f> dst;

	src.push_back(Point2f(250, 655));
	dst.push_back(Point2f(380, 436));

	src.push_back(Point2f(145, 337));
	dst.push_back(Point2f(334, 165));

	src.push_back(Point2f(108, 534));
	dst.push_back(Point2f(288, 317));

	src.push_back(Point2f(382, 531));
	dst.push_back(Point2f(521, 382));

	Mat correctedImgR;
	Mat H;
	Mat inverseH;


	//计算透视变换参数并进行变换
	H = getPerspectiveTransform(src, dst);
	warpPerspective(imgR, correctedImgR, H, imgR.size());

	//化为灰度图像
	cvtColor(imgL, grayImgL, CV_BGR2GRAY);
	cvtColor(correctedImgR, grayImgR, CV_BGR2GRAY);
	//提取特征点
	MoravecPoints(grayImgL, featPoints, 15, 100, 10);

	//绘制在左像片上提取的特征点
	imgL.copyTo(temp);
	for (auto it = featPoints.begin(); it != featPoints.end(); it++)
	{
		circle(temp, *it, 5, Scalar(0, 0, 255), 1, FILLED);
	}
	imwrite("Moravec.jpg", temp);

	//进行相关系数匹配
	correlationCoefMatching(grayImgL, grayImgR, featPoints, matchedPoints, 15, 0.95);
	
	//变回变换前的匹配点坐标
	vector<pairedPoint2i> matchedOriPoints;		
	inverseH = getPerspectiveTransform(dst, src);
	transToOriCoor(inverseH, matchedPoints, matchedOriPoints);	
	
	//最小二乘匹配
	cvtColor(imgR, grayImgR, CV_BGR2GRAY);
	vector<pairedPoint2f> refinedMatchedPoints;
	LeastSquaresMatching(grayImgL, grayImgR, matchedOriPoints, refinedMatchedPoints, 9, 13, 0.97);

	//绘制匹配点
	showMatchedPoints(imgL, imgR, comparation, refinedMatchedPoints);
	imshow("match", comparation);
	imwrite("match.jpg", comparation);

	cout << endl;
	for (auto it = refinedMatchedPoints.begin(); it < refinedMatchedPoints.end(); it++)
	{
		cout << it->first << " " << it->second << endl;
	}

	waitKey(0);



#endif // MY_DEBUG

}

