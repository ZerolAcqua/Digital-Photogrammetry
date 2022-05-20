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
	//��������img����ȡͼ��
	Mat imgL, imgR;
	//�ҶȻ���ͼ��
	Mat grayImgL, grayImgR;
	Mat temp;
	//��ȡ��������
	vector<Point2i>featPoints;
	//ƥ���
	vector<pairedPoint2i>matchedPoints;
	Mat comparation;

#ifndef MY_DEBUG
	/*------------------------------------------��һ��Ӱ��--------------------------------------------*/
	imgL = imread("data/����Ӱ���������1.jpg");
	imgR = imread("data/����Ӱ���������2.jpg");
	if (imgL.empty()|| imgR.empty())
	{
		return 0;
	}

	//��Ϊ�Ҷ�ͼ��
	cvtColor(imgL, grayImgL, CV_BGR2GRAY);
	cvtColor(imgR, grayImgR, CV_BGR2GRAY);
	//��ȡ������
	MoravecPoints(grayImgL, featPoints, 15, 1000, 5);

	//����������Ƭ����ȡ��������
	imgL.copyTo(temp);
	for (auto it = featPoints.begin(); it != featPoints.end(); it++)
	{
		circle(temp, *it, 5, Scalar(0, 0, 255), 1, FILLED);
	}
	imwrite("Moravec.jpg", temp);

	//�������ϵ��ƥ��
	correlationCoefMatching(grayImgL, grayImgR, featPoints, matchedPoints, 10, 0.90);

	//��С����ƥ��
	vector<pairedPoint2f> refinedMatchedPoints;
	LeastSquaresMatching(grayImgL, grayImgR, matchedPoints, refinedMatchedPoints);

	//����ƥ���
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
	/*------------------------------------------�ڶ���Ӱ��--------------------------------------------*/
	imgL = imread("data/�Ϳշ�ͧӰ���������1.jpg");
	imgR = imread("data/�Ϳշ�ͧӰ���������2.jpg");
	if (imgL.empty() || imgR.empty())
	{
		return 0;
	}

	//��ʼ�Ķ�ͬ����ȷ����ŵ�͸�ӱ任����
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


	//����͸�ӱ任���������б任
	H = getPerspectiveTransform(src, dst);
	warpPerspective(imgR, correctedImgR, H, imgR.size());

	//��Ϊ�Ҷ�ͼ��
	cvtColor(imgL, grayImgL, CV_BGR2GRAY);
	cvtColor(correctedImgR, grayImgR, CV_BGR2GRAY);
	//��ȡ������
	MoravecPoints(grayImgL, featPoints, 15, 100, 10);

	//����������Ƭ����ȡ��������
	imgL.copyTo(temp);
	for (auto it = featPoints.begin(); it != featPoints.end(); it++)
	{
		circle(temp, *it, 5, Scalar(0, 0, 255), 1, FILLED);
	}
	imwrite("Moravec.jpg", temp);

	//�������ϵ��ƥ��
	correlationCoefMatching(grayImgL, grayImgR, featPoints, matchedPoints, 15, 0.95);
	
	//��ر任ǰ��ƥ�������
	vector<pairedPoint2i> matchedOriPoints;		
	inverseH = getPerspectiveTransform(dst, src);
	transToOriCoor(inverseH, matchedPoints, matchedOriPoints);	
	
	//��С����ƥ��
	cvtColor(imgR, grayImgR, CV_BGR2GRAY);
	vector<pairedPoint2f> refinedMatchedPoints;
	LeastSquaresMatching(grayImgL, grayImgR, matchedOriPoints, refinedMatchedPoints, 9, 13, 0.97);

	//����ƥ���
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

