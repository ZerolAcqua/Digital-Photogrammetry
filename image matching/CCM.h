#pragma once
#include <opencv2/core/core.hpp>   
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/imgproc/types_c.h>

/**
*	���ϵ��ƥ���һЩ����
*	Author:@Hu Lechen
*/

using namespace std;
using namespace cv;

typedef pair<Point2i, Point2i> pairedPoint2i;
typedef pair<Point2f, Point2f> pairedPoint2f;

/**
*	@brief
*	������Ӱ���������ƥ������Ӱ��
*
*	@param imgL					��Ӱ��
*	@param imgR					��Ӱ��
*	@param combinedImg			ƴ�ϵ�Ӱ��
*	@param pairedPoints			ƥ��������
*
*	@return ����״̬:			0��������|	-1����ͼ���С��ƥ��
*/
int showMatchedPoints(const Mat& imgL, const Mat& imgR, Mat& combinedImg, const vector<pairedPoint2i>& pairedPoints);

/**
*	@brief
*	������Ӱ���������ƥ������Ӱ��
*
*	@param imgL					��Ӱ��
*	@param imgR					��Ӱ��
*	@param combinedImg			ƴ�ϵ�Ӱ��
*	@param pairedPoints			ƥ��������
*
*	@return ����״̬:			0��������|	-1����ͼ���С��ƥ��
*/
int showMatchedPoints(const Mat& imgL, const Mat& imgR, Mat& combinedImg, const vector<pairedPoint2f>& pairedPoints);

/**
*	@brief 
*	��ȡӰ��������㣬�����������������
*
*	@param img					�Ҷ�Ӱ��
*	@param points				�����������
*	@param intereWindowSize		��Ȥֵ���ڵĴ�С��ʵ�ʼ���ʱ�Ĵ��ڴ�С���в�ͬ��
*	@param threshold			��Ȥֵ������ֵ
*	@param depressWindowSize	���ƾֲ�����󴰿ڴ�С��ȡ0�򲻽�������
* 
*	@return ����״̬:			0��������|	-1ͼ����CV_8UC1���ͻҶ�ͼ��|	-2���ڲ��������� 
*/
int MoravecPoints(const Mat& img, vector<Point2i>&points,int intereWindowSize=5,int threshold=0,int depressWindowSize = 5);


/**
*	@brief
*	����������ڵ���Ӱ���������ϵ��
*
*	@param imgLWindow			Ŀ��Ӱ�������δ���
*	@param imgR					��Ӱ��
*	@param points				ƥ�������Ӱ���е�����
*	@param maxR					���ϵ��
*
*	@return ����״̬:			0��������|	-1ͼ����CV_8UC1���ͻҶ�ͼ��|	-2���ڲ���������
*/
int singleWindowMatching(const Mat& imgLWindow, const Mat& imgR, Point2i& point, double& maxR);

/**
*	@brief
*	����������ڵ���Ӱ���������ϵ��
*
*	@param imgLWindow			Ŀ��Ӱ�������δ���
*	@param imgR					��Ӱ��
*	@param maxR					���ϵ��
*
*	@return ����״̬:			0��������|	-1ͼ����CV_8UC1���ͻҶ�ͼ��|	-2���ڲ���������
*/
int singleWindowMatching(const Mat& imgLWindow, const Mat& imgR, double& maxR);

/**
*	@brief
*	������Ӱ���������ƥ������Ӱ��
*
*	@param imgL					��Ҷ�Ӱ��
*	@param imgR					�һҶ�Ӱ��
*	@param points				�����������
*	@param pairedPoints				�����������
*	@param targetWindowSize		Ŀ�괰�ڴ�С��ʵ�ʼ���ʱ�Ĵ��ڴ�С���в�ͬ��
*
*	@return ����״̬:			0��������|	-1ͼ����CV_8UC1���ͻҶ�ͼ��|	-2����ͼ���С��ƥ��
*/
int correlationCoefMatching(const Mat& imgL, const Mat& imgR, const vector<Point2i>& points, vector<pairedPoint2i>& pairedPoints,
	int targetWindowSize = 5, double threshold = 0.85);

/**
*	@brief
*	��͸�ӱ任���ƥ�����꣬�ָ����任֮ǰ
*	@param inverseH				�ӱ任�󵽱任ǰ��͸�ӱ任����
*	@param pairedPoints			�任�������
*	@param matchedNewPoints		�任ǰ������
*	@return ����״̬:			0��������
*/
int transToOriCoor(const Mat& inverseH, vector<pairedPoint2i>& matchedPoints, vector<pairedPoint2i>& matchedOriPoints);