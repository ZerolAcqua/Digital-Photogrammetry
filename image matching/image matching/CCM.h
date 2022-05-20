#pragma once
#include <opencv2/core/core.hpp>   
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/imgproc/types_c.h>

/**
*	相关系数匹配的一些函数
*	Author:@Hu Lechen
*/

using namespace std;
using namespace cv;

typedef pair<Point2i, Point2i> pairedPoint2i;
typedef pair<Point2f, Point2f> pairedPoint2f;

/**
*	@brief
*	根据左影像的特征点匹配两张影像
*
*	@param imgL					左影像
*	@param imgR					右影像
*	@param combinedImg			拼合的影像
*	@param pairedPoints			匹配点的坐标
*
*	@return 运行状态:			0运行正常|	-1左右图像大小不匹配
*/
int showMatchedPoints(const Mat& imgL, const Mat& imgR, Mat& combinedImg, const vector<pairedPoint2i>& pairedPoints);

/**
*	@brief
*	根据左影像的特征点匹配两张影像
*
*	@param imgL					左影像
*	@param imgR					右影像
*	@param combinedImg			拼合的影像
*	@param pairedPoints			匹配点的坐标
*
*	@return 运行状态:			0运行正常|	-1左右图像大小不匹配
*/
int showMatchedPoints(const Mat& imgL, const Mat& imgR, Mat& combinedImg, const vector<pairedPoint2f>& pairedPoints);

/**
*	@brief 
*	提取影像的特征点，并返回特征点的坐标
*
*	@param img					灰度影像
*	@param points				特征点的坐标
*	@param intereWindowSize		兴趣值窗口的大小（实际计算时的窗口大小略有不同）
*	@param threshold			兴趣值经验阈值
*	@param depressWindowSize	抑制局部非最大窗口大小，取0则不进行抑制
* 
*	@return 运行状态:			0运行正常|	-1图像不是CV_8UC1类型灰度图像|	-2窗口不是正方形 
*/
int MoravecPoints(const Mat& img, vector<Point2i>&points,int intereWindowSize=5,int threshold=0,int depressWindowSize = 5);


/**
*	@brief
*	计算给定窗口到右影像的最大相关系数
*
*	@param imgLWindow			目标影像正方形窗口
*	@param imgR					右影像
*	@param points				匹配点在右影像中的坐标
*	@param maxR					相关系数
*
*	@return 运行状态:			0运行正常|	-1图像不是CV_8UC1类型灰度图像|	-2窗口不是正方形
*/
int singleWindowMatching(const Mat& imgLWindow, const Mat& imgR, Point2i& point, double& maxR);

/**
*	@brief
*	计算给定窗口到右影像的最大相关系数
*
*	@param imgLWindow			目标影像正方形窗口
*	@param imgR					右影像
*	@param maxR					相关系数
*
*	@return 运行状态:			0运行正常|	-1图像不是CV_8UC1类型灰度图像|	-2窗口不是正方形
*/
int singleWindowMatching(const Mat& imgLWindow, const Mat& imgR, double& maxR);

/**
*	@brief
*	根据左影像的特征点匹配两张影像
*
*	@param imgL					左灰度影像
*	@param imgR					右灰度影像
*	@param points				特征点的坐标
*	@param pairedPoints				特征点对坐标
*	@param targetWindowSize		目标窗口大小（实际计算时的窗口大小略有不同）
*
*	@return 运行状态:			0运行正常|	-1图像不是CV_8UC1类型灰度图像|	-2左右图像大小不匹配
*/
int correlationCoefMatching(const Mat& imgL, const Mat& imgR, const vector<Point2i>& points, vector<pairedPoint2i>& pairedPoints,
	int targetWindowSize = 5, double threshold = 0.85);

/**
*	@brief
*	将透视变换后的匹配坐标，恢复到变换之前
*	@param inverseH				从变换后到变换前的透视变换矩阵
*	@param pairedPoints			变换后的坐标
*	@param matchedNewPoints		变换前的坐标
*	@return 运行状态:			0运行正常
*/
int transToOriCoor(const Mat& inverseH, vector<pairedPoint2i>& matchedPoints, vector<pairedPoint2i>& matchedOriPoints);