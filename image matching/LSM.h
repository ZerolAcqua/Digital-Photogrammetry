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
*	根据已有的匹配数据，进行最小两乘匹配,坐标是以相关系数匹配点为中心的相对坐标
*
*	@param imgL					左灰度影像
*	@param imgR					右灰度影像
*	@param pairedPoints			匹配的点对坐标
*	@param targetWindowSize		目标窗口大小
*	@param searchWindowSize		搜索窗口大小
* 	@param threshold			相关系数退出条件
*	@param iterationTimes		最大迭代次数
*
*	@return 运行状态:			0运行正常|	-1图像不是CV_8UC1类型灰度图像|	-2窗口不是正方形|	-3窗口太小无法进行最小二乘平差
*/
int LeastSquaresMatching(const Mat& imgL, const Mat& imgR, vector<pairedPoint2i>& pairedPoints, vector<pairedPoint2f>& refinedPairedPoints,
	int targetWindowSize = 5, int searchWindowSize = 13, double threshold = 0.98, int iterationTimes = 10);