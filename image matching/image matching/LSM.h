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
*	�������е�ƥ�����ݣ�������С����ƥ��,�����������ϵ��ƥ���Ϊ���ĵ��������
*
*	@param imgL					��Ҷ�Ӱ��
*	@param imgR					�һҶ�Ӱ��
*	@param pairedPoints			ƥ��ĵ������
*	@param targetWindowSize		Ŀ�괰�ڴ�С
*	@param searchWindowSize		�������ڴ�С
* 	@param threshold			���ϵ���˳�����
*	@param iterationTimes		����������
*
*	@return ����״̬:			0��������|	-1ͼ����CV_8UC1���ͻҶ�ͼ��|	-2���ڲ���������|	-3����̫С�޷�������С����ƽ��
*/
int LeastSquaresMatching(const Mat& imgL, const Mat& imgR, vector<pairedPoint2i>& pairedPoints, vector<pairedPoint2f>& refinedPairedPoints,
	int targetWindowSize = 5, int searchWindowSize = 13, double threshold = 0.98, int iterationTimes = 10);