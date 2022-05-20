#include "CCM.h"

int showMatchedPoints(const Mat& imgL, const Mat& imgR, Mat& combinedImg, const vector<pairedPoint2i>& pairedPoints)
{
	Scalar temp;
	int i = 1;
	if (imgL.size != imgR.size)
	{
		return -1;
	}
	combinedImg = Mat(imgL.rows, imgL.cols * 2, CV_8UC3);
	imgL.copyTo(combinedImg(Rect(Point2i(0, 0), Size(imgL.cols, imgL.rows))));
	imgR.copyTo(combinedImg(Rect(Point2i(imgL.cols, 0), Size(imgL.cols, imgL.rows))));

	for (auto it = pairedPoints.begin(); it != pairedPoints.end(); it++)
	{
		temp = Scalar((i * 10 * 10 * 10) % 205 + 50, (i * 10 * 10) % 205 + 50, 255);
		line(combinedImg, it->first, it->second + Point2i(imgL.cols, 0), temp, 1);
		circle(combinedImg, it->first, 5, temp, 1, FILLED);
		circle(combinedImg, it->second + Point2i(imgL.cols, 0), 5, temp, 1, FILLED);
		i++;
	}
	return 0;
}

int showMatchedPoints(const Mat& imgL, const Mat& imgR, Mat& combinedImg, const vector<pairedPoint2f>& pairedPoints)
{
	Scalar temp;
	int i = 1;
	if (imgL.size != imgR.size)
	{
		return -1;
	}
	combinedImg = Mat(imgL.rows, imgL.cols * 2, CV_8UC3);
	imgL.copyTo(combinedImg(Rect(Point2i(0, 0), Size(imgL.cols, imgL.rows))));
	imgR.copyTo(combinedImg(Rect(Point2i(imgL.cols, 0), Size(imgL.cols, imgL.rows))));

	for (auto it = pairedPoints.begin(); it != pairedPoints.end(); it++)
	{
		temp = Scalar((i * 10 * 10 * 10) % 205 + 50, (i * 10 * 10) % 205 + 50, 255);
		line(combinedImg, it->first, it->second + Point2f(imgL.cols, 0), temp, 1);
		circle(combinedImg, it->first, 5, temp, 1, FILLED);
		circle(combinedImg, it->second + Point2f(imgL.cols, 0), 5, temp, 1, FILLED);
		i++;
	}
	return 0;
}


int MoravecPoints(const Mat& img, vector<Point2i>& points, int intereWindowSize, int threshold, int depressWindowSize)
{
	if (img.type() != CV_8UC1)
	{
		return -1;
	}

	points.clear();
	//计算各像元兴趣值
	int k = intereWindowSize / 2;
	int v1, v2, v3, v4 = 0;	//四个方向的灰度差平方和
	double iv = 0;			//最小的平方和作为兴趣值
	Mat intereMat(img.rows, img.cols, CV_32SC1, Scalar(0));	//存储兴趣值

	//两边各留k个像素
	for (int row = k; row < img.rows - k; row++)
	{
		for (int col = k; col < img.cols - k; col++)
		{
			v1 = v2 = v3 = v4 = 0;
			for (int i = -k; i < k; i++)
			{
				v1 += pow(img.at<uchar>(row, col + i) - img.at<uchar>(row, col + i + 1), 2);
				v2 += pow(img.at<uchar>(row + 1, col + i) - img.at<uchar>(row + i + 1, col + i + 1), 2);
				v3 += pow(img.at<uchar>(row + i, col) - img.at<uchar>(row + i + 1, col), 2);
				v4 += pow(img.at<uchar>(row - i, col + i) - img.at<uchar>(row - i - 1, col + i + 1), 2);
			}
			//改变，除以（2*k）进行归一化，效果如何？
			iv = MIN(v1, MIN(v2, MIN(v3, v4))) / (2 * k);
			//小于阈值取0
			intereMat.at<int>(row, col) = iv < threshold ? 0 : iv;
		}
	}
	//抑制局部非最大
	int maxIntere = 0;		//最大兴趣值
	int curIntere = 0;		//当前兴趣值
	bool flag = false;		//是否找到了最大值标识

	if (depressWindowSize <= 0)
	{
		//不做抑制
		for (int i = 0; i < img.rows; i++)
		{
			for (int j = 0; j < img.cols; j++)
			{
				curIntere = intereMat.at<int>(i, j);
				if (curIntere == 0)
					continue;
				points.push_back(Point2f(j, i));

			}
		}
	}
	else {
		//做抑制
		for (int i = 0; i < img.rows / depressWindowSize; i++)
		{
			for (int j = 0; j < img.cols / depressWindowSize; j++)
			{
				maxIntere = 0;
				curIntere = 0;
				flag = false;
				for (int row = 0; row < depressWindowSize && flag == false; row++)
				{
					//找到最大值
					for (int col = 0; col < depressWindowSize; col++)
					{
						curIntere = intereMat.at<int>(i * depressWindowSize + row, j * depressWindowSize + col);
						if (curIntere == 0)
							continue;
						if (curIntere > maxIntere)
							maxIntere = curIntere;
					}
					//去掉非最大值
					for (int col = 0; col < depressWindowSize; col++)
					{
						curIntere = intereMat.at<int>(i * depressWindowSize + row, j * depressWindowSize + col);
						if (curIntere == maxIntere && maxIntere != 0)
						{
							//一个区域只选一个最大值
							points.push_back(Point2f(j * depressWindowSize + col, i * depressWindowSize + row));
							flag = true;
							break;
						}
					}
				}
			}
		}
	}
	return 0;
}

int singleWindowMatching(const Mat& imgLWindow, const Mat& imgR, Point2i& point, double& maxR)
{
	//需要是灰度
	if (imgLWindow.type() != CV_8UC1 || imgR.type() != CV_8UC1)
	{
		return -1;
	}
	//窗口得是正方形
	if (imgLWindow.rows != imgLWindow.cols)
	{
		return -2;
	}

	int windowSize = imgLWindow.rows;
	int rows = imgR.rows;
	int cols = imgR.cols;

	Mat tempWindow;

	double R = 0;			//相关系数
	double Rnum = 0;		//相关系数的分子
	double Rden = 1;		//相关系数的分母
	maxR = -2;				//当前最大的相关系数
	Point2i matchedPoint;	//当前最相关的匹配点

	int Sg = 0;		//Sum(g)
	int Sg_ = 0;	//Sum(g')
	int Sg2 = 0;	//Sum(g^2)
	int Sg_2 = 0;	//Sum(g'^2)
	int Sgg_ = 0;	//Sum(g*g')

	//目标窗口与相邻的两块窗口进行匹配时，会有一定的重复计算，
	//以下变量是为了减少计算量
	int preSg_ = 0;		//(i,j-1)窗口的第一列的Sum(g')
	int preSg_2 = 0;	//(i,j-1)窗口的第一列的Sum(g'^2)
	int curSg_ = 0;		//(i,j)窗口的最后一列的Sum(g')
	int curSg_2 = 0;	//(i,j)窗口的最后一列的Sum(g'^2)
	int tempPreSg_;
	int tempPreSg_2 = 0;

	//目标窗口的灰度是固定的，先进行计算
	for (int i = 0; i < windowSize; i++)
	{
		for (int j = 0; j < windowSize; j++)
		{
			Sg += imgLWindow.at<uchar>(i, j);
			Sg2 += pow(imgLWindow.at<uchar>(i, j), 2);
		}
	}

	for (int row = 0; row <= rows - windowSize; row++)
	{
		//每行的第一个窗口需要全部计算
		tempWindow = imgR(Rect(Point2i(0, row), Size(windowSize, windowSize)));
		Sg_ = 0;
		Sg_2 = 0;
		Sgg_ = 0;
		preSg_ = 0;
		preSg_2 = 0;
		for (int i = 0; i < windowSize; i++)
		{
			preSg_ += tempWindow.at<uchar>(i, 0);
			preSg_2 += pow(tempWindow.at<uchar>(i, 0), 2);
			for (int j = 0; j < windowSize; j++)
			{
				Sg_ += tempWindow.at<uchar>(i, j);
				Sg_2 += pow(tempWindow.at<uchar>(i, j), 2);
				Sgg_ += imgLWindow.at<uchar>(i, j) * tempWindow.at<uchar>(i, j);
			}
		}
		Rnum = Sgg_ - double(Sg * Sg_) / (windowSize * windowSize);
		Rden = sqrt((Sg2 - double(Sg * Sg) / (windowSize * windowSize)) * ((Sg_2 - double(Sg_ * Sg_) / (windowSize * windowSize))));
		R = Rnum / Rden;
		if (R > maxR)
		{
			maxR = R;
			point = Point(windowSize / 2, row + windowSize / 2);
		}

		//同一行之后的窗口
		for (int col = 1; col <= cols - windowSize; col++)
		{
			tempWindow = imgR(Rect(Point2i(col, row), Size(windowSize, windowSize)));
			Sgg_ = 0;
			tempPreSg_ = 0;
			tempPreSg_2 = 0;
			curSg_ = 0;
			curSg_2 = 0;
			for (int i = 0; i < windowSize; i++)
			{
				tempPreSg_ += tempWindow.at<uchar>(i, 0);
				tempPreSg_2 += pow(tempWindow.at<uchar>(i, 0), 2);
				curSg_ += tempWindow.at<uchar>(i, windowSize - 1);
				curSg_2 += pow(tempWindow.at<uchar>(i, windowSize - 1), 2);
				for (int j = 0; j < windowSize; j++)
				{
					Sgg_ += imgLWindow.at<uchar>(i, j) * tempWindow.at<uchar>(i, j);
				}
			}
			Sg_ += curSg_ - preSg_;
			Sg_2 += curSg_2 - preSg_2;
			preSg_ = tempPreSg_;
			preSg_2 = tempPreSg_2;

			Rnum = Sgg_ - double(Sg * Sg_) / (windowSize* windowSize);
			Rden = sqrt((Sg2 - double(Sg * Sg) / (windowSize * windowSize)) * ((Sg_2 - double(Sg_ * Sg_) / (windowSize * windowSize))));
			R = Rnum / Rden;
			if (R > maxR)
			{
				maxR = R;
				point = Point(col + windowSize / 2, row + windowSize / 2);
			}
		}
	}
	return 0;
}

int singleWindowMatching(const Mat& imgLWindow, const Mat& imgR, double& maxR)
{
	//需要是灰度
	if (imgLWindow.type() != CV_8UC1 || imgR.type() != CV_8UC1)
	{
		return -1;
	}
	//窗口得是正方形
	if (imgLWindow.rows != imgLWindow.cols)
	{
		return -2;
	}

	int windowSize = imgLWindow.rows;
	int rows = imgR.rows;
	int cols = imgR.cols;

	Mat tempWindow;

	double R = 0;			//相关系数
	double Rnum = 0;		//相关系数的分子
	double Rden = 1;		//相关系数的分母
	maxR = -2;				//当前最大的相关系数
	Point2i matchedPoint;	//当前最相关的匹配点

	int Sg = 0;		//Sum(g)
	int Sg_ = 0;	//Sum(g')
	int Sg2 = 0;	//Sum(g^2)
	int Sg_2 = 0;	//Sum(g'^2)
	int Sgg_ = 0;	//Sum(g*g')

	//目标窗口与相邻的两块窗口进行匹配时，会有一定的重复计算，
	//以下变量是为了减少计算量
	int preSg_ = 0;		//(i,j-1)窗口的第一列的Sum(g')
	int preSg_2 = 0;	//(i,j-1)窗口的第一列的Sum(g'^2)
	int curSg_ = 0;		//(i,j)窗口的最后一列的Sum(g')
	int curSg_2 = 0;	//(i,j)窗口的最后一列的Sum(g'^2)
	int tempPreSg_;
	int tempPreSg_2 = 0;

	//目标窗口的灰度是固定的，先进行计算
	for (int i = 0; i < windowSize; i++)
	{
		for (int j = 0; j < windowSize; j++)
		{
			Sg += imgLWindow.at<uchar>(i, j);
			Sg2 += pow(imgLWindow.at<uchar>(i, j), 2);
		}
	}

	for (int row = 0; row <= rows - windowSize; row++)
	{
		//每行的第一个窗口需要全部计算
		tempWindow = imgR(Rect(Point2i(0, row), Size(windowSize, windowSize)));
		Sg_ = 0;
		Sg_2 = 0;
		Sgg_ = 0;
		preSg_ = 0;
		preSg_2 = 0;
		for (int i = 0; i < windowSize; i++)
		{
			preSg_ += tempWindow.at<uchar>(i, 0);
			preSg_2 += pow(tempWindow.at<uchar>(i, 0), 2);
			for (int j = 0; j < windowSize; j++)
			{
				Sg_ += tempWindow.at<uchar>(i, j);
				Sg_2 += pow(tempWindow.at<uchar>(i, j), 2);
				Sgg_ += imgLWindow.at<uchar>(i, j) * tempWindow.at<uchar>(i, j);
			}
		}
		Rnum = Sgg_ - double(Sg * Sg_) / (windowSize * windowSize);
		Rden = sqrt((Sg2 - double(Sg * Sg) / (windowSize * windowSize)) * ((Sg_2 - double(Sg_ * Sg_) / (windowSize * windowSize))));
		R = Rnum / Rden;
		if (R > maxR)
		{
			maxR = R;
		}

		//同一行之后的窗口
		for (int col = 1; col <= cols - windowSize; col++)
		{
			tempWindow = imgR(Rect(Point2i(col, row), Size(windowSize, windowSize)));
			Sgg_ = 0;
			tempPreSg_ = 0;
			tempPreSg_2 = 0;
			curSg_ = 0;
			curSg_2 = 0;
			for (int i = 0; i < windowSize; i++)
			{
				tempPreSg_ += tempWindow.at<uchar>(i, 0);
				tempPreSg_2 += pow(tempWindow.at<uchar>(i, 0), 2);
				curSg_ += tempWindow.at<uchar>(i, windowSize - 1);
				curSg_2 += pow(tempWindow.at<uchar>(i, windowSize - 1), 2);
				for (int j = 0; j < windowSize; j++)
				{
					Sgg_ += imgLWindow.at<uchar>(i, j) * tempWindow.at<uchar>(i, j);
				}
			}
			Sg_ += curSg_ - preSg_;
			Sg_2 += curSg_2 - preSg_2;
			preSg_ = tempPreSg_;
			preSg_2 = tempPreSg_2;

			Rnum = Sgg_ - double(Sg * Sg_) / (windowSize * windowSize);
			Rden = sqrt((Sg2 - double(Sg * Sg) / (windowSize * windowSize)) * ((Sg_2 - double(Sg_ * Sg_) / (windowSize * windowSize))));
			R = Rnum / Rden;
			if (R > maxR)
			{
				maxR = R;
			}
		}
	}
	return 0;
}

int correlationCoefMatching(const Mat& imgL, const Mat& imgR, const vector<Point2i>& points, vector<pairedPoint2i>& pairedPoints, int targetWindowSize, double threshold)
{
	//需要是灰度
	if (imgL.type() != CV_8UC1 || imgR.type() != CV_8UC1)
	{
		return -1;
	}
	//虽然图像的大小差别不大，基本没有影响，但为了稳定性，要求图像的大小相同
	if (imgL.size != imgR.size)
	{
		return -2;
	}
	pairedPoints.clear();
	int k = targetWindowSize / 2;
	int rows = imgL.rows;
	int cols = imgL.cols;

	Mat targetWindow;
	Mat tempWindow;
	double R = 0;			//相关系数
	Point2i matchedPoint;	//当前最相关的匹配点
	Point2i predict(0, 0);	//预测窗口位置的平移向量（即左右视差和上下视差）
	Mat searchWindow;		//搜索窗口
	int searchSize = 30;	//搜索窗口大小
	searchSize /= 2;
	int weight = 0;
	Rect inter;				//矩形交集

	for (auto itPoint = points.begin(); itPoint < points.end(); itPoint++)
	{
		//剔除靠近图像边缘的特征，一是窗口大小超出了图像范围，二是框幅式影像本身边缘的变形就较大
		if (itPoint->x < k || cols - itPoint->x < k || itPoint->y < k || rows - itPoint->y < k)
		{
			continue;
		}
		if (pairedPoints.size() > 3)
		{
			//用已匹配的点进行预测
			targetWindow = imgL(Rect(Point2i(itPoint->x - k, itPoint->y - k), Size(2 * k, 2 * k)));
			//并确定搜索区域
			if (itPoint->x - predict.x - searchSize < 0 || itPoint->y - predict.y - searchSize<0
				|| itPoint->x - predict.x + searchSize > cols || itPoint->y - predict.y + searchSize  >rows)
			{
				inter = Rect(Point2i(itPoint->x - k, itPoint->y - k), Size(2 * k, 2 * k)) & Rect(Point2i(0, 0), Size(imgL.cols, imgL.rows));
				if (inter.empty())
					continue;
				targetWindow = imgL(inter);
			}
			else
			{
				searchWindow = imgR(Rect(Point2i(itPoint->x - predict.x - searchSize, itPoint->y - predict.y - searchSize), Size(searchSize * 2, searchSize * 2)));
			}
			//对当前窗口进行匹配
			singleWindowMatching(targetWindow, searchWindow, matchedPoint, R);
			matchedPoint.x += itPoint->x - predict.x - searchSize;
			matchedPoint.y += itPoint->y - predict.y - searchSize;
			//相关系数太小则匹配不成功
			if (R > threshold)
			{
				predict = (predict * double(pairedPoints.size()) + *itPoint - matchedPoint) / double(pairedPoints.size() + 1);
				pairedPoints.push_back(pairedPoint2i(*itPoint, matchedPoint));
			}
		}
		else
		{
			targetWindow = imgL(Rect(Point2i(itPoint->x - k, itPoint->y - k), Size(2 * k, 2 * k)));
			//对当前窗口进行匹配
			singleWindowMatching(targetWindow, imgR, matchedPoint, R);
			//相关系数太小则匹配不成功
			if (R > threshold)
			{
				predict = (predict * int(pairedPoints.size()) + *itPoint - matchedPoint) / int(pairedPoints.size() + 1);
				pairedPoints.push_back(pairedPoint2i(*itPoint, matchedPoint));
			}
		}

	}
	return 0;
}

int transToOriCoor(const Mat& inverseH, vector<pairedPoint2i>& matchedPoints, vector<pairedPoint2i>& matchedOriPoints)
{
	double x, y;
	double x_, y_, z_;
	Point2i temp;
	matchedOriPoints.clear();
	inverseH.convertTo(inverseH, CV_64F);

	for (auto it = matchedPoints.begin(); it < matchedPoints.end(); it++)
	{
		x = it->second.x;
		y = it->second.y;
		x_ = inverseH.at<double>(0, 0) * x + inverseH.at<double>(0, 1) * y + inverseH.at<double>(0, 2);
		y_ = inverseH.at<double>(1, 0) * x + inverseH.at<double>(1, 1) * y + inverseH.at<double>(1, 2);
		z_ = inverseH.at<double>(2, 0) * x + inverseH.at<double>(2, 1) * y + inverseH.at<double>(2, 2);
		temp.x = x_ / z_;
		temp.y = y_ / z_;
		matchedOriPoints.push_back(pairedPoint2i(it->first, temp));
	}

	return 0;
}