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
	//�������Ԫ��Ȥֵ
	int k = intereWindowSize / 2;
	int v1, v2, v3, v4 = 0;	//�ĸ�����ĻҶȲ�ƽ����
	double iv = 0;			//��С��ƽ������Ϊ��Ȥֵ
	Mat intereMat(img.rows, img.cols, CV_32SC1, Scalar(0));	//�洢��Ȥֵ

	//���߸���k������
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
			//�ı䣬���ԣ�2*k�����й�һ����Ч����Σ�
			iv = MIN(v1, MIN(v2, MIN(v3, v4))) / (2 * k);
			//С����ֵȡ0
			intereMat.at<int>(row, col) = iv < threshold ? 0 : iv;
		}
	}
	//���ƾֲ������
	int maxIntere = 0;		//�����Ȥֵ
	int curIntere = 0;		//��ǰ��Ȥֵ
	bool flag = false;		//�Ƿ��ҵ������ֵ��ʶ

	if (depressWindowSize <= 0)
	{
		//��������
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
		//������
		for (int i = 0; i < img.rows / depressWindowSize; i++)
		{
			for (int j = 0; j < img.cols / depressWindowSize; j++)
			{
				maxIntere = 0;
				curIntere = 0;
				flag = false;
				for (int row = 0; row < depressWindowSize && flag == false; row++)
				{
					//�ҵ����ֵ
					for (int col = 0; col < depressWindowSize; col++)
					{
						curIntere = intereMat.at<int>(i * depressWindowSize + row, j * depressWindowSize + col);
						if (curIntere == 0)
							continue;
						if (curIntere > maxIntere)
							maxIntere = curIntere;
					}
					//ȥ�������ֵ
					for (int col = 0; col < depressWindowSize; col++)
					{
						curIntere = intereMat.at<int>(i * depressWindowSize + row, j * depressWindowSize + col);
						if (curIntere == maxIntere && maxIntere != 0)
						{
							//һ������ֻѡһ�����ֵ
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
	//��Ҫ�ǻҶ�
	if (imgLWindow.type() != CV_8UC1 || imgR.type() != CV_8UC1)
	{
		return -1;
	}
	//���ڵ���������
	if (imgLWindow.rows != imgLWindow.cols)
	{
		return -2;
	}

	int windowSize = imgLWindow.rows;
	int rows = imgR.rows;
	int cols = imgR.cols;

	Mat tempWindow;

	double R = 0;			//���ϵ��
	double Rnum = 0;		//���ϵ���ķ���
	double Rden = 1;		//���ϵ���ķ�ĸ
	maxR = -2;				//��ǰ�������ϵ��
	Point2i matchedPoint;	//��ǰ����ص�ƥ���

	int Sg = 0;		//Sum(g)
	int Sg_ = 0;	//Sum(g')
	int Sg2 = 0;	//Sum(g^2)
	int Sg_2 = 0;	//Sum(g'^2)
	int Sgg_ = 0;	//Sum(g*g')

	//Ŀ�괰�������ڵ����鴰�ڽ���ƥ��ʱ������һ�����ظ����㣬
	//���±�����Ϊ�˼��ټ�����
	int preSg_ = 0;		//(i,j-1)���ڵĵ�һ�е�Sum(g')
	int preSg_2 = 0;	//(i,j-1)���ڵĵ�һ�е�Sum(g'^2)
	int curSg_ = 0;		//(i,j)���ڵ����һ�е�Sum(g')
	int curSg_2 = 0;	//(i,j)���ڵ����һ�е�Sum(g'^2)
	int tempPreSg_;
	int tempPreSg_2 = 0;

	//Ŀ�괰�ڵĻҶ��ǹ̶��ģ��Ƚ��м���
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
		//ÿ�еĵ�һ��������Ҫȫ������
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

		//ͬһ��֮��Ĵ���
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
	//��Ҫ�ǻҶ�
	if (imgLWindow.type() != CV_8UC1 || imgR.type() != CV_8UC1)
	{
		return -1;
	}
	//���ڵ���������
	if (imgLWindow.rows != imgLWindow.cols)
	{
		return -2;
	}

	int windowSize = imgLWindow.rows;
	int rows = imgR.rows;
	int cols = imgR.cols;

	Mat tempWindow;

	double R = 0;			//���ϵ��
	double Rnum = 0;		//���ϵ���ķ���
	double Rden = 1;		//���ϵ���ķ�ĸ
	maxR = -2;				//��ǰ�������ϵ��
	Point2i matchedPoint;	//��ǰ����ص�ƥ���

	int Sg = 0;		//Sum(g)
	int Sg_ = 0;	//Sum(g')
	int Sg2 = 0;	//Sum(g^2)
	int Sg_2 = 0;	//Sum(g'^2)
	int Sgg_ = 0;	//Sum(g*g')

	//Ŀ�괰�������ڵ����鴰�ڽ���ƥ��ʱ������һ�����ظ����㣬
	//���±�����Ϊ�˼��ټ�����
	int preSg_ = 0;		//(i,j-1)���ڵĵ�һ�е�Sum(g')
	int preSg_2 = 0;	//(i,j-1)���ڵĵ�һ�е�Sum(g'^2)
	int curSg_ = 0;		//(i,j)���ڵ����һ�е�Sum(g')
	int curSg_2 = 0;	//(i,j)���ڵ����һ�е�Sum(g'^2)
	int tempPreSg_;
	int tempPreSg_2 = 0;

	//Ŀ�괰�ڵĻҶ��ǹ̶��ģ��Ƚ��м���
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
		//ÿ�еĵ�һ��������Ҫȫ������
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

		//ͬһ��֮��Ĵ���
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
	//��Ҫ�ǻҶ�
	if (imgL.type() != CV_8UC1 || imgR.type() != CV_8UC1)
	{
		return -1;
	}
	//��Ȼͼ��Ĵ�С��𲻴󣬻���û��Ӱ�죬��Ϊ���ȶ��ԣ�Ҫ��ͼ��Ĵ�С��ͬ
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
	double R = 0;			//���ϵ��
	Point2i matchedPoint;	//��ǰ����ص�ƥ���
	Point2i predict(0, 0);	//Ԥ�ⴰ��λ�õ�ƽ���������������Ӳ�������Ӳ
	Mat searchWindow;		//��������
	int searchSize = 30;	//�������ڴ�С
	searchSize /= 2;
	int weight = 0;
	Rect inter;				//���ν���

	for (auto itPoint = points.begin(); itPoint < points.end(); itPoint++)
	{
		//�޳�����ͼ���Ե��������һ�Ǵ��ڴ�С������ͼ��Χ�����ǿ��ʽӰ�����Ե�ı��ξͽϴ�
		if (itPoint->x < k || cols - itPoint->x < k || itPoint->y < k || rows - itPoint->y < k)
		{
			continue;
		}
		if (pairedPoints.size() > 3)
		{
			//����ƥ��ĵ����Ԥ��
			targetWindow = imgL(Rect(Point2i(itPoint->x - k, itPoint->y - k), Size(2 * k, 2 * k)));
			//��ȷ����������
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
			//�Ե�ǰ���ڽ���ƥ��
			singleWindowMatching(targetWindow, searchWindow, matchedPoint, R);
			matchedPoint.x += itPoint->x - predict.x - searchSize;
			matchedPoint.y += itPoint->y - predict.y - searchSize;
			//���ϵ��̫С��ƥ�䲻�ɹ�
			if (R > threshold)
			{
				predict = (predict * double(pairedPoints.size()) + *itPoint - matchedPoint) / double(pairedPoints.size() + 1);
				pairedPoints.push_back(pairedPoint2i(*itPoint, matchedPoint));
			}
		}
		else
		{
			targetWindow = imgL(Rect(Point2i(itPoint->x - k, itPoint->y - k), Size(2 * k, 2 * k)));
			//�Ե�ǰ���ڽ���ƥ��
			singleWindowMatching(targetWindow, imgR, matchedPoint, R);
			//���ϵ��̫С��ƥ�䲻�ɹ�
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