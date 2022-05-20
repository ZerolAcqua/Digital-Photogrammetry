#include "LSM.h"
#include "CCM.h"

#include <iostream>

int LeastSquaresMatching(const Mat& imgL, const Mat& imgR, vector<pairedPoint2i>& pairedPoints,vector<pairedPoint2f>& refinedPairedPoints, int targetWindowSize, int searchWindowSize, double threshold, int iterationTimes)
{
	//��Ҫ�ǻҶ�
	if (imgL.type() != CV_8UC1 || imgR.type() != CV_8UC1)
	{
		return -1;
	}
	//��Ȼͼ��Ĵ�С��𲻴󣬻���û��Ӱ�죬��Ϊ���ȶ��ԣ�Ҫ��ͼ��Ĵ�С��ͬ
	if (imgL.size() != imgR.size())
	{
		return -2;
	}
	refinedPairedPoints.clear();
	int row = imgL.rows;
	int col = imgR.cols;

	//��ֵ
	double a0, a1, a2;
	double b0, b1, b2;
	double h0, h1;
	//������
	double d_a0, d_a1, d_a2;
	double d_b0, d_b1, d_b2;
	double d_h0, d_h1;
	//����ֵ
	double a0_, a1_, a2_;
	double b0_, b1_, b2_;
	double h0_, h1_;

	int targetWindowK = targetWindowSize / 2;
	int searchWindowK = searchWindowSize / 2;
	//���ڴ�СͳһΪ����
	targetWindowSize = targetWindowK * 2 + 1;
	searchWindowSize = searchWindowK * 2 + 1;

	//����������Ӱ���ϵľ�������
	Point2i targetOrigin;
	Point2i searchOrigin;

	Mat targetWindow;
	Mat searchWindow;
	Mat resampleWindow = Mat(targetWindowSize, targetWindowSize, CV_8UC1);

	int x1, y1;			//g1�ϵ����꣨ȡ��������ĵ�������꣩
	double x2, y2;		//���ξ�������g2�ϵ����꣨ȡ��������ĵ�������꣩��֮��Ҫ���в�ֵ

	//����˫���Բ�ֵ
	double w11, w12, w21, w22;
	double Deltax, Deltay;
	int x11, y11;	//��ֵ�����Ͻ�����

	//���̸���
	int num;

	//g1,g2��g2����
	double g1;
	double g2;
	double g2_x;		//g2'x
	double g2_y;		//g2'y

	//���ϵ��,�Ƚϸ���ǰ�����ϵ���Ƿ����
	double R = 0;
	double R_ = 0;

	//Ƿ�������޷�ƽ��
	if (targetWindowSize <= 1)
		return -3;


	Mat matC = Mat(targetWindowSize * targetWindowSize, 8, CV_64F);
	Mat matL = Mat(targetWindowSize * targetWindowSize, 1, CV_64F);
	Mat Delta;
	Mat matCt;
	Mat matCtC;
	Mat matCtC_i;

	//�������е�ƥ���
	for (auto it = pairedPoints.begin(); it < pairedPoints.end(); it++)
	{
		//��ʼ��
		a0 = 0, a1 = 1, a2 = 0;
		b0 = 0, b1 = 0, b2 = 1;
		h0 = 0, h1 = 1;
		a0_ = 0, a1_ = 1, a2_ = 0;
		b0_ = 0, b1_ = 0, b2_ = 1;
		h0_ = 0, h1_ = 1;
		num = 0;
		R = 0, R_ = 0;

		//ȡĿ�괰�ں���������
		targetOrigin = it->first;
		searchOrigin = it->second;

		if (targetOrigin.x - targetWindowK < 0 || targetOrigin.x + targetWindowK >= col
			|| targetOrigin.y - targetWindowK < 0 || targetOrigin.y + targetWindowK >= row)
		{
			continue;
		}
		//����������ʵ��������ν����ֻ�Ǽ��ξ����ز����Ĺ�����,ֻҪ�ܱ�֤�ز���ʱ���������������ھ���
		if (searchOrigin.x - searchWindowK < 0 || searchOrigin.x + searchWindowK >= col
			|| searchOrigin.y - searchWindowK < 0 || searchOrigin.y + searchWindowK >= row)
		{
			continue;
		}
		targetWindow = imgL(Rect(Point2i(targetOrigin.x - targetWindowK, targetOrigin.y - targetWindowK), Size(targetWindowSize, targetWindowSize)));
		searchWindow = imgR(Rect(Point2i(searchOrigin.x - searchWindowK, searchOrigin.y - searchWindowK), Size(searchWindowSize, searchWindowSize)));

		//����ѭ��
		for (int i = 0; R <= R_ && i < iterationTimes; i++)
		{
			R = R_;
			num = 0;
			
			//���²���
			a0 = a0_;
			a1 = a1_;
			a2 = a2_;
			b0 = b0_;
			b1 = b1_;
			b2 = b2_;
			h0 = h0_;
			h1 = h1_;
			matC = Mat(targetWindowSize * targetWindowSize, 8, CV_64F);
			matL = Mat(targetWindowSize * targetWindowSize, 1, CV_64F);

			//���ξ����ز���
			for (int r = 0; r < targetWindowSize; r++)
			{
				y1 = r - targetWindowK;
				for (int c = 0; c < targetWindowSize; c++)
				{
					g1 = targetWindow.at<uchar>(r, c);
					x1 = c - targetWindowK;
					//����任�������
					x2 = a0 + a1 * x1 + a2 * y1;
					y2 = b0 + b1 * x1 + b2 * y1;

					//˫���Բ�ֵ
					x11 = floor(x2);
					y11 = floor(y2);
					Deltax = x2 - x11;
					Deltay = y2 - y11;

					//����,�����Ͻǵ���ݶ�Ϊ����
					if (abs(x11) >= searchWindowK || abs(y11) >= searchWindowK)
					{
						//�޷����㵼������������
						continue;
					}
					//���������е����껯Ϊ���Ͻ�ԭ��
					x11 += searchWindowK;
					y11 += searchWindowK;

					g2_x = 0.5 * (searchWindow.at<uchar>(y11, x11 + 1) - searchWindow.at<uchar>(y11, x11 - 1));
					g2_y = 0.5 * (searchWindow.at<uchar>(y11 + 1, x11) - searchWindow.at<uchar>(y11 - 1, x11));

					//˫���Բ�ֵȨ��
					w11 = (1 - Deltax) * (1 - Deltay);
					w12 = Deltax * (1 - Deltay);
					w21 = (1 - Deltax) * Deltay;
					w22 = Deltax * Deltay;
					//��ֵ
					resampleWindow.at<uchar>(r, c) = g2 = w11 * searchWindow.at<uchar>(y11, x11) + w12 * searchWindow.at<uchar>(y11, x11 + 1)
						+ w21 * searchWindow.at<uchar>(y11 + 1, x11) + w22 * searchWindow.at<uchar>(y11 + 1, x11 + 1);

					matC.at<double>(num, 0) = 1;
					matC.at<double>(num, 1) = g2;
					matC.at<double>(num, 2) = g2_x;
					matC.at<double>(num, 3) = x1 * g2_x;
					matC.at<double>(num, 4) = y1 * g2_x;
					matC.at<double>(num, 5) = g2_y;
					matC.at<double>(num, 6) = x1 * g2_y;
					matC.at<double>(num, 7) = y1 * g2_y;
					matL.at<double>(num, 0) = g1 - (h0 + h1 * g2);

					num++;
				}
			}
			if (num < 8)
			{
				break;
			}

			//���������
			matC = matC.rowRange(0, num);
			matL = matL.rowRange(0, num);

			transpose(matC, matCt);
			matCtC = matCt * matC;
			invert(matCtC, matCtC_i);
			Delta = matCtC_i * matCt * matL;

			d_h0 = Delta.at<double>(0, 0);
			d_h1 = Delta.at<double>(1, 0);
			d_a0 = Delta.at<double>(2, 0);
			d_a1 = Delta.at<double>(3, 0);
			d_a2 = Delta.at<double>(4, 0);
			d_b0 = Delta.at<double>(5, 0);
			d_b1 = Delta.at<double>(6, 0);
			d_b2 = Delta.at<double>(7, 0);

			//����У����������ϵ��
			resampleWindow = h0 + h1 * resampleWindow;
			singleWindowMatching(targetWindow, resampleWindow, R_);
			cout << R << " ";

			//�洢�²���
			a0_ = a0 + d_a0 + a0 * d_a1 + b0 * d_a2;
			a1_ = a1 + a1 * d_a1 + b1 * d_a2;
			a2_ = a2 + a2 * d_a1 + b2 * d_a2;
			b0_ = b0 + d_b0 + a0 * d_b1 + b0 * d_b2;
			b1_ = b1 + a1 * d_b1 + b1 * d_b2;
			b2_ = b2 + a2 * d_b1 + b2 * d_b2;
			h0_ = h0 + d_h0 + h0 * d_h1;
			h1_ = h1 + h1 * d_h1;
		}
		cout << endl;
		if (R < threshold)
		{
			//δ�ܴﵽ����Ҫ��
			continue;
		}
		//�������ƥ��ĵ�λ
		Point2f img1, img2;
		//����Ƭ��ƫ��
		double g1_x = 0;
		double g1_y = 0;
		double sXg1_x_2 = 0;	//sum(x*g1_x^2)
		double sYg1_y_2 = 0;	//sum(Y*g1_y^2)
		double sg1_x_2 = 0;		//sum(g1_x^2)
		double sg1_y_2 = 0;		//sum(g1_x^2)

		//���ڴ�С������
		if (targetWindowK <= 1)
		{
			img1 = Point2f(it->first);
		}
		else
		{
			for (int r = 1; r < targetWindowSize - 1; r++)
			{
				for (int c = 1; c < targetWindowSize - 1; c++)
				{
					g1_x = 0.5 * (targetWindow.at<uchar>(r, c + 1) - targetWindow.at<uchar>(r, c - 1));
					g1_y = 0.5 * (targetWindow.at<uchar>(r + 1, c) - targetWindow.at<uchar>(r - 1, c));

					sXg1_x_2 += (c - targetWindowK) * g1_x * g1_x;
					sYg1_y_2 += (r - targetWindowK) * g1_y * g1_y;
					sg1_x_2 += g1_x * g1_x;
					sg1_y_2 += g1_y * g1_y;
				}
			}
			img1.x = sXg1_x_2 / sg1_x_2;
			img1.y = sYg1_y_2 / sg1_y_2;
		}
		img2.x = a0 + a1 * img1.x + a2 * img1.y;
		img2.y = b0 + b1 * img1.x + b2 * img1.y;

		//�ɴ��ڵĹ�һ���꣬��ΪӰ������
		img1 = img1 + Point2f(targetOrigin);
		img2 = img2 + Point2f(searchOrigin);
		refinedPairedPoints.push_back(pairedPoint2f(img1, img2));
	}
	return 0;
}
