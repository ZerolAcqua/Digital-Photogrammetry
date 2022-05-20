#include <stdio.h>
#include <istream>
#include <fstream> 
#include <iostream>

#include <opencv2/core/core.hpp>   
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/imgproc/types_c.h>

#include "ImagePair.h"

#define PI 3.1415926535

using namespace std;
using namespace cv;

#define USE_FIRST 0

int main()
{
#ifdef USE_FIRST 
	//��һ��Ӱ��
	{
		//��ȡӰ��
		ImagePair oriImgPair;
		oriImgPair.mMatLImg = imread("����Ӱ��/��һ��(ģ��)/left.bmp");
		oriImgPair.mMatRImg = imread("����Ӱ��/��һ��(ģ��)/right.bmp");
		//��ȡ�ڷ�λԪ��
		{
			ifstream infile;   //������

			double jo, io;
			double a11, a12, a21, a22;
			double b11, b12, b21, b22;
			double f;
			double pixelSize;

			infile.open("����Ӱ��/��һ��(ģ��)/left.iop", ios::in);
			if (!infile.is_open())
			{
				cout << "position file not found." << endl;
			}
			else {
				char temp[128];

				infile >> jo >> io;
				infile.getline(temp, 128);		//��������
				infile >> a11 >> a12;
				infile >> a21 >> a22;
				infile.getline(temp, 128);		//��������
				infile >> b11 >> b12;
				infile >> b21 >> b22;
				infile.getline(temp, 128);		//��������
				infile >> f >> pixelSize;

				infile.close();   //�ر��ļ�
			}

			oriImgPair.mMatPix2ImgRot = { {	b11 * pixelSize,	b12 * pixelSize	},
											{	b21 * pixelSize,	b22 * pixelSize	} };
			oriImgPair.mMatPix2ImgBar = Matrix({ (-b11 * jo - b12 * io) * pixelSize	,
													(-b21 * jo - b22 * io) * pixelSize }).transpose();

			oriImgPair.mMatImg2PixRot = { {	a11 / pixelSize,	a12 / pixelSize	},
											{	a21 / pixelSize,	a22 / pixelSize	} };
			oriImgPair.mMatImg2PixBar = Matrix({ jo,io }).transpose();

			oriImgPair.mf = f;
		}
		//��ȡ�ⷽλԪ��
		{
			ifstream infile;   //������

			double Xs, Ys, Zs, phi, omega, kappa;

			infile.open("����Ӱ��/��һ��(ģ��)/left_right_�ⷽλԪ��.txt", ios::in);
			if (!infile.is_open())
			{
				cout << "position file not found." << endl;
			}
			else {
				char temp[128];
				infile.getline(temp, 128);		//����ע��

				infile >> Xs >> Ys >> Zs >> phi >> omega >> kappa;
				double L[6] = { Xs, Ys, Zs, phi / 180 * PI, omega / 180 * PI, kappa / 180 * PI };
				oriImgPair.setLEOP(L, 6);

				infile >> Xs >> Ys >> Zs >> phi >> omega >> kappa;
				double R[6] = { Xs, Ys, Zs, phi / 180 * PI, omega / 180 * PI, kappa / 180 * PI };
				oriImgPair.setREOP(R, 6);

				infile.close();   //�ر��ļ�
			}
		}

		ImagePair newImgPair = oriImgPair.getEpipolarImgPair();


		imwrite("����Ӱ��/left1.bmp", newImgPair.mMatLImg);
		imwrite("����Ӱ��/right1.bmp", newImgPair.mMatRImg);


		//imshow("��Ӱ��", oriImgPair.mMatLImg);
		//imshow("��Ӱ��", oriImgPair.mMatLImg);
		//waitKey();

		//ImagePair newImgPair = oriImgPair.getEpipolarImgPair();
		//oriImgPair.~ImagePair();

		//imshow("��Ӱ��", newImgPair.mMatLImg);
		//imshow("��Ӱ��", newImgPair.mMatLImg);
		//waitKey();

		//cout << newImgPair.mMatImg2PixRot;
	}
#endif

#ifndef USE_FIRST 
	//�ڶ���Ӱ��
	{
		//��ȡӰ��
		ImagePair oriImgPair;
		oriImgPair.mMatLImg = imread("����Ӱ��/�ڶ���(���˻����)/Images/IMG_0182.jpg");
		oriImgPair.mMatRImg = imread("����Ӱ��/�ڶ���(���˻����)/Images/IMG_0183.jpg");
		//��ȡ�ڷ�λԪ��
		{
			ifstream infile;   //������

			double x0, y0;

			double f;
			double pixelSize;
			int w, h;

			infile.open("����Ӱ��/�ڶ���(���˻����)/camera.txt", ios::in);
			if (!infile.is_open())
			{
				cout << "position file not found." << endl;
			}
			else {
				char temp[128];

				infile.getline(temp, 128);		//����ע��
				infile >> temp >> f;
				infile >> temp >> x0 >> y0;
				infile >> temp >> pixelSize;
				infile >> temp >> w;
				infile >> temp >> h;

				infile.close();   //�ر��ļ�
			}

			oriImgPair.mMatPix2ImgRot = {	{		pixelSize,		0			},
											{		0,				pixelSize	} };
			oriImgPair.mMatPix2ImgBar = Matrix({	-w / 2.0 * pixelSize - x0	,
													-h / 2.0 * pixelSize - y0 }).transpose();

			oriImgPair.mMatImg2PixRot = {	{		1/pixelSize,	0			},
											{		0,				1/pixelSize	} };
			oriImgPair.mMatImg2PixBar = Matrix({	x0 / pixelSize + w / 2.0	,
													y0 / pixelSize + h / 2 }).transpose();

			oriImgPair.mf = f;
		}
		//��ȡ�ⷽλԪ��
		{
			ifstream infile;   //������

			double Xs, Ys, Zs, phi, omega, kappa;

			infile.open("����Ӱ��/�ڶ���(���˻����)/IMG_0182.jpg.aop", ios::in);
			if (!infile.is_open())
			{
				cout << "position file not found." << endl;
			}
			else {
				char temp[128];
				infile.getline(temp, 128);		//��������
				infile >> Xs >> Ys >> Zs;
				infile >> phi >> omega >> kappa;
				double L[6] = { Xs, Ys, Zs, phi, omega , kappa };
				oriImgPair.setLEOP(L, 6);
				infile.close();   //�ر��ļ�
			}

			infile.open("����Ӱ��/�ڶ���(���˻����)/IMG_0183.jpg.aop", ios::in);
			if (!infile.is_open())
			{
				cout << "position file not found." << endl;
			}
			else {
				char temp[128];
				infile.getline(temp, 128);		//��������
				infile >> Xs >> Ys >> Zs;
				infile >> phi >> omega >> kappa;
				double R[6] = { Xs, Ys, Zs, phi, omega , kappa };
				oriImgPair.setREOP(R, 6);
				infile.close();   //�ر��ļ�
			}

		}

		ImagePair newImgPair = oriImgPair.getEpipolarImgPair();


		imwrite("����Ӱ��/left2.bmp", newImgPair.mMatLImg);
		imwrite("����Ӱ��/right2.bmp", newImgPair.mMatRImg);
	}
#endif

}