#include <opencv2/core/core.hpp>   
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/imgproc/types_c.h>
#include "Matrix.h"
using namespace std;
using namespace cv;

#pragma once
class ImagePair
{
public:
	//����Ӱ��
	Mat mMatLImg;
	Mat mMatRImg;
	//�ⷽλԪ��EOP[6]: Xs Ys Zs phi omega kappa
	double mLEOP[6] = { 0 };
	double mREOP[6] = { 0 };
	bool mLHaveEularAngle = false;	// ŷ������Ϣ�Ƿ���Ч����Ϊ�е�Ӱ�����ת��û�м���
	bool mRHaveEularAngle = false;	// ŷ������Ϣ�Ƿ���Ч����Ϊ�е�Ӱ�����ת��û�м���
	Matrix mMatLeftR = Matrix::eye(3);		//��ת����
	Matrix mMatRightR = Matrix::eye(3);

	//�ڷ�λԪ�أ�**��Ϊ����Ӱ���ڷ�λԪ����ͬ**
	// f
	// x0,y0���ڶ�����󣬱�Ϊ0�����ؿ���
	//����������,ֱ��תΪ��������Ϊԭ���Ӱ�����ꡣ
	double mf = 1;
	//�Ƿ�Ϊ����Ӱ��
	bool mIsEpipolarImg = false;
	//�ڶ������,ֱ�����������������������Ϊԭ���Ӱ�������ת��
	//ImgCoor = mMatPix2ImgRot * PixCoor + mMatPix2ImgBar
	Matrix mMatPix2ImgRot = Matrix::eye(2);
	Matrix mMatPix2ImgBar = Matrix::zeros(2,1);
	//PixCoor = mMatImg2PixRot * ImgCoor + mMatImg2PixBar
	Matrix mMatImg2PixRot = Matrix::eye(2);
	Matrix mMatImg2PixBar = Matrix::zeros(2, 1);


public:
	/**
	*	@brief ���к���
	*	������Ӱ���ⷽλԪ��
	*
	*	@param lEOP		��Ӱ���ⷽλԪ��
	*	@param size		���鳤��
	* 
	*	@return			true�������óɹ�
	*					false��������ʧ��
	*/
	bool setLEOP(double* lEOP,int size);
	/**
	*	@brief ���к���
	*	������Ӱ���ⷽλԪ��
	*
	*	@param rEOP		��Ӱ���ⷽλԪ��
	*	@param size		���鳤��
	*
	*	@return			true�������óɹ�
	*					false��������ʧ��
	*/
	bool setREOP(double* rEOP, int size);
	/**
	*	@brief ���к���
	*   ��ȡ����Ӱ���
	*
	*	@return ImagePair	����Ӱ���
	*/
	ImagePair getEpipolarImgPair();
	/**
	*	@brief ���к���
	*   ����������ת��ΪӰ������
	*
	*	@param pixCoor	�����������
	*	@return Matrix	Ӱ���������
	*/
	Matrix pix2Img(const Matrix& pixCoor);
	/**
	*	@brief ���к���
	*   ��Ӱ������ת��Ϊ��������
	*
	*	@param imgCoor	Ӱ���������
	*	@param Matrix	�����������
	*/
	Matrix img2Pix(const Matrix& imgCoor);

	/**
	*	@brief ���к���
	*   ��ԭʼӰ������ת��Ϊ����Ӱ������
	*
	*	@param oriImgCoor	ԭʼӰ������
	* 	@param f			ԭʼӰ������
	* 	@param fn			����Ӱ������
	* 	@param M			��ת����
	*	
	*	@return Matrix		����Ӱ������
	*/
	Matrix oriImg2EpipolarImg(const Matrix& oriImgCoor,double f,double fn,const Matrix& M);

	/**
	*	@brief ���к���
	*   ������Ӱ������ת��ΪԭʼӰ������
	*
	*	@param epipolarImgCoor	����Ӱ������
	* 	@param f				ԭʼӰ������
	* 	@param fn				����Ӱ������
	* 	@param M				��ת����
	*	@return	Matrix			ԭʼӰ������
	*/
	Matrix epipolarImg2OriImg(const Matrix& epipolarImgCoor, double f, double fn, const Matrix& M);

private:
	/**
	*	@brief ˽�к���
	*   ���غ���Ӱ���
	*
	*	@return ImagePair	���غ���Ӱ���
	*/
	ImagePair calculateEpipolarOP();
	
	/**
	*	@brief ˽�к���
	*   ����Ӱ��ָ���������������ֵ������ֵ����
	* 
	*	@param img		����ֵ��Ӱ��
	*	@param pixCoor	ָ����������
	*
	*	@return Vec3b	��������ֵ����
	*/
	Vec3b interpPix(const Mat& img, Matrix pixCoor);

	/**
	*	@brief ˽�к���
	*   �����ĸ�ƽ���������Сx�����y���꣨�Ծ�����ʽ��
	*
	*	@param coor1	ƽ���������
	*	@param coor2	ƽ���������
	*	@param coor3	ƽ���������
	*	@param coor4	ƽ���������
	*
	*	@return Matrix	��С��x������С��y��ϵ��������
	*/
	static Matrix findMin(const Matrix& coor1, const Matrix& coor2, const Matrix& coor3, const Matrix& coor4);

	/**
	*	@brief ˽�к���
	*   �����ĸ�ƽ����������x�����y���꣨�Ծ�����ʽ��
	*
	*	@param coor1	ƽ���������
	*	@param coor2	ƽ���������
	*	@param coor3	ƽ���������
	*	@param coor4	ƽ���������
	*
	*	@return Matrix	����x��������y��ϵ��������
	*/
	static Matrix findMax(const Matrix& coor1, const Matrix& coor2, const Matrix& coor3, const Matrix& coor4);
};

