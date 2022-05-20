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
	//左右影像
	Mat mMatLImg;
	Mat mMatRImg;
	//外方位元素EOP[6]: Xs Ys Zs phi omega kappa
	double mLEOP[6] = { 0 };
	double mREOP[6] = { 0 };
	bool mLHaveEularAngle = false;	// 欧拉角信息是否有效，因为有的影像对旋转角没有计算
	bool mRHaveEularAngle = false;	// 欧拉角信息是否有效，因为有的影像对旋转角没有计算
	Matrix mMatLeftR = Matrix::eye(3);		//旋转矩阵
	Matrix mMatRightR = Matrix::eye(3);

	//内方位元素，**认为左右影像内方位元素相同**
	// f
	// x0,y0在内定向处理后，变为0，不必考虑
	//即像素坐标,直接转为以像主点为原点的影像坐标。
	double mf = 1;
	//是否为核线影像
	bool mIsEpipolarImg = false;
	//内定向矩阵,直接完成像素坐标与以像主点为原点的影像坐标的转换
	//ImgCoor = mMatPix2ImgRot * PixCoor + mMatPix2ImgBar
	Matrix mMatPix2ImgRot = Matrix::eye(2);
	Matrix mMatPix2ImgBar = Matrix::zeros(2,1);
	//PixCoor = mMatImg2PixRot * ImgCoor + mMatImg2PixBar
	Matrix mMatImg2PixRot = Matrix::eye(2);
	Matrix mMatImg2PixBar = Matrix::zeros(2, 1);


public:
	/**
	*	@brief 公有函数
	*	设置左影像外方位元素
	*
	*	@param lEOP		左影像外方位元素
	*	@param size		数组长度
	* 
	*	@return			true——设置成功
	*					false——设置失败
	*/
	bool setLEOP(double* lEOP,int size);
	/**
	*	@brief 公有函数
	*	设置右影像外方位元素
	*
	*	@param rEOP		右影像外方位元素
	*	@param size		数组长度
	*
	*	@return			true——设置成功
	*					false——设置失败
	*/
	bool setREOP(double* rEOP, int size);
	/**
	*	@brief 公有函数
	*   获取核线影像对
	*
	*	@return ImagePair	核线影像对
	*/
	ImagePair getEpipolarImgPair();
	/**
	*	@brief 公有函数
	*   将像素坐标转化为影像坐标
	*
	*	@param pixCoor	像素坐标矩阵
	*	@return Matrix	影像坐标矩阵
	*/
	Matrix pix2Img(const Matrix& pixCoor);
	/**
	*	@brief 公有函数
	*   将影像坐标转化为像素坐标
	*
	*	@param imgCoor	影像坐标矩阵
	*	@param Matrix	像素坐标矩阵
	*/
	Matrix img2Pix(const Matrix& imgCoor);

	/**
	*	@brief 公有函数
	*   将原始影像坐标转化为核线影像坐标
	*
	*	@param oriImgCoor	原始影像坐标
	* 	@param f			原始影像主距
	* 	@param fn			核线影像主距
	* 	@param M			旋转矩阵
	*	
	*	@return Matrix		核线影像坐标
	*/
	Matrix oriImg2EpipolarImg(const Matrix& oriImgCoor,double f,double fn,const Matrix& M);

	/**
	*	@brief 公有函数
	*   将核线影像坐标转化为原始影像坐标
	*
	*	@param epipolarImgCoor	核线影像坐标
	* 	@param f				原始影像主距
	* 	@param fn				核线影像主距
	* 	@param M				旋转矩阵
	*	@return	Matrix			原始影像坐标
	*/
	Matrix epipolarImg2OriImg(const Matrix& epipolarImgCoor, double f, double fn, const Matrix& M);

private:
	/**
	*	@brief 私有函数
	*   返回核线影像对
	*
	*	@return ImagePair	返回核线影像对
	*/
	ImagePair calculateEpipolarOP();
	
	/**
	*	@brief 私有函数
	*   返回影像指定像素坐标的像素值，含插值过程
	* 
	*	@param img		待插值的影像
	*	@param pixCoor	指定像素坐标
	*
	*	@return Vec3b	返回像素值向量
	*/
	Vec3b interpPix(const Mat& img, Matrix pixCoor);

	/**
	*	@brief 私有函数
	*   返回四个平面坐标的最小x坐标和y坐标（以矩阵形式）
	*
	*	@param coor1	平面坐标矩阵
	*	@param coor2	平面坐标矩阵
	*	@param coor3	平面坐标矩阵
	*	@param coor4	平面坐标矩阵
	*
	*	@return Matrix	最小的x，和最小的y组合的坐标矩阵
	*/
	static Matrix findMin(const Matrix& coor1, const Matrix& coor2, const Matrix& coor3, const Matrix& coor4);

	/**
	*	@brief 私有函数
	*   返回四个平面坐标的最大x坐标和y坐标（以矩阵形式）
	*
	*	@param coor1	平面坐标矩阵
	*	@param coor2	平面坐标矩阵
	*	@param coor3	平面坐标矩阵
	*	@param coor4	平面坐标矩阵
	*
	*	@return Matrix	最大的x，和最大的y组合的坐标矩阵
	*/
	static Matrix findMax(const Matrix& coor1, const Matrix& coor2, const Matrix& coor3, const Matrix& coor4);
};

