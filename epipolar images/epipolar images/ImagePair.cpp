#include "ImagePair.h"

bool ImagePair::setLEOP(double* lEOP, int size)
{
	if (size < 6)
	{
		return false;
	}
	else
	{
		memcpy(this->mLEOP, lEOP, 6 * sizeof(double));
		this->mLHaveEularAngle = true;
		double phi = this->mLEOP[3], omega = this->mLEOP[4], kappa = this->mLEOP[5];
		Matrix R_phi = {	{	cos(phi),		0,				-sin(phi)	},
							{	0,				1,				0			},
							{	sin(phi),		0,				cos(phi)	} };

		Matrix R_omega = {	{	1,				0,				0			},
							{	0,				cos(omega),		-sin(omega)	},
							{	0,				sin(omega),		cos(omega)	} };

		Matrix R_kappa = {	{	cos(kappa),		-sin(kappa),	0			},
							{	sin(kappa),		cos(kappa),		0			},
							{	0,				0,				1			} };

		this->mMatLeftR = R_phi * R_omega * R_kappa;
		return true;
	}
}
bool ImagePair::setREOP(double* rEOP, int size)
{
	if (size < 6)
	{
		return false;
	}
	else
	{
		memcpy(this->mREOP, rEOP, 6 * sizeof(double));
		this->mRHaveEularAngle = true;
		double phi = this->mREOP[3], omega = this->mREOP[4], kappa = this->mREOP[5];
		Matrix R_phi = {	{	cos(phi),		0,				-sin(phi)	},
							{	0,				1,				0			},
							{	sin(phi),		0,				cos(phi)	} };

		Matrix R_omega = {	{	1,				0,				0			},
							{	0,				cos(omega),		-sin(omega)	},
							{	0,				sin(omega),		cos(omega)	} };

		Matrix R_kappa = {	{	cos(kappa),		-sin(kappa),	0			},
							{	sin(kappa),		cos(kappa),		0			},
							{	0,				0,				1			} };

		this->mMatRightR = R_phi * R_omega * R_kappa;
		return true;
	}
}

ImagePair ImagePair::getEpipolarImgPair()
{
	// 若已经是核线影像了，则不进行处理
	if (this->mIsEpipolarImg == true)
	{
		return *this;
	}

	//------------------ 计算核线影像相应的外方位元素 ------------------
	ImagePair epipolar = calculateEpipolarOP();

	//------------------ 核线影像重采样 ------------------
	double f = this->mf;
	double fn = epipolar.mf;
	int epiRows = epipolar.mMatLImg.rows;
	int epiCols = epipolar.mMatLImg.cols;
	Matrix epiPixCoor = Matrix::zeros(2, 1);
	Matrix oriPixCoor = Matrix::zeros(2, 1);
	//----- 左片 -----
	Matrix ML = epipolar.mMatLeftR.transpose() * this->mMatLeftR;

	for (int r = 0; r < epiRows; r++)
	{
		for (int c = 0; c < epiCols; c++)
		{
			//计算坐标

			//epiPixCoor[0][0] = c;
			//epiPixCoor[1][0] = epiCols - 1 - r;
			//cout << epiPixCoor << endl;
			//epiPixCoor = epipolar.pix2Img(epiPixCoor);
			//cout << epiPixCoor << endl;
			//this->epipolarImg2OriImg(epiPixCoor, f, fn, ML);
			//cout << epiPixCoor << endl;
			//oriPixCoor = this->img2Pix(epiPixCoor);
			//cout << oriPixCoor << endl;
			//epipolar.mMatLImg.at<Vec3b>(r, c) = interpPix(this->mMatLImg, oriPixCoor);

			epiPixCoor[0][0] = c;
			epiPixCoor[1][0] = epiRows - 1 - r;
			oriPixCoor = this->img2Pix(this->epipolarImg2OriImg(epipolar.pix2Img(epiPixCoor), f, fn, ML));
			epipolar.mMatLImg.at<Vec3b>(r, c) = interpPix(this->mMatLImg, oriPixCoor);
		}
	}

	//----- 右片 -----
	Matrix MR = epipolar.mMatRightR.transpose() * this->mMatRightR;

	for (int r = 0; r < epiRows; r++)
	{
		for (int c = 0; c < epiCols; c++)
		{
			//计算坐标
			epiPixCoor[0][0] = c;
			epiPixCoor[1][0] = epiRows - 1 - r;
			oriPixCoor = this->img2Pix(this->epipolarImg2OriImg(epipolar.pix2Img(epiPixCoor), f, fn, MR));
			epipolar.mMatRImg.at<Vec3b>(r, c) = interpPix(this->mMatRImg, oriPixCoor);
		}
	}


	return epipolar;
}

ImagePair ImagePair::calculateEpipolarOP()
{
	ImagePair epipolar;
	epipolar.mIsEpipolarImg = true;

	/*********************************************************************
	* 整个“水平”影像的生成过程中，使用的都是phi-kappa-omega转角系统，
	* 以方便地将影像行旋转到与基线平行的“水平”平面上。因此，这个过程中
	* 还涉及利用旋转矩阵计算欧拉角的过程，以获得最适合的omega角。
	*
	* “归一化平面”，就是指利用“水平”影像生成核线的过程中，指定的与
	*  基线平行的平面，通过计算这个归一化平面上的像片的外方位元素，就可
	* 以建立原始影像与“水平”影像之间的关系，从而进行重采样过程。
	*
	* 对于已换算到以像主点为原点的像点坐标，归一化平面影像与原始影像像点
	* 坐标存在以下关系：
	*
	*
	*	-			-										-		-
	*	|	x_n		|										|	x	|
	*	|	y_n		|	=	lambda/lambda_n * R_n^T * R*	|	y	|
	*	|	-f_n	|										|	-f	|
	*	-			-										-		-
	*
	*********************************************************************/

	//------------------ 归一化平面参数计算 ------------------
	double pixSize =
		sqrt((this->mMatPix2ImgRot * this->mMatPix2ImgRot.transpose())[0][0]);	//像元大小

	double Bx = mREOP[0] - mLEOP[0];
	double By = mREOP[1] - mLEOP[1];
	double Bz = mREOP[2] - mLEOP[2];
	double omega_n = 0; // TODO:计算最佳omega角
	double phi_n = atan(Bz / Bx);
	double kappa_n = atan(By / sqrt(Bx * Bx + Bz * Bz));
	Matrix R_phi = {	{	cos(phi_n),		0,				-sin(phi_n)		},
						{	0,				1,				0				},
						{	sin(phi_n),		0,				cos(phi_n)		} };

	Matrix R_kappa = { {	cos(kappa_n),	-sin(kappa_n),	0				},
						{	sin(kappa_n),	cos(kappa_n),	0				},
						{	0,				0,				1				} };

	Matrix R_omega = { {	1,				0,				0				},
						{	0,				cos(omega_n),	-sin(omega_n)	},
						{	0,				sin(omega_n),	cos(omega_n)	} };

	//内方位元素
	epipolar.mMatPix2ImgRot = Matrix::eye(2) * pixSize;
	epipolar.mMatImg2PixRot = Matrix::eye(2) / pixSize;
	epipolar.mf = this->mf;
	//外方位线元素
	epipolar.mLEOP[0] = this->mLEOP[0];
	epipolar.mLEOP[1] = this->mLEOP[1];
	epipolar.mLEOP[2] = this->mLEOP[2];
	epipolar.mREOP[0] = this->mREOP[0];
	epipolar.mREOP[1] = this->mREOP[1];
	epipolar.mREOP[2] = this->mREOP[2];
	//TODO：外方位角元素欧拉角
	epipolar.mLHaveEularAngle = epipolar.mRHaveEularAngle = false;
	epipolar.mLEOP[3] = epipolar.mLEOP[4] = epipolar.mLEOP[5] = 0;
	epipolar.mREOP[3] = epipolar.mREOP[4] = epipolar.mREOP[5] = 0;
	epipolar.mMatLeftR = R_phi * R_kappa * R_omega;
	epipolar.mMatRightR = epipolar.mMatLeftR;


	//------------------ 核线影像的尺寸计算 ------------------
	//原始影像下左、上左、上右、下右角点坐标
	Matrix cornerDL = Matrix::zeros(2, 1), cornerUL = Matrix::zeros(2, 1)
		, cornerUR = Matrix::zeros(2, 1), cornerDR = Matrix::zeros(2, 1);

	//----- 左片 -----
	//计算左片关系式 M = R_n^T * R
	Matrix ML = epipolar.mMatLeftR.transpose() * this->mMatLeftR;
	int rowL = this->mMatLImg.rows;
	int colL = this->mMatLImg.cols;

	//左原始影像下左、上左、上右、下右角点像素坐标
	cornerDL[0][0] = 0, cornerDL[1][0] = 0;
	cornerUL[0][0] = 0, cornerUL[1][0] = rowL - 1;
	cornerUR[0][0] = colL - 1, cornerUR[1][0] = rowL - 1;
	cornerDR[0][0] = colL - 1, cornerDR[1][0] = 0;

	//左原始影像下左、上左、上右、下右角点影像坐标
	cornerDL = this->pix2Img(cornerDL);
	cornerUL = this->pix2Img(cornerUL);
	cornerUR = this->pix2Img(cornerUR);
	cornerDR = this->pix2Img(cornerDR);

	//左核线影像下左、上左、上右、下右角点影像坐标
	cornerDL = oriImg2EpipolarImg(cornerDL, this->mf, epipolar.mf, ML);
	cornerUL = oriImg2EpipolarImg(cornerUL, this->mf, epipolar.mf, ML);
	cornerUR = oriImg2EpipolarImg(cornerUR, this->mf, epipolar.mf, ML);
	cornerDR = oriImg2EpipolarImg(cornerDR, this->mf, epipolar.mf, ML);

	//左核线影像下左、上左、上右、下右角点像素坐标
	cornerDL = epipolar.img2Pix(cornerDL);
	cornerUL = epipolar.img2Pix(cornerUL);
	cornerUR = epipolar.img2Pix(cornerUR);
	cornerDR = epipolar.img2Pix(cornerDR);

	//左核线影像的下左、上右角点像素坐标
	Matrix leftEpipolarPixDL, leftEpipolarPixUR;
	leftEpipolarPixDL = findMin(cornerDL, cornerUL, cornerUR, cornerDR);
	leftEpipolarPixUR = findMax(cornerDL, cornerUL, cornerUR, cornerDR);

	//----- 右片 -----
	//计算右片关系式 M = R_n^T * R
	Matrix MR = epipolar.mMatRightR.transpose() * this->mMatRightR;
	int rowR = this->mMatRImg.rows;
	int colR = this->mMatRImg.cols;


	//右原始影像下左、上左、上右、下右角点像素坐标
	cornerDL[0][0] = 0, cornerDL[1][0] = 0;
	cornerUL[0][0] = 0, cornerUL[1][0] = rowR - 1;
	cornerUR[0][0] = colR - 1, cornerUR[1][0] = rowR - 1;
	cornerDR[0][0] = colR - 1, cornerDR[1][0] = 0;

	//右原始影像下左、上左、上右、下右角点影像坐标
	cornerDL = this->pix2Img(cornerDL);
	cornerUL = this->pix2Img(cornerUL);
	cornerUR = this->pix2Img(cornerUR);
	cornerDR = this->pix2Img(cornerDR);

	//右核线影像下左、上左、上右、下右角点影像坐标
	cornerDL = oriImg2EpipolarImg(cornerDL, this->mf, epipolar.mf, MR);
	cornerUL = oriImg2EpipolarImg(cornerUL, this->mf, epipolar.mf, MR);
	cornerUR = oriImg2EpipolarImg(cornerUR, this->mf, epipolar.mf, MR);
	cornerDR = oriImg2EpipolarImg(cornerDR, this->mf, epipolar.mf, MR);

	//右核线影像下左、上左、上右、下右角点像素坐标
	cornerDL = epipolar.img2Pix(cornerDL);
	cornerUL = epipolar.img2Pix(cornerUL);
	cornerUR = epipolar.img2Pix(cornerUR);
	cornerDR = epipolar.img2Pix(cornerDR);

	//右核线影像的下左、上右角点像素坐标
	Matrix rightEpipolarPixDL, rightEpipolarPixUR;
	rightEpipolarPixDL = findMin(cornerDL, cornerUL, cornerUR, cornerDR);
	rightEpipolarPixUR = findMax(cornerDL, cornerUL, cornerUR, cornerDR);


	//----- 核线影像最终范围与内方位元素 -----
	//核线影像的下左、上右角点像素坐标
	Matrix epipolarPixDL, epipolarPixUR;
	epipolarPixDL = findMin(leftEpipolarPixDL, leftEpipolarPixUR, rightEpipolarPixDL, rightEpipolarPixUR);
	epipolarPixUR = findMax(leftEpipolarPixDL, leftEpipolarPixUR, rightEpipolarPixDL, rightEpipolarPixUR);
	//影像大小
	int rows = int(epipolarPixUR[1][0] - epipolarPixDL[1][0]) + 1;
	int cols = int(epipolarPixUR[0][0] - epipolarPixDL[0][0]) + 1;
	epipolar.mMatLImg = Mat(rows, cols, CV_8UC3, Scalar(0, 0, 0));
	epipolar.mMatRImg = Mat(rows, cols, CV_8UC3, Scalar(0, 0, 0));

	//调整内方位元素（像主点）
	double jo = -int(epipolarPixDL[0][0]);
	double io = -int(epipolarPixDL[1][0]);
	epipolar.mMatPix2ImgBar = Matrix({ -jo * pixSize,-io * pixSize }).transpose();
	epipolar.mMatImg2PixBar = Matrix({ jo,io }).transpose();

	return epipolar;
}




Matrix ImagePair::pix2Img(const Matrix& pixCoor)
{	
	return 	mMatPix2ImgRot * pixCoor + mMatPix2ImgBar;
}
Matrix ImagePair::img2Pix(const Matrix& imgCoor)
{
	return 	mMatImg2PixRot * imgCoor + mMatImg2PixBar;
}
Matrix ImagePair::oriImg2EpipolarImg(const Matrix& oriImgCoor, double f, double fn, const Matrix& M)
{
	Matrix temp = Matrix::zeros(2, 1);
	double x = oriImgCoor[0][0], y = oriImgCoor[1][0];
	temp[0][0] = -fn *	(M[0][0] * x + M[0][1] * y - M[0][2] * f) /
						(M[2][0] * x + M[2][1] * y - M[2][2] * f);
	temp[1][0] = -fn *	(M[1][0] * x + M[1][1] * y - M[1][2] * f) /
						(M[2][0] * x + M[2][1] * y - M[2][2] * f);
	return temp;
}
Matrix ImagePair::epipolarImg2OriImg(const Matrix& epipolarImgCoor, double f, double fn, const Matrix& M)
{
	Matrix temp = Matrix::zeros(2, 1);
	double x = epipolarImgCoor[0][0], y = epipolarImgCoor[1][0];
	temp[0][0] = -f *	(M[0][0] * x + M[1][0] * y - M[2][0] * fn) /
						(M[0][2] * x + M[1][2] * y - M[2][2] * fn);
	temp[1][0] = -f *	(M[0][1] * x + M[1][1] * y - M[2][1] * fn) /
						(M[0][2] * x + M[1][2] * y - M[2][2] * fn);
	return temp;
}











Matrix ImagePair::findMin(const Matrix& coor1, const Matrix& coor2, const Matrix& coor3, const Matrix& coor4)
{
	Matrix temp = Matrix::zeros(2, 1);
	temp[0][0] = MIN(coor1[0][0], MIN(coor2[0][0], MIN(coor3[0][0], coor4[0][0])));
	temp[1][0] = MIN(coor1[1][0], MIN(coor2[1][0], MIN(coor3[1][0], coor4[1][0])));
	return temp;
}

Matrix ImagePair::findMax(const Matrix& coor1, const Matrix& coor2, const Matrix& coor3, const Matrix& coor4)
{
	Matrix temp = Matrix::zeros(2, 1);
	temp[0][0] = MAX(coor1[0][0], MAX(coor2[0][0], MAX(coor3[0][0], coor4[0][0])));
	temp[1][0] = MAX(coor1[1][0], MAX(coor2[1][0], MAX(coor3[1][0], coor4[1][0])));
	return temp;
}

Vec3b ImagePair::interpPix(const Mat& img,Matrix pixCoor)
{
	int cols = img.cols;
	int rows = img.rows;

	//opencv图像坐标系的原点在左上角
	double x = pixCoor[0][0];
	double y = rows - 1 - pixCoor[1][0];
	int ix = floor(x);
	int iy = floor(y);
	double dx = ix - floor(x);
	double dy = iy - floor(y);

	if (x < 0 || x >= cols - 1
		|| y < 0 || y >= rows - 1)
	{
		return Vec3b(0, 0, 0);
	}
	else
	{
		return (1 - dx) * (1 - dy) * Vec3f(img.at<Vec3b>(iy, ix))
			+ dx * (1 - dy) * Vec3f(img.at<Vec3b>(iy, ix + 1))
			+ (1 - dx) * dy * Vec3f(img.at<Vec3b>(iy + 1, ix))
			+ dx * dy * Vec3f(img.at<Vec3b>(iy + 1, ix + 1));
	}
}