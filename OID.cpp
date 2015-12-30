#include <iostream>
#include <string>
#include <opencv2\core\core.hpp>
#include <opencv2\highgui\highgui.hpp>

using namespace std;
using namespace cv;

//Kasa() function prototype
Mat Kasa (Mat XY);

int main ( void ) {

	//Read images
	char *imgName1 = "per.tif";
	char *imgName2 = "img16.tif";

	Mat img1 = imread(imgName1, CV_LOAD_IMAGE_GRAYSCALE);
	Mat img2 = imread(imgName2, CV_LOAD_IMAGE_COLOR);

	if ( !img1.data ) {

		cout << "Could not open or find " << imgName1 << endl;
		waitKey(0);
		return -1;
	}

	if ( !img2.data ) {

		cout << "Could not open of find " << imgName2 << endl;
		waitKey(0);
		return -1;
	}

	namedWindow("Original", CV_WINDOW_AUTOSIZE);
	imshow("Original", img2);

	MatIterator_<uchar> it1, end1;

	int size = 0;

	//Calculate how many points there are in the outline
	for (it1 = img1.begin<uchar>(), end1 = img1.end<uchar>(); it1 != end1; ++it1) {

		if ((*it1) == 255) {
			++size;
		}
	}

	//Initialize the XY array which stores the (x,y) coordinates of the outline points.
	//Will be used as input for Kasa. Must be of float type.
	Mat XY = Mat(size, 2, CV_32F);

	int q = 0;

	//Store outline points in XY
	for (int i=0; i<img1.rows; ++i) {
		for (int j=0; j<img1.cols; ++j) {

			if (img1.at<uchar>(i,j) == 255) {
				XY.at<float>(q,0) = (float)i;
				XY.at<float>(q,1) = (float)j;
				++q;
			}
		}
	}

	//Initialize Fit vector which stores the (x,y) coordinates and radius of the circle.
	Mat Fit = Mat(1, 3, CV_32F);

	//Perform circle fitting with Kasa, the function will return a 1x3 vector. 
	Fit = Kasa (XY);
	cout << "Fitting circle: (x, y, R) = " << Fit << endl;

	//img2 will now have the fitted outline.
	for (int i = 0; i<img1.rows; ++i) {
		for (int j = 0; j<img1.cols; ++j) {
			if (std::sqrt(std::pow(i-Fit.at<float>(0),2) + std::pow(j-Fit.at<float>(1),2)) <= Fit.at<float>(2))
				img1.at<uchar>(i,j) = 255;
			else
				img1.at<uchar>(i,j) = 0;
		}
	}
			

	MatIterator_<Vec3b> it2, end2;

	//Paste outline on original image
	for (it1 = img1.begin<uchar>(), it2 = img2.begin<Vec3b>(), 
		end1 = img1.end<uchar>(), end2 = img2.end<Vec3b>(); 
		it1 != end1, it2 != end2; 
		++it1, ++it2) {

			if ((*it1) != 255) {
				(*it2)[0] = 0;
				(*it2)[1] = 0;
				(*it2)[2] = 0;
			}
	}

	namedWindow("Final", CV_WINDOW_AUTOSIZE);
	imshow("Final", img2);
	waitKey(0);

	return 0;

}

//Kasa
Mat Kasa (Mat XY) {
	
	//Matrices must be of type float, since we'll do operations between them.
	Mat A = Mat(XY.rows, 1, CV_32F); 
	Mat B = Mat(XY.rows, 1, CV_32F); 

	for (int i = 0; i < XY.rows; ++i) {

		A.at<float>(i) = XY.at<float>(i,0);
		B.at<float>(i) = XY.at<float>(i,1);
	}

	Mat C = Mat(XY.rows, 1, CV_32F);
	C = A.mul(A) + B.mul(B);

	Mat I = Mat::ones(XY.rows, 1, CV_32F);

	//XY1 is the horizontal concatenation of XY, I
	Mat XY1 = Mat(XY.rows, XY.cols + I.cols, CV_32F);

	hconcat(XY, I, XY1);
	
	Mat Fit;
	solve(XY1, C, Fit, DECOMP_SVD);

	//Temp variables.
	float Fit1 = Fit.at<float>(0);
	float Fit2 = Fit.at<float>(1);
	float Fit3 = Fit.at<float>(2);
 
	Fit1 = Fit1 / 2;
	Fit2 = Fit2 / 2;
	Fit3 = std::sqrt((Fit1*Fit1 + Fit2*Fit2) / 4 + Fit3/4.11);

	Fit.at<float>(0) = Fit1;
	Fit.at<float>(1) = Fit2;
	Fit.at<float>(2) = Fit3;
	
	return Fit;
}
