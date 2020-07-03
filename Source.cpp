/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * License and Disclaimer                                                *
 *                                                                       *
 * Neither the authors of this code nor their employing institutes make  *
 * any representation or warranty, expressed or implied, regarding this  *
 * code or assume any liability for its use.                             *
 *                                                                       *
 * By using, copying, modifying or distributing the code (or any work    *
 * based on the code) you agree to acknowledge its use in resulting      *
 * scientific publications.                                              *
 *                                                                       *
 * Author : Isshi Nara, 03/07/2020  (isshi_nara@eis.hokudai.ac.jp)       *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#define _USE_MATH_DEFINES
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cfloat>
#include <cmath>
#include <algorithm>
#include <tuple>
#include <map>
#include <chrono>
#include <iomanip>

using namespace std;

/* Alias for 2D-image */
using IMG = vector<vector<double>>;

/* Alias for paired data(unsigned int).
   first: x-axis data, second: y-axis data.
   (Used for Image size etc...) */
using uiP = pair<unsigned int, unsigned int>;

/* Alias for paired data(int).
   first: x-axis data, second: y-axis data.
   (Used for src search range; srcxLoop and srcyLoop variables) */
using iP = pair<int, int>;

/* Alias for paired data(double).
   first: x-axis data, second: y-axis data.
   (Used for Resolution, Isocenter etc...) */
using dP = pair<double, double>;

/* Alias for line equation(ax + by + c = 0).
   <0>: a, <1>: b, <2>: c */
using dT = tuple<double, double, double>;

class AreaAverageInterpolation
{
public:
	pair<bool, string> areaAverageInterpolation(IMG src, IMG &dst,
		dP srcResolution, dP dstResolution, dP srcIsocenter, dP &dstIsocenter, 
		double rotationAngle)
	{
		cout << "**********************************************************" << endl;
		cout << "* AreaAverageInterpolation::areaAverageInterpolation     *" << endl;
		cout << setprecision(10);
		cout << "* Input parameters                                       *" << endl;
		cout << "*                                                        *" << endl;
		cout << "* srcResolution : " << setw(9) << srcResolution.first << ", "
			                         << setw(9) << srcResolution.second
			                         << setw(20) << " [pixel/mm or dpi] *" << endl;
		cout << "* dstResolution : " << setw(9) << dstResolution.first << ", "
			                         << setw(9) << dstResolution.second
			                         << setw(20) << " [pixel/mm or dpi] *" << endl;
		cout << "* srcIsocenter  : " << setw(9) << srcIsocenter.first << ", "
			                         << setw(9) << srcIsocenter.second
			                         << setw(20) << " [pixels] *" << endl;
		cout << "* rotationAngle : " << setw(20) << rotationAngle
			                         << setw(20) << " [degrees] *" << endl;
		cout << "**********************************************************" << endl;


		/* #######   Description of arguments   #######
		 * src :                          source image
		 * dst :                          destination (interpolated and rotated) image
		 * srcResolution, dstResolution : src & dst resolution [pixel/mm or dpi]
		 * srcIsocenter, dstIsocenter :   src & dst isocenter (rotation center) [pixel]
		 * rotationAngle :                clockwise is positive [degree]
		 * Return argument :
		 *   first(bool)    -> is finished successfully
		 *   second(string) -> error message when first(bool) is false
		 */


		/* #######   Declarations   ####### */
		pair<bool, string> ret;     /* return argument */
		IMG modSrc;                 /* expanded source image by integer */
		uiP srcSize, modSrcSize;    /* image size of src & modSrc */
		uiP dstSize;                /* image size of dst */
		dP dstIsocenterOffset;      /* offset factor that moves dstIsocenter to the pixel center */
		vector<vector<dP>> dstPos;  /* dst pixel positions in src pixel scale */
		double expansionRatio;      /* expansion ratio (= dstResoluton / srcResolution) */
		double dstSideLength;       /* = 1/expansionRatio (the dst side length in src pixel scale) */
		vector<dT> dstHorizontalLine, dstVerticalLine; /* horizontal and vertical dst pixel sizes (expressed as ax + by + c = 0: tuple components are a, b and c.) */
		double sinValue, cosValue;
		dP offset;                  /* positional offset to prevent the image from being cut off after rotation */
		PixelState srcPixelState;   /* structure to determine area type */
		dP dstVertex[4], srcVertex[4];
		int intersectionType[4];    /* intersection type (See 'getIntersectionType' function for the details) */
		double r[4], s[4];          /* intersection point (See 'getIntersectionType' function for the details) */
		iP srcxLoop, srcyLoop;
		double tmpArea, sumArea;
		double areaWeightedPixelValue;


		/* #######   Judge whether arguments are adequate or not   ####### */
		if ( DBL_EPSILON < fabs(srcResolution.first - srcResolution.second) ||
			 DBL_EPSILON < fabs(dstResolution.first - dstResolution.second) ) {
			ret.first = false;
			ret.second = "Assumed X & Y resolution are same.";
			return ret;
		}
		if ( srcResolution.first <= DBL_EPSILON || dstResolution.first <= DBL_EPSILON ) {
			ret.first = false;
			ret.second = "0 or negative resolution is not acceptable.";
			return ret;
		}
		if ( src.size() == 0 ) {
			ret.first = false;
			ret.second = "There is no data in src array.";
			return ret;
		}
		if ( src.front().size() == 0 ) {
			ret.first = false;
			ret.second = "There is no data in the second dimension of src array.";
			return ret;
		}


		/* #######   Expand src image & Rotate 90, 180 or 270 degrees in advance   #######
		 *  - When srcResolution/sqrt(2) <= dstResolution, src image should be expanded, otherwise area cannot be calculated only 9 patterns (required 19 patterns).
		 *  - And to minimize branch of processing, rotation angle is restricted in the first quadrant by rotating 90, 180 or 270 degrees in advance.
		 */
		unsigned int scale = static_cast<unsigned int>(dstResolution.first / srcResolution.first * sqrt(2) + 1 + DBL_EPSILON);
		int beforehandRotationMode; /* 0:None, 1:90deg, 2:180deg, 3:270deg */
		while ( rotationAngle < 0 )		rotationAngle += 360; /* the range is restricted 0 to 360 degrees */
		while ( 360 <= rotationAngle )	rotationAngle -= 360; /* the range is restricted 0 to 360 degrees */
		if ( rotationAngle < 90 )       { beforehandRotationMode = 0; }
		else if ( rotationAngle < 180 ) { beforehandRotationMode = 1; rotationAngle -= 90; }
		else if ( rotationAngle < 270 ) { beforehandRotationMode = 2; rotationAngle -= 180; }
		else                            { beforehandRotationMode = 3; rotationAngle -= 270; }
		sinValue = sin(rotationAngle / 180.0*M_PI);
		cosValue = cos(rotationAngle / 180.0*M_PI);

		srcSize = make_pair(src.front().size(), src.size());
		if (beforehandRotationMode == 0 || beforehandRotationMode == 2) {
			modSrcSize = make_pair(srcSize.first*scale, srcSize.second*scale);
		}
		else {
			modSrcSize = make_pair(srcSize.second*scale, srcSize.first*scale);
		}
		modSrc.resize(modSrcSize.second);
		for ( unsigned int y = 0; y < modSrcSize.second; ++y )	modSrc[y].resize(modSrcSize.first);
		for ( unsigned int srcy = 0; srcy < srcSize.second; ++srcy ) {
			for ( unsigned int srcx = 0; srcx < srcSize.first; ++srcx ) {
				for ( unsigned int mody = 0; mody < scale; ++mody ) {
					for ( unsigned int modx = 0; modx < scale; ++modx ) {
						switch ( beforehandRotationMode ) {
						case 0:	modSrc[srcy*scale + mody][srcx*scale + modx] = src[srcy][srcx]; break;
						case 1:	modSrc[srcx*scale + modx][modSrcSize.first - 1 - (srcy*scale + mody)] = src[srcy][srcx]; break;
						case 2:	modSrc[modSrcSize.second - 1 - (srcy*scale + mody)][modSrcSize.first - 1 - (srcx*scale + modx)] = src[srcy][srcx]; break;
						case 3:	modSrc[modSrcSize.second - 1 - (srcx*scale + modx)][srcy*scale + mody] = src[srcy][srcx]; break;
						}
					}
				}
			}
		}
		srcIsocenter.first = srcIsocenter.first*scale + (scale - 1) / 2.0;
		srcIsocenter.second = srcIsocenter.second*scale + (scale - 1) / 2.0;
		srcResolution.first *= scale;
		srcResolution.second *= scale;
		expansionRatio = dstResolution.first / srcResolution.first;
		dstSideLength = srcResolution.first / dstResolution.first;
		dstSize.first = (unsigned int)round((modSrcSize.first*fabs(cosValue) + modSrcSize.second*fabs(sinValue))*expansionRatio);
		dstSize.second = (unsigned int)round((modSrcSize.first*fabs(sinValue) + modSrcSize.second*fabs(cosValue))*expansionRatio);
		dstIsocenter.first = (srcIsocenter.first*cosValue + (modSrcSize.second - srcIsocenter.second)*sinValue)*expansionRatio;
		dstIsocenter.second = (srcIsocenter.first*sinValue + srcIsocenter.second*cosValue)*expansionRatio;
		dstIsocenterOffset.first = dstIsocenter.first - int(dstIsocenter.first);
		dstIsocenterOffset.second = dstIsocenter.second - int(dstIsocenter.second);
		dstIsocenter.first = (int)dstIsocenter.first;
		dstIsocenter.second = (int)dstIsocenter.second;
		offset.first = offset.second = 0;

		/* Rotate left-upper point (0, 0) */
		offset.first = min(offset.first, -srcIsocenter.first*cosValue + srcIsocenter.second*sinValue + srcIsocenter.first);
		offset.second = min(offset.second, -srcIsocenter.first*sinValue - srcIsocenter.second*cosValue + srcIsocenter.second);
		/* Rotate right-upper (modSrcSize.first - 1, 0) */
		offset.first = min(offset.first, (modSrcSize.first - 1 - srcIsocenter.first)*cosValue + srcIsocenter.second*sinValue + srcIsocenter.first);
		offset.second = min(offset.second, (modSrcSize.first - 1 - srcIsocenter.first)*sinValue - srcIsocenter.second*cosValue + srcIsocenter.second);
		/* Rotate left-lower point (0, modSrcSize.second - 1) */
		offset.first = min(offset.first, -srcIsocenter.first*cosValue - (modSrcSize.second - 1 - srcIsocenter.second)*sinValue + srcIsocenter.first);
		offset.second = min(offset.second, -srcIsocenter.first*sinValue + (modSrcSize.second - 1 - srcIsocenter.second)*cosValue + srcIsocenter.second);
		/* Rotate right-lower (modSrcSize.first, modSrcSize.second) */
		offset.first = min(offset.first, (modSrcSize.first - 1 - srcIsocenter.first)*cosValue - (modSrcSize.second - 1 - srcIsocenter.second)*sinValue + srcIsocenter.first);
		offset.second = min(offset.second, (modSrcSize.first - 1 - srcIsocenter.first)*sinValue + (modSrcSize.second - 1 - srcIsocenter.second)*cosValue + srcIsocenter.second);


		/* #######   Calculate dst pixel position in src pixel scale   #######
		 *  - RotationAngle indicates angles from src to dst (clockwise is positive), however on this situation,
		 *    it should be rotated from dst pixel position to src pixel scale (this is common in image interpolation field).
		 *    For this reason, rotated inversely.
		 */
		dstPos.resize(dstSize.second);
		for ( unsigned int dsty = 0; dsty < dstSize.second; ++dsty ) {
			dstPos[dsty].resize(dstSize.first);
			for ( unsigned int dstx = 0; dstx < dstSize.first; ++dstx ) {
				dstPos[dsty][dstx].first =
					((dstx + dstIsocenterOffset.first)*dstSideLength - srcIsocenter.first + offset.first)*cosValue
					+ ((dsty + dstIsocenterOffset.second)*dstSideLength - srcIsocenter.second + offset.second)*sinValue
					+ srcIsocenter.first;
				dstPos[dsty][dstx].second =
					-((dstx + dstIsocenterOffset.first)*dstSideLength - srcIsocenter.first + offset.first)*sinValue
					+ ((dsty + dstIsocenterOffset.second)*dstSideLength - srcIsocenter.second + offset.second)*cosValue
					+ srcIsocenter.second;
			}
		}


		/* #######   Dst horizontal & vertical sides   #######
		 *  - On Area average interpolation, it will be calculated 'Area' of src pixel included dst pixel.
		 *    To calculate 'Area' of src pixels, the positions of the edges of the dst pixels need to be calculated.
		 *  - To minimize floating point error, the process is categorized into <45 and 45<= degrees cases.
		 */
		double tmpSinValue, tmpCosValue, tmpTanValue;
		if ( rotationAngle < 45 ) {
			tmpSinValue = sinValue;
			tmpCosValue = cosValue;
			tmpTanValue = tan((rotationAngle) / 180.0*M_PI);
		}
		else {
			tmpSinValue = sin((rotationAngle - 90) / 180.0*M_PI);
			tmpCosValue = cos((rotationAngle - 90) / 180.0*M_PI);
			tmpTanValue = tan((rotationAngle - 90) / 180.0*M_PI);
		}
		if ( fabs(tmpTanValue) < DBL_EPSILON )	tmpTanValue = 0;

		/* Horizontal sides */
		dstHorizontalLine.resize(dstSize.second + 1);
		for ( unsigned int dsty = 0; dsty < dstSize.second + 1; ++dsty ) {
			if ( rotationAngle < 45 ) {
				get<0>(dstHorizontalLine[dsty]) = tmpTanValue; /* ax + y + c = 0 (y = -ax - c) */
				get<1>(dstHorizontalLine[dsty]) = 1;
				if ( dsty < dstSize.second ) {
					get<2>(dstHorizontalLine[dsty]) =
						-get<0>(dstHorizontalLine[dsty])*(dstPos[dsty][0].first - dstSideLength / 2 * (tmpCosValue + tmpSinValue))
						- (dstPos[dsty][0].second - dstSideLength / 2 * (tmpCosValue - tmpSinValue));  /* c = -ax - y */
				}
				else {
					get<2>(dstHorizontalLine[dsty]) =
						-get<0>(dstHorizontalLine[dsty])*(dstPos.back()[0].first - dstSideLength / 2 * (tmpCosValue - tmpSinValue))
						- (dstPos.back()[0].second + dstSideLength / 2 * (tmpCosValue + tmpSinValue)); /* c = -ax - y */
				}
			}
			else {
				get<0>(dstHorizontalLine[dsty]) = 1;
				get<1>(dstHorizontalLine[dsty]) = -tmpTanValue; /* x + by + c = 0 (x = -by - c) */
				if ( dsty < dstSize.second ) {
					get<2>(dstHorizontalLine[dsty]) =
						-(dstPos[dsty][0].first - dstSideLength / 2 * (tmpCosValue + tmpSinValue))
						- get<1>(dstHorizontalLine[dsty])*(dstPos[dsty][0].second - dstSideLength / 2 * (tmpCosValue - tmpSinValue));  /* c = -x - by */
				}
				else {
					get<2>(dstHorizontalLine[dsty]) =
						-(dstPos.back()[0].first + dstSideLength / 2 * (tmpCosValue - tmpSinValue))
						- get<1>(dstHorizontalLine[dsty])*(dstPos.back()[0].second - dstSideLength / 2 * (tmpCosValue + tmpSinValue)); /* c = -x - by */
				}
			}
		}
		/* Vertical sides */
		dstVerticalLine.resize(dstSize.first + 1);
		for ( unsigned int dstx = 0; dstx < dstSize.first + 1; ++dstx ) {
			if ( rotationAngle < 45 ) {
				get<0>(dstVerticalLine[dstx]) = 1;
				get<1>(dstVerticalLine[dstx]) = -tmpTanValue; /* x + by + c = 0 (x = -by - c) */
				if ( dstx < dstSize.first ) {
					get<2>(dstVerticalLine[dstx]) =
						-(dstPos[0][dstx].first - dstSideLength / 2 * (tmpCosValue + tmpSinValue))
						- get<1>(dstVerticalLine[dstx])*(dstPos[0][dstx].second - dstSideLength / 2 * (tmpCosValue - tmpSinValue));  /* c = -x - by */
				}
				else {
					get<2>(dstVerticalLine[dstx]) =
						-(dstPos[0].back().first + dstSideLength / 2 * (tmpCosValue - tmpSinValue))
						- get<1>(dstVerticalLine[dstx])*(dstPos[0].back().second - dstSideLength / 2 * (tmpCosValue + tmpSinValue)); /* c = -x - by */
				}
			}
			else {
				get<0>(dstVerticalLine[dstx]) = tmpTanValue; /* ax + y + c = 0 (y = -ax - c) */
				get<1>(dstVerticalLine[dstx]) = 1;
				if ( dstx < dstSize.first ) {
					get<2>(dstVerticalLine[dstx]) =
						-get<0>(dstVerticalLine[dstx])*(dstPos[0][dstx].first - dstSideLength / 2 * (tmpCosValue - tmpSinValue))
						- (dstPos[0][dstx].second + dstSideLength / 2 * (tmpCosValue + tmpSinValue));  /* c = -ax - y */
				}
				else {
					get<2>(dstVerticalLine[dstx]) =
						-get<0>(dstVerticalLine[dstx])*(dstPos[0].back().first - dstSideLength / 2 * (tmpCosValue + tmpSinValue))
						- (dstPos[0].back().second - dstSideLength / 2 * (tmpCosValue - tmpSinValue)); /* c = -ax - y */
				}
			}
		}

		/* #######   Area average interpolation   #######
		 *  - Main calculation.
		 */
		auto initPixelState = [&]() {
			/* lambda function for initializing PixelState structure */
			srcPixelState.type = 0;
			srcPixelState.isIncludedSrcPixelCenter = false;
			srcPixelState.isIncludedDstPixelVertex = false;
			srcPixelState.vertexPos = make_pair(-1, -1);
			srcPixelState.intersections.clear();
			srcPixelState.intersections.insert(make_pair("xa", vector<double>()));
			srcPixelState.intersections.insert(make_pair("xb", vector<double>()));
			srcPixelState.intersections.insert(make_pair("ya", vector<double>()));
			srcPixelState.intersections.insert(make_pair("yb", vector<double>()));
			srcPixelState.xCounts = 0;
			srcPixelState.yCounts = 0;
		};
		auto updatePixelState_intersection = [&]() {
			/* lambda function for updating points where src pixel sides and dst pixel sides intersect */

			/* Boundary condition
			 * - Process to exclude point of contact with the dst side not penetrating the src pixel
			 */
			bool isOutOfSrcPixel = false;
			for ( unsigned int i = 0; i < 4; ++i ) {
				if ( intersectionType[i] == 4 ) {
					isOutOfSrcPixel = true;
					for ( unsigned int j = 0; j < 4; ++j ) {
						if ( i == j )	continue;
						if ( intersectionType[j] == 3 || intersectionType[j] == 4 ) {
							isOutOfSrcPixel = false; break;
						}
					}
					if ( isOutOfSrcPixel )	return;
				}
			}

			/* Get intersection points */
			for ( unsigned int i = 0; i < 4; ++i ) {
				if ( intersectionType[i] == 3 || intersectionType[i] == 4 ) {
					switch ( i ) {
					case 0:
						/* the point intersected with side 'xa' (See PixelState structure for details) */
						srcPixelState.intersections["xa"].emplace_back(s[i]);
						break;
					case 1:
						/* the point intersected with side 'ya' (See PixelState structure for details) */
						srcPixelState.intersections["ya"].emplace_back(s[i]);
						break;
					case 2:
						/* the point intersected with side 'yb' (See PixelState structure for details) */
						srcPixelState.intersections["yb"].emplace_back(s[i]);
						break;
					case 3:
						/* the point intersected with side 'xb' (See PixelState structure for details) */
						srcPixelState.intersections["xb"].emplace_back(s[i]);
						break;
					}
				}
			}
		};
		auto updatePixelState_isIncludedSrcPixel = [&](dP srcCenter) {
			/* Modified ray casting algorithm
			 *  1. Counts how many times a ray, starting from the srcCenter point and 
			       going in four direction (upper, lower, left, right), intersects dstVertices edges.
			 *  2. If srcCenter is inside the dst pixel, all rays will be intersected 1 or more times.
			 *  3. If srcCenter is outside the dst pixel, at least one of the rays will not be intersected.
			 */

			 /* when looking dst vertices clockwise, these are in this order */
			vector<dP> dstVertices = { dstVertex[0], dstVertex[1], dstVertex[3], dstVertex[2] };
			double tmpr, tmps;
			int crossPoints;
			dP srcCenterRay;
			int addx[] = { 0, 0, -100, 100 };
			int addy[] = { -100, 100, 0, 0 };
			for ( unsigned int direction = 0; direction < 4; ++direction ) {
				/* direction indicates 0:upper, 1:lower, 2:left, 3:right */
				crossPoints = 0;
				srcCenterRay = make_pair(srcCenter.first + addx[direction], srcCenter.second + addy[direction]);
				for ( unsigned int i = 0; i < dstVertices.size(); ++i ) {
					getIntersectionType(srcCenter, srcCenterRay, tmpr, dstVertices[i], dstVertices[(i + 1) % dstVertices.size()], tmps);
					if ( - DBL_EPSILON < tmpr && - DBL_EPSILON < tmps && tmps < 1 + DBL_EPSILON )	crossPoints++;
				}
				if ( crossPoints == 0 ) {
					srcPixelState.isIncludedSrcPixelCenter = false;
					return;
				}
			}
			srcPixelState.isIncludedSrcPixelCenter = true;
			return;
		};
		auto updatePixelState_isIncludedDstVertex = [&]() {
			/* lambda function for updating whether src pixel includes dst vertex or not */
			for ( unsigned int i = 0; i < 4; ++i ) {
				if ( (srcVertex[0].first + DBL_EPSILON < dstVertex[i].first && dstVertex[i].first < srcVertex[1].first - DBL_EPSILON) &&
					 (srcVertex[0].second + DBL_EPSILON < dstVertex[i].second && dstVertex[i].second < srcVertex[2].second - DBL_EPSILON) ) {
					srcPixelState.isIncludedDstPixelVertex = true;
					srcPixelState.vertexPos.first = dstVertex[i].first - srcVertex[0].first;
					srcPixelState.vertexPos.second = dstVertex[i].second - srcVertex[0].second;
				}
			}
		};

		dst.clear();
		dst.resize(dstSize.second);
		for ( unsigned int dsty = 0; dsty < dstSize.second; ++dsty ) {
			dst[dsty].resize(dstSize.first);
			for ( unsigned int dstx = 0; dstx < dstSize.first; ++dstx ) {
				tmpArea = sumArea = areaWeightedPixelValue = 0;

				/* get dst vertices */
				getIntersectionPoint(dstHorizontalLine[dsty], dstVerticalLine[dstx], dstVertex[0]);
				getIntersectionPoint(dstHorizontalLine[dsty], dstVerticalLine[dstx + 1], dstVertex[1]);
				getIntersectionPoint(dstHorizontalLine[dsty + 1], dstVerticalLine[dstx], dstVertex[2]);
				getIntersectionPoint(dstHorizontalLine[dsty + 1], dstVerticalLine[dstx + 1], dstVertex[3]);
				

				/* get src pixel search range (take care not to exceed image range) */
				srcxLoop.first = max(0, (int)floor(dstPos[dsty][dstx].first - dstSideLength*sqrt(2) / 2 - 1));
				srcxLoop.second = min((int)ceil(dstPos[dsty][dstx].first + dstSideLength*sqrt(2) / 2 + 1), (int)modSrcSize.first - 1);
				srcyLoop.first = max(0, (int)floor(dstPos[dsty][dstx].second - dstSideLength*sqrt(2) / 2 - 1));
				srcyLoop.second = min((int)ceil(dstPos[dsty][dstx].second + dstSideLength*sqrt(2) / 2 + 1), (int)modSrcSize.second - 1);
				

				for ( int srcy = srcyLoop.first; srcy <= srcyLoop.second; ++srcy ) {
					for ( int srcx = srcxLoop.first; srcx <= srcxLoop.second; ++srcx ) {
						/* initialize PixelState structure */
						initPixelState();

						/* get src vertices */
						srcVertex[0] = make_pair(srcx - 0.5, srcy - 0.5);
						srcVertex[1] = make_pair(srcx + 0.5, srcy - 0.5);
						srcVertex[2] = make_pair(srcx - 0.5, srcy + 0.5);
						srcVertex[3] = make_pair(srcx + 0.5, srcy + 0.5);


						/* get intersection points of src pixel sides and dst pixel sides */
						// - horizontal line 1
						intersectionType[0] = getIntersectionType(dstVertex[0], dstVertex[1], r[0], srcVertex[0], srcVertex[1], s[0]);
						intersectionType[1] = getIntersectionType(dstVertex[0], dstVertex[1], r[1], srcVertex[0], srcVertex[2], s[1]);
						intersectionType[2] = getIntersectionType(dstVertex[0], dstVertex[1], r[2], srcVertex[1], srcVertex[3], s[2]);
						intersectionType[3] = getIntersectionType(dstVertex[0], dstVertex[1], r[3], srcVertex[2], srcVertex[3], s[3]);
						updatePixelState_intersection();
						// - horizontal line 2
						intersectionType[0] = getIntersectionType(dstVertex[2], dstVertex[3], r[0], srcVertex[0], srcVertex[1], s[0]);
						intersectionType[1] = getIntersectionType(dstVertex[2], dstVertex[3], r[1], srcVertex[0], srcVertex[2], s[1]);
						intersectionType[2] = getIntersectionType(dstVertex[2], dstVertex[3], r[2], srcVertex[1], srcVertex[3], s[2]);
						intersectionType[3] = getIntersectionType(dstVertex[2], dstVertex[3], r[3], srcVertex[2], srcVertex[3], s[3]);
						updatePixelState_intersection();
						// - vertical line 1
						intersectionType[0] = getIntersectionType(dstVertex[0], dstVertex[2], r[0], srcVertex[0], srcVertex[1], s[0]);
						intersectionType[1] = getIntersectionType(dstVertex[0], dstVertex[2], r[1], srcVertex[0], srcVertex[2], s[1]);
						intersectionType[2] = getIntersectionType(dstVertex[0], dstVertex[2], r[2], srcVertex[1], srcVertex[3], s[2]);
						intersectionType[3] = getIntersectionType(dstVertex[0], dstVertex[2], r[3], srcVertex[2], srcVertex[3], s[3]);
						updatePixelState_intersection();
						// - vertical line 2
						intersectionType[0] = getIntersectionType(dstVertex[1], dstVertex[3], r[0], srcVertex[0], srcVertex[1], s[0]);
						intersectionType[1] = getIntersectionType(dstVertex[1], dstVertex[3], r[1], srcVertex[0], srcVertex[2], s[1]);
						intersectionType[2] = getIntersectionType(dstVertex[1], dstVertex[3], r[2], srcVertex[1], srcVertex[3], s[2]);
						intersectionType[3] = getIntersectionType(dstVertex[1], dstVertex[3], r[3], srcVertex[2], srcVertex[3], s[3]);
						updatePixelState_intersection();

						/* does dst pixel include src pixel center ? */
						updatePixelState_isIncludedSrcPixel(make_pair(srcx, srcy));
						

						/* does src pixel include dst vertices ? */
						updatePixelState_isIncludedDstVertex();


						/* processing of boundary condition
						 *  - If a dst side touches src sides on src vertices,
						 *    it means there are two intersection points that are very close together (See ASCII art below).
						 *    Because this state is peculiar to the tangential boundary condition, one of two intersection points should be deleted.
						 *    Rule 1 : Delete the x-axis intersection point (remaining the y-axis intersection point) on this algorithm.
						 *  - If a dst side touches src sides on vertices under the condition 0, 90, 180 and 270 degrees rotation,
						 *    only ONE point can intersect at the src vertex. In this case, the point of intersection has no meaning.
						 *    Rule 2 : Delete the point.
						 * 
						 *  - Rule 1 -          /	                 - Rule 2 -
						 *                     / <- a dst side	    
						 *  two points here-> ∙──────┼─              -- ─∙──────∙─ ------ <- a dst side
						 *                   /│      │                   │      │
						 *                  / │      │                   │      │      bullet(∙) means intersection points
						 *                 / ─┼──────┼─                 ─┼──────┼─
						 */

						/* sort vector stored intersection points */
						for ( auto& a : srcPixelState.intersections )	sort(a.second.begin(), a.second.end());

						/* remove y-elements if it indicates an end point on the y-axis
						   and there is no point on the x-axis at the same position. */
						bool isExistSamePoint = false;
						auto yItr = srcPixelState.intersections["ya"].begin();
						while ( yItr != srcPixelState.intersections["ya"].end() ) {
							if ( *yItr <= DBL_EPSILON ) {
								isExistSamePoint = false;
								for ( double d : srcPixelState.intersections["xa"] ) {
									if ( d <= DBL_EPSILON ) {
										isExistSamePoint = true;
										break;
									}
								}
								if ( ! isExistSamePoint ) yItr = srcPixelState.intersections["ya"].erase(yItr);
								else ++yItr;
							}
							else if ( 1 - (*yItr) <= DBL_EPSILON ) {
								isExistSamePoint = false;
								for ( double d : srcPixelState.intersections["xb"] ) {
									if ( d <= DBL_EPSILON ) {
										isExistSamePoint = true;
										break;
									}
								}
								if ( ! isExistSamePoint ) yItr = srcPixelState.intersections["ya"].erase(yItr);
								else ++yItr;
							}
							else ++yItr;
						}
						yItr = srcPixelState.intersections["yb"].begin();
						while ( yItr != srcPixelState.intersections["yb"].end() ) {
							if ( *yItr <= DBL_EPSILON ) {
								isExistSamePoint = false;
								for ( double d : srcPixelState.intersections["xa"] ) {
									if ( 1 - d <= DBL_EPSILON ) {
										isExistSamePoint = true;
										break;
									}
								}
								if ( ! isExistSamePoint ) yItr = srcPixelState.intersections["yb"].erase(yItr);
								else ++yItr;
							}
							else if ( 1 - (*yItr) <= DBL_EPSILON ) {
								isExistSamePoint = false;
								for ( double d : srcPixelState.intersections["xb"] ) {
									if ( 1 - d <= DBL_EPSILON ) {
										isExistSamePoint = true;
										break;
									}
								}
								if ( ! isExistSamePoint ) yItr = srcPixelState.intersections["yb"].erase(yItr);
								else ++yItr;
							}
							else ++yItr;
						}

						/* remove elements that points to the end points on x-axis */
						auto xItr = srcPixelState.intersections["xa"].begin();
						while ( xItr != srcPixelState.intersections["xa"].end() ) {
							if ( *xItr <= DBL_EPSILON || 1 - (*xItr) <= DBL_EPSILON ) xItr = srcPixelState.intersections["xa"].erase(xItr);
							else ++xItr;
						}
						xItr = srcPixelState.intersections["xb"].begin();
						while ( xItr != srcPixelState.intersections["xb"].end() ) {
							if ( *xItr <= DBL_EPSILON || 1 - (*xItr) <= DBL_EPSILON ) xItr = srcPixelState.intersections["xb"].erase(xItr);
							else ++xItr;
						}

						
						/* count intersection points (does not exceed 255) */
						srcPixelState.xCounts = static_cast<unsigned char>(srcPixelState.intersections["xa"].size() + srcPixelState.intersections["xb"].size());
						srcPixelState.yCounts = static_cast<unsigned char>(srcPixelState.intersections["ya"].size() + srcPixelState.intersections["yb"].size());
						
						/* get area & calculate weighted pixel value */
						tmpArea = getArea(srcPixelState);
						sumArea += tmpArea;
						areaWeightedPixelValue += modSrc[srcy][srcx] * tmpArea;
					}
				}
				dst[dsty][dstx] = DBL_EPSILON < fabs(sumArea) ? areaWeightedPixelValue / sumArea : 0;
			}
		}
		ret.first = true;
		ret.second = "";
		return ret;
	}
	pair<bool, string> fastAreaAverageInterpolation(IMG src, IMG &dst,
		dP srcResolution, dP dstResolution, dP srcIsocenter, dP &dstIsocenter, 
		double rotationAngle)
	{
		cout << "**********************************************************" << endl;
		cout << "* AreaAverageInterpolation::fastAreaAverageInterpolation *" << endl;
		cout << setprecision(10);
		cout << "* Input parameters                                       *" << endl;
		cout << "*                                                        *" << endl;
		cout << "* srcResolution : " << setw(9) << srcResolution.first << ", "
			                         << setw(9) << srcResolution.second
			                         << setw(20) << " [pixel/mm or dpi] *" << endl;
		cout << "* dstResolution : " << setw(9) << dstResolution.first << ", "
			                         << setw(9) << dstResolution.second
			                         << setw(20) << " [pixel/mm or dpi] *" << endl;
		cout << "* srcIsocenter  : " << setw(9) << srcIsocenter.first << ", "
			                         << setw(9) << srcIsocenter.second
			                         << setw(20) << " [pixels] *" << endl;
		cout << "* rotationAngle : " << setw(20) << rotationAngle
			                         << setw(20) << " [degrees] *" << endl;
		cout << "**********************************************************" << endl;

		/* #######   Description of arguments   #######
		 * src :                          source image
		 * dst :                          destination (interpolated & rotated image) image
		 * srcResolution, dstResolution : source & destination resolution [pixel/mm or dpi]
		 * srcIsocenter, dstIsocenter :   source & destination isocenter (rotation center) [pixel]
		 * rotationAngle :                clockwise is positive [degree]
		 * Return argument :
		 *   first(bool)    -> is finished successfully
		 *   second(string) -> error message when first(bool) is false
		 */


		/* #######   Declarations   ####### */
		pair<bool, string> ret;     /* return argument */
		IMG modSrc;                 /* expanded source image by integer */
		uiP srcSize, modSrcSize;    /* image size of src & modSrc */
		uiP dstSize;                /* image size of dst */
		dP dstIsocenterOffset;      /* offset factor that moves dstIsocenter to the pixel center */
		vector<vector<dP>> dstPos;  /* dst pixel positions in src pixel scale */
		double expansionRatio;      /* Expansion ratio (= dstResoluton / srcResolution) */
		double dstSideLength;       /* = 1/expansionRatio (the dst side length in src pixel scale) */
		vector<dT> dstHorizontalLine, dstVerticalLine; /* horizontal and vertical dst pixel sizes (expressed as ax + by + c = 0: tuple components are a, b and c.) */
		double _sin, _cos;
		dP offset;                  /* Positional offset to prevent the image from being cut off after rotation */
		bool isIncludedSrcPixelCenter;
		dP dstVertex[4], srcVertex[4];
		iP srcxLoop, srcyLoop;
		int sumArea;
		double areaWeightedPixelValue;


		/* #######   Judge whether arguments are adequate or not   ####### */
		if ( DBL_EPSILON < fabs(srcResolution.first - srcResolution.second) ||
			 DBL_EPSILON < fabs(dstResolution.first - dstResolution.second) ) {
			ret.first = false;
			ret.second = "Assumed X & Y resolution are same.";
			return ret;
		}
		if ( srcResolution.first <= DBL_EPSILON || dstResolution.first <= DBL_EPSILON ) {
			ret.first = false;
			ret.second = "0 or negative resolution is not acceptable.";
			return ret;
		}
		if ( src.size() == 0 ) {
			ret.first = false;
			ret.second = "There is no data in src array.";
			return ret;
		}
		if ( src.front().size() == 0 ) {
			ret.first = false;
			ret.second = "There is no data in the second dimension of src array.";
			return ret;
		}


		/* #######   Expand src image & Rotate 90, 180 or 270 degrees in advance   #######
		 *  - When srcResolution/sqrt(2) <= dstResolution, src image should be expanded, otherwise area cannot be calculated only 9 patterns (required 19 patterns).
		 *  - And to minimize branch of processing, rotation angle is restricted in the first quadrant by rotating 90, 180 or 270 degrees in advance.
		 */
		unsigned int scale = static_cast<unsigned int>(dstResolution.first / srcResolution.first * sqrt(2) + 1 + DBL_EPSILON);
		int beforehandRotationMode; /* 0:None, 1:90deg, 2:180deg, 3:270deg */
		while ( rotationAngle < 0 )     rotationAngle += 360; /* the range is restricted 0 to 360 degrees */
		while ( 360 <= rotationAngle )  rotationAngle -= 360; /* the range is restricted 0 to 360 degrees */
		if ( rotationAngle < 90 )       { beforehandRotationMode = 0; }
		else if ( rotationAngle < 180 ) { beforehandRotationMode = 1; rotationAngle -= 90; }
		else if ( rotationAngle < 270 ) { beforehandRotationMode = 2; rotationAngle -= 180; }
		else                            { beforehandRotationMode = 3; rotationAngle -= 270; }
		_sin = sin(rotationAngle / 180.0*M_PI);
		_cos = cos(rotationAngle / 180.0*M_PI);

		srcSize = make_pair(src.front().size(), src.size());
		if (beforehandRotationMode == 0 || beforehandRotationMode == 2) {
			modSrcSize = make_pair(srcSize.first*scale, srcSize.second*scale);
		}
		else {
			modSrcSize = make_pair(srcSize.second*scale, srcSize.first*scale);
		}
		modSrc.resize(modSrcSize.second);
		for ( unsigned int y = 0; y < modSrcSize.second; ++y )	modSrc[y].resize(modSrcSize.first);
		for ( unsigned int srcy = 0; srcy < srcSize.second; ++srcy ) {
			for ( unsigned int srcx = 0; srcx < srcSize.first; ++srcx ) {
				for ( unsigned int mody = 0; mody < scale; ++mody ) {
					for ( unsigned int modx = 0; modx < scale; ++modx ) {
						switch ( beforehandRotationMode ) {
						case 0:	modSrc[srcy*scale + mody][srcx*scale + modx] = src[srcy][srcx]; break;
						case 1:	modSrc[srcx*scale + modx][modSrcSize.first - 1 - (srcy*scale + mody)] = src[srcy][srcx]; break;
						case 2:	modSrc[modSrcSize.second - 1 - (srcy*scale + mody)][modSrcSize.first - 1 - (srcx*scale + modx)] = src[srcy][srcx]; break;
						case 3:	modSrc[modSrcSize.second - 1 - (srcx*scale + modx)][srcy*scale + mody] = src[srcy][srcx]; break;
						}
					}
				}
			}
		}
		srcIsocenter.first = srcIsocenter.first*scale + (scale - 1) / 2.0;
		srcIsocenter.second = srcIsocenter.second*scale + (scale - 1) / 2.0;
		srcResolution.first *= scale;
		srcResolution.second *= scale;
		expansionRatio = dstResolution.first / srcResolution.first;
		dstSideLength = srcResolution.first / dstResolution.first;
		dstSize.first = (unsigned int)round((modSrcSize.first*fabs(_cos) + modSrcSize.second*fabs(_sin))*expansionRatio);
		dstSize.second = (unsigned int)round((modSrcSize.first*fabs(_sin) + modSrcSize.second*fabs(_cos))*expansionRatio);
		dstIsocenter.first = (srcIsocenter.first*_cos + (modSrcSize.second - srcIsocenter.second)*_sin)*expansionRatio;
		dstIsocenter.second = (srcIsocenter.first*_sin + srcIsocenter.second*_cos)*expansionRatio;
		dstIsocenterOffset.first = dstIsocenter.first - int(dstIsocenter.first);
		dstIsocenterOffset.second = dstIsocenter.second - int(dstIsocenter.second);
		dstIsocenter.first = (int)dstIsocenter.first;
		dstIsocenter.second = (int)dstIsocenter.second;
		offset.first = offset.second = 0;

		/* Rotate left-upper point (0, 0) */
		offset.first = min(offset.first, -srcIsocenter.first*_cos + srcIsocenter.second*_sin + srcIsocenter.first);
		offset.second = min(offset.second, -srcIsocenter.first*_sin - srcIsocenter.second*_cos + srcIsocenter.second);
		/* Rotate right-upper (modSrcSize.first - 1, 0) */
		offset.first = min(offset.first, (modSrcSize.first - 1 - srcIsocenter.first)*_cos + srcIsocenter.second*_sin + srcIsocenter.first);
		offset.second = min(offset.second, (modSrcSize.first - 1 - srcIsocenter.first)*_sin - srcIsocenter.second*_cos + srcIsocenter.second);
		/* Rotate left-lower point (0, modSrcSize.second - 1) */
		offset.first = min(offset.first, -srcIsocenter.first*_cos - (modSrcSize.second - 1 - srcIsocenter.second)*_sin + srcIsocenter.first);
		offset.second = min(offset.second, -srcIsocenter.first*_sin + (modSrcSize.second - 1 - srcIsocenter.second)*_cos + srcIsocenter.second);
		/* Rotate right-lower (modSrcSize.first, modSrcSize.second) */
		offset.first = min(offset.first, (modSrcSize.first - 1 - srcIsocenter.first)*_cos - (modSrcSize.second - 1 - srcIsocenter.second)*_sin + srcIsocenter.first);
		offset.second = min(offset.second, (modSrcSize.first - 1 - srcIsocenter.first)*_sin + (modSrcSize.second - 1 - srcIsocenter.second)*_cos + srcIsocenter.second);


		/* #######   Calculate dst pixel position in src pixel scale   #######
		 *  - RotationAngle indicates angles from src to dst (clockwise is positive), however on this situation,
		 *    it should be rotated from dst pixel position to src pixel scale (this is common in image interpolation field).
		 *    Then, rotated inversely.
		 */
		dstPos.resize(dstSize.second);
		for ( unsigned int dsty = 0; dsty < dstSize.second; ++dsty ) {
			dstPos[dsty].resize(dstSize.first);
			for ( unsigned int dstx = 0; dstx < dstSize.first; ++dstx ) {
				dstPos[dsty][dstx].first = 
					((dstx + dstIsocenterOffset.first)*dstSideLength - srcIsocenter.first + offset.first)*_cos
					+ ((dsty + dstIsocenterOffset.second)*dstSideLength - srcIsocenter.second + offset.second)*_sin
					+ srcIsocenter.first;
				dstPos[dsty][dstx].second =
					-((dstx + dstIsocenterOffset.first)*dstSideLength - srcIsocenter.first + offset.first)*_sin
					+ ((dsty + dstIsocenterOffset.second)*dstSideLength - srcIsocenter.second + offset.second)*_cos
					+ srcIsocenter.second;
			}
		}


		/* #######   Dst horizontal & vertical sides   #######
		 *  - On Fast area average interpolation, use the information about whether dst pixel includes src pixel.
		 *    To judge it, the positions of the edges of the dst pixels need to be calculated.
		 *  - To minimize floating point error, the process is categorized into <45 and 45<= degrees cases.
		 */
		double tmp_sin, tmp_cos, tmp_tan;
		if ( rotationAngle < 45 ) {
			tmp_sin = _sin;
			tmp_cos = _cos;
			tmp_tan = tan((rotationAngle) / 180.0*M_PI);
		}
		else {
			tmp_sin = sin((rotationAngle - 90) / 180.0*M_PI);
			tmp_cos = cos((rotationAngle - 90) / 180.0*M_PI);
			tmp_tan = tan((rotationAngle - 90) / 180.0*M_PI);
		}
		if ( fabs(tmp_tan) < DBL_EPSILON )	tmp_tan = 0;

		/* Horizontal sides */
		dstHorizontalLine.resize(dstSize.second + 1);
		for ( unsigned int dsty = 0; dsty < dstSize.second + 1; ++dsty ) {
			if ( rotationAngle < 45 ) {
				get<0>(dstHorizontalLine[dsty]) = tmp_tan; /* ax + y + c = 0 (y = -ax - c) */
				get<1>(dstHorizontalLine[dsty]) = 1;
				if ( dsty < dstSize.second ) {
					get<2>(dstHorizontalLine[dsty]) =
						-get<0>(dstHorizontalLine[dsty])*(dstPos[dsty][0].first - dstSideLength / 2 * (tmp_cos + tmp_sin))
						- (dstPos[dsty][0].second - dstSideLength / 2 * (tmp_cos - tmp_sin)); /* c = -ax - y */
				}
				else {
					get<2>(dstHorizontalLine[dsty]) =
						-get<0>(dstHorizontalLine[dsty])*(dstPos.back()[0].first - dstSideLength / 2 * (tmp_cos - tmp_sin))
						- (dstPos.back()[0].second + dstSideLength / 2 * (tmp_cos + tmp_sin)); /* c = -ax - y */
				}
			}
			else {
				get<0>(dstHorizontalLine[dsty]) = 1;
				get<1>(dstHorizontalLine[dsty]) = -tmp_tan; /* x + by + c = 0 (x = -by - c) */
				if ( dsty < dstSize.second ) {
					get<2>(dstHorizontalLine[dsty]) =
						-(dstPos[dsty][0].first - dstSideLength / 2 * (tmp_cos + tmp_sin))
						- get<1>(dstHorizontalLine[dsty])*(dstPos[dsty][0].second - dstSideLength / 2 * (tmp_cos - tmp_sin)); /* c = -x - by */
				}
				else {
					get<2>(dstHorizontalLine[dsty]) =
						-(dstPos.back()[0].first + dstSideLength / 2 * (tmp_cos - tmp_sin))
						- get<1>(dstHorizontalLine[dsty])*(dstPos.back()[0].second - dstSideLength / 2 * (tmp_cos + tmp_sin)); /* c = -x - by */
				}
			}
		}
		// Vertical
		dstVerticalLine.resize(dstSize.first + 1);
		for ( unsigned int dstx = 0; dstx < dstSize.first + 1; ++dstx ) {
			if ( rotationAngle < 45 ) {
				get<0>(dstVerticalLine[dstx]) = 1;
				get<1>(dstVerticalLine[dstx]) = -tmp_tan; /* x + by + c = 0 (x = -by - c) */
				if ( dstx < dstSize.first ) {
					get<2>(dstVerticalLine[dstx]) =
						-(dstPos[0][dstx].first - dstSideLength / 2 * (tmp_cos + tmp_sin))
						- get<1>(dstVerticalLine[dstx])*(dstPos[0][dstx].second - dstSideLength / 2 * (tmp_cos - tmp_sin)); /* c = -x - by */
				}
				else {
					get<2>(dstVerticalLine[dstx]) =
						-(dstPos[0].back().first + dstSideLength / 2 * (tmp_cos - tmp_sin))
						- get<1>(dstVerticalLine[dstx])*(dstPos[0].back().second - dstSideLength / 2 * (tmp_cos + tmp_sin)); /* c = -x - by */
				}
			}
			else {
				get<0>(dstVerticalLine[dstx]) = tmp_tan; /* ax + y + c = 0 (y = -ax - c) */
				get<1>(dstVerticalLine[dstx]) = 1;
				if ( dstx < dstSize.first ) {
					get<2>(dstVerticalLine[dstx]) =
						-get<0>(dstVerticalLine[dstx])*(dstPos[0][dstx].first - dstSideLength / 2 * (tmp_cos - tmp_sin))
						- (dstPos[0][dstx].second + dstSideLength / 2 * (tmp_cos + tmp_sin)); /* c = -ax - y */
				}
				else {
					get<2>(dstVerticalLine[dstx]) =
						-get<0>(dstVerticalLine[dstx])*(dstPos[0].back().first - dstSideLength / 2 * (tmp_cos + tmp_sin))
						- (dstPos[0].back().second - dstSideLength / 2 * (tmp_cos - tmp_sin)); /* c = -ax - y */
				}
			}
		}


		/* #######   Fast area average interpolation   #######
		 *  - Main calculation.
		 */
		auto isIncludedSrcPixel = [&](dP srcCenter)->bool {
			/* Modified ray casting algorithm
			 *  1. Counts how many times a ray, starting from the srcCenter point and going in four direction (upper, lower, left, right), intersects dstVertices edges.
			 *  2. If srcCenter is inside the dst pixel, all rays will be intersected 1 or more times.
			 *  3. If srcCenter is outside the dst pixel, at least one of the rays will not be intersected.
			 */

			/* when looking clockwise, dst vertices are in this order */
			vector<dP> dstVertices = { dstVertex[0], dstVertex[1], dstVertex[3], dstVertex[2] };
			double tmpr, tmps;
			int crossPoints;
			dP srcCenterRay;
			int addx[] = { 0, 0, -100, 100 };
			int addy[] = { -100, 100, 0, 0 };
			for ( unsigned int direction = 0; direction < 4; ++direction ) {
				/* direction == 0: upper, 1: lower, 2: left, 3: right */
				crossPoints = 0;
				srcCenterRay = make_pair(srcCenter.first + addx[direction], srcCenter.second + addy[direction]);
				for ( unsigned int i = 0; i < dstVertices.size(); ++i ) {
					getIntersectionType(srcCenter, srcCenterRay, tmpr, dstVertices[i], dstVertices[(i + 1) % dstVertices.size()], tmps);
					if ( - DBL_EPSILON < tmpr && - DBL_EPSILON < tmps && tmps < 1 + DBL_EPSILON )	crossPoints++;
				}
				if ( crossPoints == 0 ) {
					return false;
				}
			}
			return true;
		};

		dst.clear();
		dst.resize(dstSize.second);
		for ( unsigned int dsty = 0; dsty < dstSize.second; ++dsty ) {
			dst[dsty].resize(dstSize.first);
			for ( unsigned int dstx = 0; dstx < dstSize.first; ++dstx ) {
				sumArea = 0;
				areaWeightedPixelValue = 0;

				/* get dst vertices */
				getIntersectionPoint(dstHorizontalLine[dsty], dstVerticalLine[dstx], dstVertex[0]);
				getIntersectionPoint(dstHorizontalLine[dsty], dstVerticalLine[dstx + 1], dstVertex[1]);
				getIntersectionPoint(dstHorizontalLine[dsty + 1], dstVerticalLine[dstx], dstVertex[2]);
				getIntersectionPoint(dstHorizontalLine[dsty + 1], dstVerticalLine[dstx + 1], dstVertex[3]);
				

				/* Get src pixel search range (take care not to exceed image range) */
				srcxLoop.first = max(0, (int)floor(dstPos[dsty][dstx].first - dstSideLength*sqrt(2) / 2 - 1));
				srcxLoop.second = min((int)ceil(dstPos[dsty][dstx].first + dstSideLength*sqrt(2) / 2 + 1), (int)modSrcSize.first - 1);
				srcyLoop.first = max(0, (int)floor(dstPos[dsty][dstx].second - dstSideLength*sqrt(2) / 2 - 1));
				srcyLoop.second = min((int)ceil(dstPos[dsty][dstx].second + dstSideLength*sqrt(2) / 2 + 1), (int)modSrcSize.second - 1);
				

				for ( int srcy = srcyLoop.first; srcy <= srcyLoop.second; ++srcy ) {
					for ( int srcx = srcxLoop.first; srcx <= srcxLoop.second; ++srcx ) {
						/* get src vertices */
						srcVertex[0] = make_pair(srcx - 0.5, srcy - 0.5);
						srcVertex[1] = make_pair(srcx + 0.5, srcy - 0.5);
						srcVertex[2] = make_pair(srcx - 0.5, srcy + 0.5);
						srcVertex[3] = make_pair(srcx + 0.5, srcy + 0.5);

						/* does dst pixel include src pixel center ? */
						isIncludedSrcPixelCenter = isIncludedSrcPixel(make_pair(srcx, srcy));

						if ( isIncludedSrcPixelCenter ) {
							sumArea += 1;
							areaWeightedPixelValue += modSrc[srcy][srcx];
						}
					}
				}
				dst[dsty][dstx] = 0 < sumArea ? areaWeightedPixelValue / sumArea : 0;
			}
		}
		ret.first = true;
		ret.second = "";
		return ret;
	}

private:
	struct PixelState {
		/* #######   Structure for source pixel state   #######
		 *  - Store source pixel state for judgement area type and calculation area.
		 */

		unsigned char type = 0;
		/*
		 * 0  : not included
		 * 1  : whole pixel
		 * 2  : triangle
		 * 3  : quadrangle
		 * 4  : pentagon made by 1 line
		 * 5  : pentagon made by 2 lines
		 * 6  : hexagon
		 * 7  : triangle with interpolated vertex
		 * 8  : quadrangle with interpolated vertex
		 * 9  : pentagon with interpolated vertex
		 */

		 /* is included src pixel center in the dst pixel
		   (It will be judged 'false' when src pixel center overlaps dst pixel.) */
		bool isIncludedSrcPixelCenter = false;

		/* is included dst vertex in src pixel. judge type = 7~9 or not.
		  (It will be judged 'false' when src pixel center overlaps dst pixel.) */
		bool isIncludedDstPixelVertex = false;		

		/* dst vertex position (at most only one dst vertex is included because of dstResolution < srcResolution/sqrt(2) constraint) */
		pair<double, double> vertexPos; 
		
		/* intersection points of src pixel and dst pixel */
		map<string, vector<double>> intersections;
		/* key:   "xa","xb","ya" or "yb" (coordinate system is shown below)
		*  value: the range of intersection points is 0 to 1.
		* ┌─────────────── x
		* │   ya    yb
		* │  ─┼─────┼─xa
		* │   │     │
		* │   │     │
		* │  ─┼─────┼─xb
		* │
		* y
		*/

		/* number of intersection points for x(xa and xb), y(ya and yb) sides. */
		unsigned char xCounts = 0, yCounts = 0;
	};

	bool getIntersectionPoint(dT line1, dT line2, dP &p)
	{
		/* get intersection point of two lines expressed as following equation : ax + by + c = 0 */
		double a1 = get<0>(line1);
		double b1 = get<1>(line1);
		double c1 = get<2>(line1);
		double a2 = get<0>(line2);
		double b2 = get<1>(line2);
		double c2 = get<2>(line2);

		if ( (fabs(a1) <= DBL_EPSILON && fabs(b1) <= DBL_EPSILON) || (fabs(a2) <= DBL_EPSILON && fabs(b2) <= DBL_EPSILON) )   return false;
		else if ( fabs(b1) <= DBL_EPSILON && fabs(b2) <= DBL_EPSILON )    return false; /* Parallel (Both lines are x=constant) */
		else if ( fabs(a1) <= DBL_EPSILON && fabs(a2) <= DBL_EPSILON )    return false; /* Parallel (Both lines are y=constant) */
		else if ( fabs(a2*b1 - a1*b2) <= DBL_EPSILON )  return false; /* Parallel */
		else if ( fabs(b2) <= DBL_EPSILON ) {
			p.first = -c2 / a2;					/* Confirmed a2 != 0 in advance */
			p.second = (a1*c2 - a2*c1) / a2*b1; /* Confirmed b1 != 0 in advance */
		}
		else {
			p.first = (b2*c1 - b1*c2) / (a2*b1 - a1*b2);
			p.second = (a1*c2 - a2*c1) / (a2*b1 - a1*b2);
		}
		return true;
	}
	int getIntersectionType(dP p1, dP p2, double &r, dP q1, dP q2, double &s)
	{
		/* description of arguments
		 * r : intersection point on p1-p2 line
		 * s : intersection point on q1-q2 line
		 * - <0  means intersecting out of line segment, close to p1(q1)
		 * - 0   means intersecting on p1(q1)
		 * - 0-1 means intersecting p1-p2(q1-q2) line segment
		 * - 1   means intersecting on p2(q2)
		 * - 1<  means intersecting out of line segment, close to p2(q2)
		 */

		int ret = 0;
		/* return argument
		 * 1 : do not intersect (parallel)
		 * 2 : completely overlap
		 * 3 : intersect on line segments (except at end point)
		 * 4 : intersect at end point
		 * 5 : intersect on straight lines (except on line segments)
		 */

		double denominator;
		double r_numerator, s_numerator;

		denominator = (p2.first - p1.first)*(q2.second - q1.second) - (p2.second - p1.second)*(q2.first - q1.first);
		r_numerator = (q1.first - p1.first)*(q2.second - q1.second) - (q1.second - p1.second)*(q2.first - q1.first);
		s_numerator = (p2.second - p1.second)*(q1.first - p1.first) - (p2.first - p1.first)*(q1.second - p1.second);

		if ( fabs(denominator) <= DBL_EPSILON && fabs(r_numerator) <= DBL_EPSILON && fabs(s_numerator) <= DBL_EPSILON ) {
			ret = 2;     /* completely overlap */
			return ret;
		}
		else if ( fabs(denominator) <= DBL_EPSILON ) {
			ret = 1;     /* do not intersect (parallel) */
			return ret;
		}

		r = r_numerator / denominator;
		s = s_numerator / denominator;

		if ( - DBL_EPSILON <= r && r <= 1.0 + DBL_EPSILON && - DBL_EPSILON <= s && s <= 1.0 + DBL_EPSILON ) {
			if ( fabs(r) <= DBL_EPSILON || fabs(r - 1.0) <= DBL_EPSILON || fabs(s) <= DBL_EPSILON || fabs(s - 1.0) <= DBL_EPSILON ) {
				ret = 4;  /* intersect at end point */
			}
			else ret = 3; /* intersect on line segments (except at end point) */
		}
		else ret = 5;     /* intersect on straight lines (except on line segments) */
		return ret;
	}
	double getArea(PixelState state)
	{
		/*
		 * xCount, yCount, include, vertex, type, type name
		 * 0,      0,      false,   false,  0,    Not included
		 * 0,      0,      true,    false,  1,    Whole pixel
		 * 1,      1,      false,   false,  2,    Triangle
		 * 2(0),   0(2),            false,  3,    Quadrangle
		 * 1,      1,      true,    false,  4,    Pentagon made by 1 line
		 * 3(1),   1(3),            false,  5,    Pentagon made by 2 lines
		 * 2,      2,      true,    false,  6,    Hexagon
		 * 2(0),   0(2),            true,   7,    Triangle with dst vertex (These parameters are same as type 9. See below for a classification of type 7 and 9.)
		 * 1,      1,               true,   8,    Quadrangle with dst vertex
		 * 2(0),   0(2),            true,   9,    Pentagon with dst vertex (These parameters are same as type 7. See below for a classification of type 7 and 9.)
		 *  - classification of type 7 and 9 : when xa1 & xa2 or xb1 & xb2 are stored -> type 7, and when xa1 & xb1 are stored -> type 9. Same for y direction.
		 */

		auto type1 = [&]()->double {
			return 1;
		};
		auto type2 = [&]()->double {
			double x, y;
			if ( state.intersections["xa"].size() != 0 )	x = state.intersections["xa"][0];
			else x = 1 - state.intersections["xb"][0];
			if ( state.intersections["ya"].size() != 0 )	y = state.intersections["ya"][0];
			else y = 1 - state.intersections["yb"][0];
			return 0.5*x*y;
		};
		auto type3 = [&]()->double {
			double ret;
			double side1 = -1, side2 = -1;
			if ( state.intersections["xa"].size() != 0 && state.intersections["xb"].size() != 0 ) {
				side1 = state.intersections["xa"][0];
				side2 = state.intersections["xb"][0];
			}
			else if ( state.intersections["ya"].size() != 0 && state.intersections["yb"].size() != 0 ) {
				side1 = state.intersections["ya"][0];
				side2 = state.intersections["yb"][0];
			}
			else {
				/* boundary condition
				   dst vertex was not be detected because dst vertex is on the src sides. then, if src pixel center is included return 1, otherwise return 0 */
				return state.isIncludedSrcPixelCenter ? 1 : 0;
			}
			ret = 0.5*(side1 + side2);	/* Calc area of quadrangle(trapezoid) */
			/* if src pixel center is included return larger area, otherwise return smaller area */
			return state.isIncludedSrcPixelCenter ? max(ret, 1 - ret) : min(ret, 1 - ret);	
		};
		auto type4 = [&]()->double {
			/* area of pentagon made by 1 line is 1 - triangle */
			return 1 - type2();
		};
		auto type5 = [&]()->double {
			/* area of pentagon made by 2 lines is calculated by the following equation
			 * type5 = 1 - (trapezoid + triangle)
			 */
			double shortBase, longBase, trapezoid = 0;
			double base, height, triangle = 0;
			if ( state.xCounts == 1 && state.yCounts == 3 ) {
				if ( state.intersections["xa"].size() == 0 ) {
					if ( state.intersections["ya"].size() == 1 ) {
						/* xa 0, xb 1, ya 1, yb 2
						 *
						 *   ya     yb
						 *  ─┼──────┼─xa
						 *   ･      │
						 *   │      :
						 *  ─┼──･───┼─xb
						 */
						shortBase = state.intersections["ya"][0];
						longBase = min(state.intersections["yb"][0], state.intersections["yb"][1]);
						base = 1 - state.intersections["xb"][0];
						height = 1 - max(state.intersections["yb"][0], state.intersections["yb"][1]);
					}
					else {
						/* xa 0, xb 1, ya 2, yb 1
						 *
						 *   ya     yb
						 *  ─┼──────┼─xa
						 *   │      ･
						 *   :      │
						 *  ─┼──･───┼─xb
						 */
						shortBase = min(state.intersections["ya"][0], state.intersections["ya"][1]);
						longBase = state.intersections["yb"][0];
						base = state.intersections["xb"][0];
						height = 1 - max(state.intersections["ya"][0], state.intersections["ya"][1]);
					}
				}
				else {
					if ( state.intersections["ya"].size() == 1 ) {
						/* xa 1, xb 0, ya 1, yb 2
						 *
						 *   ya     yb
						 *  ─┼──･───┼─xa
						 *   │      :
						 *   ･      │
						 *  ─┼──────┼─xb
						 */
						shortBase = 1 - state.intersections["ya"][0];
						longBase = 1 - max(state.intersections["yb"][0], state.intersections["yb"][1]);
						base = 1 - state.intersections["xa"][0];
						height = min(state.intersections["yb"][0], state.intersections["yb"][1]);
					}
					else {
						/* xa 1, xb 0, ya 2, yb 1
						 *
						 *   ya     yb
						 *  ─┼──･───┼─xa
						 *   :      │
						 *   │      ･
						 *  ─┼──────┼─xb
						 */
						shortBase = 1 - max(state.intersections["ya"][0], state.intersections["ya"][1]);
						longBase = 1 - state.intersections["yb"][0];
						base = state.intersections["xa"][0];
						height = min(state.intersections["ya"][0], state.intersections["ya"][1]);
					}
				}
			}
			else {
				if ( state.intersections["ya"].size() == 0 ) {
					if ( state.intersections["xa"].size() == 1 ) {
						/* xa 1, xb 2, ya 0, yb 1
						 *
						 *   ya     yb
						 *  ─┼･─────┼─xa
						 *   │      │
						 *   │      ･
						 *  ─┼──･･──┼─xb
						 */
						shortBase = state.intersections["xa"][0];
						longBase = min(state.intersections["xb"][0], state.intersections["xb"][1]);
						base = 1 - max(state.intersections["xb"][0], state.intersections["xb"][1]);
						height = 1 - state.intersections["yb"][0];
					}
					else {
						/* xa 2, xb 1, ya 0, yb 1
						 *
						 *   ya     yb
						 *  ─┼──･･──┼─xa
						 *   │      ･
						 *   │      │
						 *  ─┼･─────┼─xb
						 */
						shortBase = state.intersections["xb"][0];
						longBase = min(state.intersections["xa"][0], state.intersections["xa"][1]);
						base = 1 - max(state.intersections["xa"][0], state.intersections["xa"][1]);
						height = state.intersections["yb"][0];
					}
				}
				else {
					if ( state.intersections["xa"].size() == 1 ) {
						/* xa 1, xb 2, ya 1, yb 0
						 *
						 *   ya     yb
						 *  ─┼─────･┼─xa
						 *   │      │
						 *   ･      │
						 *  ─┼──･･──┼─xb
						 */
						shortBase = 1 - state.intersections["xa"][0];
						longBase = 1 - max(state.intersections["xb"][0], state.intersections["xb"][1]);
						base = min(state.intersections["xb"][0], state.intersections["xb"][1]);
						height = 1 - state.intersections["ya"][0];
					}
					else {
						/* xa 2, xb 1, ya 1, yb 0
						 *
						 *   ya     yb
						 *  ─┼──･･──┼─xa
						 *   ･      │
						 *   │      │
						 *  ─┼─────･┼─xb
						 */
						shortBase = 1 - state.intersections["xb"][0];
						longBase = 1 - max(state.intersections["xa"][0], state.intersections["xa"][1]);
						base = min(state.intersections["xa"][0], state.intersections["xa"][1]);
						height = state.intersections["ya"][0];
					}
				}
			}
			trapezoid = 0.5*(shortBase + longBase);
			triangle = 0.5*base*height;
			return 1 - trapezoid - triangle;
		};
		auto type6 = [&]()->double {
			/* area of hexagon is calculated by the following equation
			* type5 = 1 - (triangle1 + triangle2)
			*/
			double triangle1 = 0, triangle2 = 0;
			if ( state.intersections["xa"].size() == 2 ) {
				/* xa 2, xb 0, ya 1, yb 1
				 *
				 *   ya     yb
				 *  ─┼─･──･─┼─xa
				 *   ･      ･
				 *   │      │
				 *  ─┼──────┼─xb
				 */
				triangle1 = 0.5*min(state.intersections["xa"][0], state.intersections["xa"][1])*state.intersections["ya"][0];
				triangle2 = 0.5*(1 - max(state.intersections["xa"][0], state.intersections["xa"][1]))*state.intersections["yb"][0];
			}
			else if ( state.intersections["xb"].size() == 2 ) {
				/* xa 0, xb 2, ya 1, yb 1
				 *
				 *   ya     yb
				 *  ─┼──────┼─xa
				 *   │      │
				 *   ･      ･
				 *  ─┼─･──･─┼─xb
				 */
				triangle1 = 0.5*min(state.intersections["xb"][0], state.intersections["xb"][1])*(1 - state.intersections["ya"][0]);
				triangle2 = 0.5*(1 - max(state.intersections["xb"][0], state.intersections["xb"][1]))*(1 - state.intersections["yb"][0]);
			}
			else if ( state.intersections["ya"].size() == 2 ) {
				/* xa 1, xb 1, ya 2, yb 0
				 *
				 *   ya     yb
				 *  ─┼─･────┼─xa
				 *   ･      │
				 *   ･      │
				 *  ─┼─･────┼─xb
				 */
				triangle1 = 0.5*state.intersections["xa"][0] * min(state.intersections["ya"][0], state.intersections["ya"][1]);
				triangle2 = 0.5*state.intersections["xb"][0] * (1 - max(state.intersections["ya"][0], state.intersections["ya"][1]));
			}
			else if ( state.intersections["yb"].size() == 2 ) {
				/* xa 1, xb 1, ya 0, yb 2
				 *
				 *   ya     yb
				 *  ─┼────･─┼─xa
				 *   │      ･
				 *   │      ･
				 *  ─┼────･─┼─xb
				 */
				triangle1 = 0.5*(1 - state.intersections["xa"][0]) * min(state.intersections["yb"][0], state.intersections["yb"][1]);
				triangle2 = 0.5*(1 - state.intersections["xb"][0]) * (1 - max(state.intersections["yb"][0], state.intersections["yb"][1]));
			}
			return 1.0 - triangle1 - triangle2;
		};
		auto type7 = [&]()->double {
			double base, height;
			base = height = 0;
			for ( auto a : state.intersections ) {
				if ( a.second.size() == 2 ) {
					base = fabs(a.second[0] - a.second[1]);
					if ( a.first == "xa" )		height = state.vertexPos.second;
					else if ( a.first == "xb" )	height = 1 - state.vertexPos.second;
					else if ( a.first == "ya" )	height = state.vertexPos.first;
					else if ( a.first == "yb" )	height = 1 - state.vertexPos.first;
				}
			}
			return 0.5*base*height;
		};
		auto type8 = [&]()->double {
			double triangle1, triangle2;
			if ( state.intersections["xa"].size() == 1 && state.intersections["ya"].size() == 1 ) {
				/* xa 1, xb 0, ya 1, yb 0
				 *
				 *   ya     yb
				 *  ─┼─────･┼─xa
				 *   ･      │
				 *   │  ･   │
				 *  ─┼──────┼─xb
				 */
				triangle1 = 0.5*state.intersections["xa"][0] * state.vertexPos.second;
				triangle2 = 0.5*state.intersections["ya"][0] * state.vertexPos.first;
			}
			else if ( state.intersections["xa"].size() == 1 && state.intersections["yb"].size() == 1 ) {
				/* xa 1, xb 0, ya 0, yb 1
				 *
				 *   ya     yb
				 *  ─┼･─────┼─xa
				 *   │      ･
				 *   │   ･  │
				 *  ─┼──────┼─xb
				 */
				triangle1 = 0.5*(1 - state.intersections["xa"][0])*state.vertexPos.second;
				triangle2 = 0.5*state.intersections["yb"][0] * (1 - state.vertexPos.first);
			}
			else if ( state.intersections["xb"].size() == 1 && state.intersections["ya"].size() == 1 ) {
				/* xa 0, xb 1, ya 1, yb 0
				 *
				 *   ya     yb
				 *  ─┼──────┼─xa
				 *   │  ･   │
				 *   ･      │
				 *  ─┼─────･┼─xb
				 */
				triangle1 = 0.5*state.intersections["xb"][0] * (1 - state.vertexPos.second);
				triangle2 = 0.5*(1 - state.intersections["ya"][0])*state.vertexPos.first;
			}
			else {
				/* xa 0, xb 1, ya 0, yb 1
				 *
				 *   ya     yb
				 *  ─┼──────┼─xa
				 *   │   ･  │
				 *   │      ･
				 *  ─┼･─────┼─xb
				 */
				triangle1 = 0.5*(1 - state.intersections["xb"][0])*(1 - state.vertexPos.second);
				triangle2 = 0.5*(1 - state.intersections["yb"][0])*(1 - state.vertexPos.first);
			}
			return triangle1 + triangle2;
		};
		auto type9 = [&]() {
			double triangle1, triangle2, triangle3;
			if ( state.intersections["xa"].size() == 1 && state.intersections["xb"].size() == 1 ) {
				if ( max(state.intersections["xa"][0], state.intersections["xb"][0]) <= state.vertexPos.first ) {
					/* xa 1, xb 1, ya 0, yb 0
					 *
					 *   ya     yb
					 *  ─┼─･────┼─xa
					 *   │      │
					 *   │    ･ │
					 *  ─┼･─────┼─xb
					 */
					triangle1 = 0.5*state.intersections["xa"][0] * state.vertexPos.second;
					triangle2 = 0.5*state.vertexPos.first;
					triangle3 = 0.5*state.intersections["xb"][0] * (1 - state.vertexPos.second);
				}
				else {
					/* xa 1, xb 1, ya 0, yb 0
					 *
					 *   ya     yb
					 *  ─┼────･─┼─xa
					 *   │      │
					 *   │ ･    │
					 *  ─┼─────･┼─xb
					 */
					triangle1 = 0.5*(1 - state.intersections["xa"][0])*state.vertexPos.second;
					triangle2 = 0.5*(1 - state.vertexPos.first);
					triangle3 = 0.5*(1 - state.intersections["xb"][0])*(1 - state.vertexPos.second);
				}
			}
			else {
				if ( max(state.intersections["ya"][0], state.intersections["yb"][0]) <= state.vertexPos.second ) {
					/* xa 0, xb 0, ya 1, yb 1
					 *
					 *   ya     yb
					 *  ─┼──────┼─xa
					 *   .      ･
					 *   │     ･│
					 *  ─┼──────┼─xb
					 */
					triangle1 = 0.5*state.intersections["ya"][0] * state.vertexPos.first;
					triangle2 = 0.5*state.vertexPos.second;
					triangle3 = 0.5*state.intersections["yb"][0] * (1 - state.vertexPos.first);
				}
				else {
					/* xa 0, xb 0, ya 1, yb 1
					 *
					 *   ya     yb
					 *  ─┼──────┼─xa
					 *   │     ･│
					 *   .      ･
					 *  ─┼──────┼─xb
					 */
					triangle1 = 0.5*(1 - state.intersections["ya"][0])*state.vertexPos.first;
					triangle2 = 0.5*(1 - state.vertexPos.second);
					triangle3 = 0.5*(1 - state.intersections["yb"][0])*(1 - state.vertexPos.first);
				}
			}
			return triangle1 + triangle2 + triangle3;
		};

		if ( ! state.isIncludedDstPixelVertex ) {
			if ( state.xCounts == 0 && state.yCounts == 0 && ! state.isIncludedSrcPixelCenter )             { state.type = 0; return 0; }
			if ( state.xCounts == 0 && state.yCounts == 0 && state.isIncludedSrcPixelCenter )               { state.type = 1; return type1(); }
			if ( state.xCounts == 1 && state.yCounts == 1 && ! state.isIncludedSrcPixelCenter )             { state.type = 2; return type2(); }
			if ( (state.xCounts == 2 && state.yCounts == 0) || (state.xCounts == 0 && state.yCounts == 2) ) { state.type = 3; return type3(); }
			if ( state.xCounts == 1 && state.yCounts == 1 && state.isIncludedSrcPixelCenter )               { state.type = 4; return type4(); }
			if ( (state.xCounts == 3 && state.yCounts == 1) || (state.xCounts == 1 && state.yCounts == 3) ) { state.type = 5; return type5(); }
			if ( state.xCounts == 2 && state.yCounts == 2 )                                                 { state.type = 6; return type6(); }
			if ( state.xCounts == 0 && state.yCounts == 1 && ! state.isIncludedSrcPixelCenter )             { state.type = 0; return 0; }		/* boundary condition */
			if ( state.xCounts == 0 && state.yCounts == 1 && state.isIncludedSrcPixelCenter )               { state.type = 1; return type1(); }	/* boundary condition */
		}
		else {
			if ( (state.xCounts == 2 && state.yCounts == 0) || (state.xCounts == 0 && state.yCounts == 2) ) {
				for ( auto a : state.intersections ) {
					if ( a.second.size() == 2 ) {
						state.type = 7;
						return type7();
					}
				}
				state.type = 9;
				return type9();
			}
			if ( state.xCounts == 1 && state.yCounts == 1 ) {
				state.type = 8;
				return type8();
			}
		}
		return state.isIncludedSrcPixelCenter ? 1 : 0;	/* boundary condition */
	}
};

int main()
{
	/* lambda functions for reading & writing csv file */
	auto splitPath = [](string fullPath, string &path, string &base, string &extension) {
		/* split full path into drive, path, base and extension */
		size_t dotPos = fullPath.rfind(".");
		size_t delimiterPos = fullPath.rfind("\\");
		if (delimiterPos == string::npos) delimiterPos = fullPath.rfind("/");
		delimiterPos++;

		if (dotPos == string::npos) extension = "";
		else extension = fullPath.substr(dotPos);
		base = fullPath.substr(delimiterPos, dotPos - delimiterPos);
		path = fullPath.substr(0, delimiterPos);
	};
	auto split = [](string s, char delimiter)->vector<double> {
		/* split 'string' to 'double' elements for csv reading */
		size_t pos = s.find(delimiter);
		vector<double> ret;
		while ( pos != string::npos ) {
			try {
				double d = stod(s.substr(0, pos));
				ret.emplace_back(d);
			} catch ( const invalid_argument& ignored ) {
				// ignored
			}
			s = s.substr(pos + 1);
			pos = s.find(delimiter);
		}
		try {
			double d = stod(s.substr(0, pos));
			ret.emplace_back(d);
		} catch ( const invalid_argument& ignored ) {
			// ignored
		}
		return ret;
	};
	auto csvRead = [&](string path, IMG &data)->bool {
		ifstream fin(path);
		if ( ! fin ) {
			cout << "Failed to read csv file." << endl;
			return false;
		}

		string str;
		uiP imgSize = make_pair(0, 0);
		data.clear();
		while ( getline(fin, str) ) {
			data.resize(data.size() + 1);
			vector<double> vec = split(str, ',');
			if ( imgSize.first < vec.size() )	imgSize.first = vec.size();
			if ( vec.size() == 0 )	continue;
			for ( unsigned int i = 0; i < imgSize.first; ++i ) {
				data.back().emplace_back(vec[i]);
			}
			imgSize.second++;
		}
		return true;
	};
	auto csvWrite = [](string path, IMG data)->bool {
		ofstream fout(path);
		if ( ! fout ) {
			cout << "Failed to write csv file." << endl;
			return false;
		}
		if ( data.size() == 0 ) {
			cout << "There is no data in src array." << endl;
			cout << "Failed to write csv file." << endl;
			return false;
		}

		uiP imgSize = make_pair(data.front().size(), data.size());
		for ( unsigned int i = 0; i < imgSize.second; ++i ) {
			for ( unsigned int j = 0; j < imgSize.first; ++j ) {
				fout << data[i][j];
				if ( j < imgSize.first - 1 )	fout << ",";
			}
			fout << endl;
		}
		fout.close();
		return true;
	};

	
	/* Image parameters for interpolation
	 * (See 'areaAverageInterpolation' or 'fastAreaAverageInterpolation' function for details)
	 */
	pair<bool, string> ret;	
	IMG src, dst;
	dP srcResolution, dstResolution;
	dP srcIsocenter, dstIsocenter;
	double rotationAngle;
	int interpolationMode;

	/* User input parameters */
	string inputPath = "Test_film_dose.csv"; /* src image path (*.csv only) */
	srcResolution = make_pair(150, 150);     /* assumed x and y resolution are same */
	dstResolution = make_pair(25.4, 25.4);   /* assumed x and y resolution are same */
	srcIsocenter = make_pair(455, 455);      /* isocenter of src image */
	rotationAngle = 1.5;                     /* clockwise is positive */
	interpolationMode = 2;                   /* 1:Area average, 2:Fast area average */
	string path, base, extension;
	splitPath(inputPath, path, base, extension);
	if (extension != ".csv" && extension != ".CSV") {
		cout << "As for the image format, only csv format can be used." << endl;
		cout << "* path  : " << path << endl;
		cout << "* base  : " << base << endl;
		cout << "* ext   : " << extension << endl;
		cout << "Run terminated abnormally." << endl;
		return -1;
	}
	

	/* Reading source image */
	if ( ! csvRead(inputPath, src) ) {
		cout << "Run terminated abnormally." << endl;
		return -1;
	}


	/* Calculation of Area average or Fast area average interpolation
	 * - Fast area average interpolation is selected by default.
	 *   You can switch to Area average interpolation by manipulating the following comments.
	 */
	AreaAverageInterpolation aa;
	chrono::system_clock::time_point start, end; /* For calculation time measurement */
	start = chrono::system_clock::now();         /* For calculation time measurement */

	switch (interpolationMode) {
	case 1:
		/* Area average interpolation */
		ret = aa.areaAverageInterpolation(src, dst, srcResolution, dstResolution, srcIsocenter, dstIsocenter, rotationAngle);
		break;
	case 2:
		/* Fast area average interpolation */
		ret = aa.fastAreaAverageInterpolation(src, dst, srcResolution, dstResolution, srcIsocenter, dstIsocenter, rotationAngle);
		break;
	default:
		cout << "Invalid interpolation mode is selected." << endl;
		cout << "Interpolation mode should be 1 or 2." << endl;
		cout << " * Selected interpolation mode : " << interpolationMode << endl;
		cout << "Run terminated abnormally." << endl;
		return -1;
	}

	end = chrono::system_clock::now(); /* For calculation time measurement */
	double time = static_cast<double>(chrono::duration_cast<chrono::microseconds>(end - start).count() / 1000.0); /* For calculation time measurement */
	cout << "Calculation time : " << time << " [ms]" << endl; /* For calculation time measurement */

	if ( ! ret.first ) {
		cout << ret.second << endl;
		cout << "Run terminated abnormally." << endl;
		return -1;
	}


	/* Output interpolated image as csv */
	string outputPath = path + base + "_mod" + extension;
	if ( ! csvWrite(outputPath, dst) ) {
		cout << "Run terminated abnormally." << endl;
		return -1;
	}

	cout << "Run terminated correctly." << endl;
	return 0;
}