#include "ImageInterpolation.h"
#include "ColorSpaces.h"
#include <math.h>

#define M_PI 3.14151926535


void sampleAndHold(const uchar input[], int xSize, int ySize, uchar output[], int newXSize, int newYSize)
{
	uchar
		*yold = new uchar[xSize * ySize](),
		*ynew = new uchar[newXSize * newYSize]();

	char
		*vold = new char[xSize * ySize / 4](),
		*uold = new char[xSize * ySize / 4](),
		*vnew = new char[newXSize * newYSize / 4](),
		*unew = new char[newXSize * newYSize / 4]();

	RGBtoYUV420(input, xSize, ySize, yold, uold, vold);

	const double
		Sh = (double)newXSize / xSize,
		Sv = (double)newYSize / ySize;

	for (int i = 0; i < newYSize; ++i)
	{
		for (int j = 0; j < newXSize; ++j)
		{
			const int
				iprime = (i - 1) / Sv + 1,
				jprime = (j - 1) / Sh + 1;

			ynew[i * newXSize + j] = yold[iprime * xSize + jprime];
		}
	}

	for (int i = 0; i < newYSize / 2; ++i)
	{
		for (int j = 0; j < newXSize / 2; ++j)
		{
			const int
				iprime = (i - 1) / Sv + 1,
				jprime = (j - 1) / Sh + 1;

			vnew[i * newXSize / 2 + j] = vold[iprime * xSize / 2  + jprime];
			unew[i * newXSize / 2 + j] = uold[iprime * xSize / 2 + jprime];
		}
	}

	YUV420toRGB(ynew, unew, vnew, newXSize, newYSize, output);

	delete[] yold;
	delete[] ynew;
	delete[] vold;
	delete[] uold;
	delete[] vnew;
	delete[] unew;
}

void bilinearInterpolate(const uchar input[], int xSize, int ySize, uchar output[], int newXSize, int newYSize)
{
	uchar
		*yold = new uchar[xSize * ySize](),
		*ynew = new uchar[newXSize * newYSize]();

	char
		*vold = new char[xSize * ySize / 4](),
		*uold = new char[xSize * ySize / 4](),
		*vnew = new char[newXSize * newYSize / 4](),
		*unew = new char[newXSize * newYSize / 4]();

	RGBtoYUV420(input, xSize, ySize, yold, uold, vold);

	const double
		Sh = (double)newXSize / xSize,
		Sv = (double)newYSize / ySize;


	for (int i = 0; i < newYSize; ++i)
	{
		for (int j = 0; j < newXSize; ++j)
		{
			const double
				a = i / Sv - floor(i / Sv),
				b = j / Sh - floor(j / Sh);

			const int
				iprime = i / Sv,
				jprime = j / Sh;

			ynew[i * newXSize + j] =
				(1 - a) * (1 - b) * yold[iprime * xSize + jprime] +
				(1 - a) * b * yold[(iprime + 1) * xSize + jprime] +
				a * (1 - b) * yold[iprime * xSize + (jprime + 1)] +
				a * b * yold[(iprime + 1) * xSize + (jprime + 1)];
		}
	}

	for (int i = 0; i < newYSize / 2; ++i)
	{
		for (int j = 0; j < newXSize / 2; ++j)
		{
			const double
				a = i / Sv - floor(i / Sv),
				b = j / Sh - floor(j / Sh);

			const int
				iprime = i / Sv,
				jprime = j / Sh;

			unew[i * newXSize / 2 + j] =
				(1 - a) * (1 - b) * uold[iprime * xSize / 2 + jprime] +
				(1 - a) * b * uold[(iprime + 1) * xSize / 2 + jprime] +
				a * (1 - b) * uold[iprime * xSize / 2 + (jprime + 1)] +
				a * b * uold[(iprime + 1) * xSize / 2 + (jprime + 1)];

			vnew[i * newXSize / 2 + j] =
				(1 - a) * (1 - b) * vold[iprime * xSize / 2 + jprime] +
				(1 - a) * b * vold[(iprime + 1) * xSize / 2 + jprime] +
				a * (1 - b) * vold[iprime * xSize / 2 + (jprime + 1)] +
				a * b * vold[(iprime + 1) * xSize / 2 + (jprime + 1)];
		}
	}


	YUV420toRGB(ynew, unew, vnew, newXSize, newYSize, output);

	delete[] yold;
	delete[] ynew;
	delete[] vold;
	delete[] uold;
	delete[] vnew;
	delete[] unew;
}


void imageRotate(const uchar input[], int xSize, int ySize, uchar output[], int m, int n, double angle)
{
	uchar
		*yold = new uchar[xSize * ySize](),
		*ynew = new uchar[xSize * ySize]();

	char
		*vold = new char[xSize * ySize / 4](),
		*uold = new char[xSize * ySize / 4](),
		*vnew = new char[xSize * ySize / 4](),
		*unew = new char[xSize * ySize / 4]();

	const double θ = M_PI * angle / 180;


	RGBtoYUV420(input, xSize, ySize, yold, uold, vold);

	for (int y = 0; y < ySize; ++y)
	{
		for (int x = 0; x < xSize; ++x)
		{
			const int
				xprime = round(x * cos(θ) - y * sin(θ) - m * cos(θ) + n * sin(θ) + m),
				yprime = round(y * cos(θ) + x * sin(θ) - m * sin(θ) - n * cos(θ) + n);

			if (xprime >= xSize || xprime < 0 || yprime >= ySize || yprime < 0)
				ynew[y * xSize + x] = 0;
			else
				ynew[y * xSize + x] = yold[yprime * xSize + xprime];
		}
	}


	for (int y = 0; y < ySize / 2; ++y)
	{
		for (int x = 0; x < xSize / 2; ++x)
		{
			const int
				xprime = round(x  * cos(θ) - y  * sin(θ) - m / 2 * cos(θ) + n / 2 * sin(θ) + m / 2),
				yprime = round(y  * cos(θ) + x  * sin(θ) - m / 2 * sin(θ) - n / 2 * cos(θ) + n / 2);
			
			if (xprime >= xSize / 2 || xprime < 0 || yprime >= ySize / 2 || yprime < 0)
			{
				unew[y * xSize / 2 + x ] = 0;
				vnew[y * xSize / 2 + x ] = 0;
			}
			else
			{
				unew[y * xSize / 2 + x] = uold[yprime * xSize / 2 + xprime];
				vnew[y * xSize / 2 + x] = vold[yprime * xSize / 2 + xprime];
			}
		}
	}


	YUV420toRGB(ynew, unew, vnew, xSize, ySize, output);

	delete[] yold;
	delete[] ynew;
	delete[] vold;
	delete[] uold;
	delete[] vnew;
	delete[] unew;
}


void imageRotateBilinear(const uchar input[], int xSize, int ySize, uchar output[], int m, int n, double angle)
{
	uchar
		*yold = new uchar[xSize * ySize](),
		*ynew = new uchar[xSize * ySize]();

	char
		*vold = new char[xSize * ySize / 4](),
		*uold = new char[xSize * ySize / 4](),
		*vnew = new char[xSize * ySize / 4](),
		*unew = new char[xSize * ySize / 4]();

	const double θ = M_PI * angle / 180;


	RGBtoYUV420(input, xSize, ySize, yold, uold, vold);

	for (int y = 0; y < ySize; ++y)
	{
		for (int x = 0; x < xSize; ++x)
		{
			const double
				xprime = x * cos(θ) - y * sin(θ) - m * cos(θ) + n * sin(θ) + m,
				yprime = y * cos(θ) + x * sin(θ) - m * sin(θ) - n * cos(θ) + n;

			if (xprime >= xSize || xprime < 0 || yprime >= ySize || yprime < 0)
			{
				ynew[y * xSize + x] = 0;
			}
			else if (xprime != floor(xprime) || yprime != floor(yprime))
			{
				const int
					xbase = floor(xprime),
					ybase = floor(yprime);

				const double
					a = yprime - ybase,
					b = xprime - xbase;

				ynew[y * xSize + x] =
					(1 - a) * (1 - b) * yold[ybase * xSize + xbase] +
					(1 - a) * b * yold[(ybase + 1) * xSize + xbase] +
					a * (1 - b) * yold[ybase * xSize + (xbase + 1)] +
					a * b * yold[(ybase + 1) * xSize + (xbase + 1)];
			}
			else
			{
				ynew[y * xSize + x] = yold[(int)(yprime * xSize + xprime)];
			}
		}
	}


	for (int y = 0; y < ySize / 2; ++y)
	{
		for (int x = 0; x < xSize / 2; ++x)
		{
			const double
				xprime = x * cos(θ) - y * sin(θ) - m / 2 * cos(θ) + n / 2 * sin(θ) + m / 2,
				yprime = y * cos(θ) + x * sin(θ) - m / 2 * sin(θ) - n / 2 * cos(θ) + n / 2;

			if (xprime >= xSize / 2 || xprime < 0 || yprime >= ySize / 2 || yprime < 0)
			{
				unew[y * xSize / 2 + x] = 0;
				vnew[y * xSize / 2 + x] = 0;
			}
			else if (xprime != floor(xprime) || yprime != floor(yprime))
			{
				const int
					xbase = floor(xprime),
					ybase = floor(yprime);

				const double
					a = yprime - ybase,
					b = xprime - xbase;

				unew[y * xSize / 2 + x] =
					(1 - a) * (1 - b) * uold[ybase * xSize / 2 + xbase] +
					(1 - a) * b * uold[(ybase + 1) * xSize / 2 + xbase] +
					a * (1 - b) * uold[ybase * xSize / 2 + (xbase + 1)] +
					a * b * uold[(ybase + 1) * xSize / 2 + (xbase + 1)];

				vnew[y * xSize / 2 + x] =
					(1 - a) * (1 - b) * vold[ybase * xSize / 2 + xbase] +
					(1 - a) * b * vold[(ybase + 1) * xSize / 2 + xbase] +
					a * (1 - a) * vold[ybase * xSize / 2 + (xbase + 1)] +
					a * b * vold[(ybase + 1) * xSize / 2 + (xbase + 1)];
			}
			else
			{
				unew[y * xSize / 2 + x] = uold[(int)(yprime * xSize / 2 + xprime)];
				vnew[y * xSize / 2 + x] = vold[(int)(yprime * xSize / 2 + xprime)];
			}
		}
	}


	YUV420toRGB(ynew, unew, vnew, xSize, ySize, output);

	delete[] yold;
	delete[] ynew;
	delete[] vold;
	delete[] uold;
	delete[] vnew;
	delete[] unew;
}