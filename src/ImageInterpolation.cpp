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

	for (int i = 0; i < newYSize / 4; ++i)
	{
		for (int j = 0; j < newXSize; ++j)
		{
			const int
				iprime = (i - 1) / Sv + 1,
				jprime = (j - 1) / Sh + 1;

			vnew[i * newXSize + j] = vold[iprime * xSize + jprime];
			unew[i * newXSize + j] = uold[iprime * xSize + jprime];
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
				b = j / Sh - floor(j / Sh),
				a = i / Sv - floor(i / Sv);

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

	for (int i = 0; i < newYSize / 4; ++i)
	{
		for (int j = 0; j < newXSize; ++j)
		{
			const double
				a = j / Sh - floor(j / Sh),
				b = i / Sv - floor(i / Sv);

			const int
				iprime = i / Sv,
				jprime = j / Sh;

			unew[i * newXSize + j] =
				(1 - a) * (1 - b) * uold[iprime * xSize + jprime] +
				(1 - a) * b * uold[(iprime + 1) * xSize + jprime] +
				a * (1 - b) * uold[iprime * xSize + (jprime + 1)] +
				a * b * uold[(iprime + 1) * xSize + (jprime + 1)];

			vnew[i * newXSize + j] =
				(1 - a) * (1 - b) * vold[iprime * xSize + jprime] +
				(1 - a) * b * vold[(iprime + 1) * xSize + jprime] +
				a * (1 - b) * vold[iprime * xSize + (jprime + 1)] +
				a * b * vold[(iprime + 1) * xSize + (jprime + 1)];
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
				xprime = x * cos(θ) - y * sin(θ) - m * cos(θ) + n * sin(θ) + m,
				yprime = y * cos(θ) + x * sin(θ) - m * sin(θ) - n * cos(θ) + n;

			if (xprime >= xSize || xprime < 0 || yprime >= ySize || yprime < 0)
				ynew[y * xSize + x] = 0;
			else
				ynew[y * xSize + x] = yold[yprime * xSize + xprime];
		}
	}


	/* No idea
	for (int y = 0; y < ySize / 4; ++y)
	{
		for (int x = 0; x < xSize; ++x)
		{
			const int
				xprime = x * cos(θ) - y * sin(θ) - m * cos(θ) + n * sin(θ) + m,
				yprime = y * cos(θ) + x * sin(θ) - m * sin(θ) - n * cos(θ) + n;

			if (xprime >= xSize || xprime < 0 || yprime >= ySize / 4 || yprime < 0)
			{
				unew[y * xSize + x] = 0;
				vnew[y * xSize + x] = 0;
			}
			else
			{
				unew[y * xSize + x] = uold[yprime * xSize + xprime];
				vnew[y * xSize + x] = vold[yprime * xSize + xprime];
			}
		}
	}*/


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


	/* Still no idea
	for (int y = 0; y < ySize / 4; ++y)
	{
		for (int x = 0; x < xSize; ++x)
		{
			const double
				xprime = x * cos(θ) - y * sin(θ) - m * cos(θ) + n * sin(θ) + m,
				yprime = y * cos(θ) + x * sin(θ) - m * sin(θ) - n * cos(θ) + n;

			if (xprime >= xSize || xprime < 0 || yprime >= ySize / 4 || yprime < 0)
			{
				unew[y * xSize + x] = 0;
				vnew[y * xSize + x] = 0;
			}
			else if (xprime != floor(xprime) || yprime != floor(yprime))
			{
				const int
					xbase = floor(xprime),
					ybase = floor(yprime);

				const double
					a = yprime - ybase,
					b = xprime - xbase;

				unew[y * xSize + x] =
					(1 - a) * (1 - b) * uold[ybase * xSize + xbase] +
					(1 - a) * b * uold[(ybase + 1) * xSize + xbase] +
					a * (1 - b) * uold[ybase * xSize + (xbase + 1)] +
					a * b * uold[(ybase + 1) * xSize + (xbase + 1)];

				vnew[y * xSize + x] =
					(1 - a) * (1 - b) * vold[ybase * xSize + xbase] +
					(1 - a) * b * vold[(ybase + 1) * xSize + xbase] +
					a * (1 - a) * vold[ybase * xSize + (xbase + 1)] +
					a * b * vold[(ybase + 1) * xSize + (xbase + 1)];
			}
			else
			{
				unew[y * xSize + x] = uold[(int)(yprime * xSize + xprime)];
				vnew[y * xSize + x] = vold[(int)(yprime * xSize + xprime)];
			}
		}
	}*/


	YUV420toRGB(ynew, unew, vnew, xSize, ySize, output);

	delete[] yold;
	delete[] ynew;
	delete[] vold;
	delete[] uold;
	delete[] vnew;
	delete[] unew;
}