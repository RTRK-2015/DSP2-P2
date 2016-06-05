
#include "ImageProcessing.h"
#include "ImageInterpolation.h"

#include <QDebug>


void imageProcessingFun(const QString& progName, QImage* const outImgs, const QImage* const inImgs, const QVector<double>& params)
{
  const int X_SIZE = inImgs->width();
  const int Y_SIZE = inImgs->height();
  /* NOTE: Calculate output image resolution and construct output image object */

  auto upto = [](int x, int y)
  {
	return (x + y - 1) & ~(y - 1);
  };


  if (progName == "Sample and hold")
  {
    /* Input image data in RGB format can be obtained with inImgs->bits() */
    /* Vertical scale factor is params[0] */
    /* Horizontal scale factor is params[1] */

    /* TO DO: Calculate output image resolution and construct output image object */
	  const int
		X_NEWSIZE = upto(X_SIZE * params[1], 4),
		Y_NEWSIZE = upto(Y_SIZE * params[0], 4);

    new (outImgs) QImage(X_NEWSIZE, Y_NEWSIZE, inImgs->format());

    /* TO DO: Perform Sample and hold interpolation  */
	sampleAndHold(
		inImgs->bits(),
		X_SIZE,
		Y_SIZE,
		outImgs->bits(),
		X_NEWSIZE,
		Y_NEWSIZE);
  }
  else if (progName == "Bilinear")
  {
    /* Input image data in RGB format can be obtained with inImgs->bits() */
    /* Vertical scale factor is params[0] */
    /* Horizontal scale factor is params[1] */

    /* TO DO: Calculate output image resolution and construct output image object */
	const int
		X_NEWSIZE = upto(X_SIZE * params[1], 4),
		Y_NEWSIZE = upto(Y_SIZE * params[0], 4);

	new (outImgs) QImage(X_NEWSIZE, Y_NEWSIZE, inImgs->format());

	/* TO DO: Perform Bilinear interpolation  */
	bilinearInterpolate(
		inImgs->bits(),
		X_SIZE,
		Y_SIZE,
		outImgs->bits(),
		X_NEWSIZE,
		Y_NEWSIZE);
  }
  else if (progName == "Rotation")
  {
    /* Input image data in RGB format can be obtained with inImgs->bits() */
    /* Rotation angle in degrees is params[0]*/
    /* Center of rotation coordinates are (XSIZE/2, YSIZE/2) */
    /* TO DO: Construct output image object */
    new (outImgs) QImage(X_SIZE, Y_SIZE, inImgs->format());

	imageRotate(
		inImgs->bits(),
		X_SIZE,
		Y_SIZE,
		outImgs->bits(),
		X_SIZE / 2,
		Y_SIZE / 2,
		params[0]);
  }
  else if (progName == "Rotation Bilinear")
  {
    /* Input image data in RGB format can be obtained with inImgs->bits() */
    /* Rotation angle in degrees is params[0]*/
    /* Center of rotation coordinates are (XSIZE/2, YSIZE/2) */

    /* TO DO: Construct output image object */
    new (outImgs) QImage(X_SIZE, Y_SIZE, inImgs->format());

    /* TO DO: Perform image rotation with bilinear interpolation */
	imageRotateBilinear(
		inImgs->bits(),
		X_SIZE,
		Y_SIZE,
		outImgs->bits(),
		X_SIZE / 2,
		Y_SIZE / 2,
		params[0]);
  }

}

