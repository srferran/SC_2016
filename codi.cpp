/**
 *
 *  Implementation of the algorithm described in "Space-Time Completion of
 *  Video", by Y. Wexler, E. Shechtman, and M. Irani and published in IEEE
 *  Transactions on Pattern Analysis and Machine Intelligence, 2007.
 *
 *  Code written by L. Garrido, february 2014.
 *
 */

#include <iostream>
#include <limits>

#include "CImg.h"
#include <vector>
#include <time.h>
using namespace std;
using namespace cimg_library;

// Window size for Wp and Vq
#define SIZE_WINDOW   5

// Maximum number of iterations to perform at each level
#define MAX_ITER      10 

// A macro to check if an integer is even 
#define EVEN(x)       (!((x) & 0x01))

// start & end used to control computing time in diferent sections.
 clock_t start;
 clock_t end;
 clock_t total;

// Data structure used to store, for each window Wp, the
// best position for window Vq
class dataWp
{
  public:
    int x, y;
    double similarity;
};

/**
 *
 *  This function prints (an error) message to screen an exits the application. 
 *
 */

void sc_error(const char *str)
{
  cout << str << endl;
  exit(1);
}

/**
 *
 * Checks if all the pixels of the mask image have values 0 or 255. A value of
 * 255 identifies a pixel that belongs to the hole (i.e. not known) and a value
 * of 0 identifies a known pixel. This function is called after reading the
 * mask image from disc to ensure that the read data is correct. A value of 255
 * in the mask imate indicates that the pixel belongs to the hole, a value of 0
 * indicates that the pixel does not belong to the hole.
 *
 */

int sc_check_mask(
    CImg<unsigned char> &mask)
{
  unsigned long offset;

  // Check if the image is a gray level image
  if (mask.spectrum() != 1)
    return 1;

  // Check if the pixel values are 0 or 255
  for(offset = 0; offset < mask.size(); offset++)
  {
    if ((mask[offset] != 0) && (mask[offset] != 255))
      return 1;
  }

  // Check that there are no hole pixels near the borders of the image.  We
  // need to perform this checking since the algorithm, when accessing a pixel,
  // does not check if we are inside the image support. By checking this
  // condition we are sure that we won't access a pixel outside the image
  // support.

  int length = SIZE_WINDOW - 1;

  // Upper part of the image
  for(int y = 0; y < length; y++)
    for(int x = 0; x < mask.width(); x++)
      if (mask(x,y) == 255)
	return 1;

  // Bottom part of the image
  for(int y = mask.height() - length; y < mask.height(); y++)
    for(int x = 0; x < mask.width(); x++)
      if (mask(x,y) == 255)
	return 1;

  // Left part of the image
  for(int y = 0; y < mask.height(); y++)
    for(int x = 0; x < length; x++)
      if (mask(x,y) == 255)
	return 1;

  // Right part of the image
  for(int y = 0; y < mask.height(); y++)
    for(int x = mask.width() - length; x < mask.width(); x++)
      if (mask(x,y) == 255)
	return 1;

  // If we arrive here everything is ok for the mask image
  return 0;
}

/**
 *
 *  Given an input ori and a mask image (that has been read from disk), this
 *  function sets to zero all pixels of ori for which the mask image has pixel
 *  value 255 (that is, the hole). The other pixels (that do not belong to the
 *  hole) remain unchanged.
 *
 */

void sc_merge_ori_mask(
    CImg<unsigned char> &ori, 
    CImg<unsigned char> &mask, 
    CImg<unsigned char> &img)
{
  // Copy the contents of ori to img
  img = ori;

  // Search in mask for all pixels that have value 255 and set the
  // corresponding img value to zero.
  for(int y = 0; y < mask.height(); y++)
  {
    for(int x = 0; x < mask.width(); x++)
    {
      if (mask(x,y) == 255) {

	for(int c = 0; c < ori.spectrum(); c++)
	  img(x,y,c) = 0; 
      }
    }
  }
}

/**
 *
 * Given a mask image compute two new binary images: the image "mask_Wp" is a
 * binary image indicating the pixels for which we center the window Wp. For
 * each window Wp we search for the best window Vq at "mask_Vq". The binary
 * image "mask_Wp" contains all pixels such that the corresponding window Wp
 * contains at least one hole pixel. The image "mask_Vq" is a binary image
 * indicating the center of the windows Vq that are candidates for the windows
 * Wp. That is, in order to fill in the hole we search for each Wp centered at
 * "mask_Wp" the best window Vq centered at "mask_Vq" according to a color
 * similarity criterion. The image "mask_Vq" are all pixels of the image such
 * that it does not contain pixels that belong to the hole. 
 *
 */

void sc_get_mask_Vq_Wp(
    CImg<unsigned char> &mask,
    CImg<unsigned char> &mask_Vq,
    CImg<unsigned char> &mask_Wp)
{
  int length = (SIZE_WINDOW - 1) / 2;

  // We initialize the images with zero
  mask_Vq.fill(0);
  mask_Wp.fill(0);

  // Loop for all the pixels of the mask image
  for(int y = 0; y < mask.height(); y++)
  {
    for(int x = 0; x < mask.width(); x++)
    {
      int count_Vq  = 0;
      int count_Wp = 0;

      // For each pixel p of the mask image, we look at the pixels 
      // that are inside a window 5x5 centered at pixel p..
      for(int k = - length; k <= length; k++)
      {
	for(int l = - length; l <= length; l++)
	{
	  int pos_y = y + k;
	  int pos_x = x + l;

          // Check if the pixel (of the window) is inside the image
	  if ((pos_y >= 0) && (pos_y < mask.height()) && (pos_x >= 0) && (pos_x < mask.width())) {
	    if (mask(pos_x, pos_y) == 0)
	      count_Vq++;

	    if (mask(pos_x,pos_y) == 255)
	      count_Wp++;
	  }
	}
      }

      // A pixel (x,y) of the image "mask_Vq" is set to 255 if no pixels of the window
      // belong to the hole 
      if (count_Vq == (SIZE_WINDOW * SIZE_WINDOW))
	mask_Vq(x,y) = 255;

      // A pixel (x,y) of the image "mask_Wp" is set to 255 if at least a pixel of the
      // window contains a pixel that is associated to the hole.
      if (count_Wp > 0)
	mask_Wp(x,y) = 255;
    }
  }
}

/**
 *
 * Given a window Wp, this function searches for the corresponding best window
 * Vq. The similarity value is computed acording to the paper that describes
 * the method.
 *
 */

void sc_search_wp(
    int wp_x,
    int wp_y,
    CImg<unsigned char> &img,
    CImg<unsigned char> &mask_Vq,
    CImg<dataWp> &matrix)
{
  int length = (SIZE_WINDOW - 1) / 2;

  double best_dist = 1e300;

  // Full search through all the image
  for(int vq_y = 0; vq_y < img.height(); vq_y++)
  {
    for(int vq_x = 0; vq_x < img.width(); vq_x++)
    {
      // We just search for Vq for the pixels labelled
      // with 255 in the "mask_Vq" image
      if (mask_Vq(vq_x,vq_y) == 255)
      {
	double dist = 0.0;

        // We compare window Wp with Vq. For that reason
	// we run through all the pixels of the window.
	for(int k = -length; k <= length; k++)
	{
	  for(int l = -length; l <= length; l++)
	  {
	    // For all the color components
	    for(int c = 0; c < img.spectrum(); c++)
	    {
	      double value = (double) img(wp_x + l, wp_y + k, c) - (double) img(vq_x + l, vq_y + k, c);
	      dist += value * value;
	    }
	  }
	}

        // This is used to compute the best candidate Vq 
	if (dist < best_dist) {
	  best_dist = dist;
	  matrix(wp_x, wp_y).x = vq_x;
	  matrix(wp_x, wp_y).y = vq_y;
	}
      }
    }
  }

  // Fill data in matrix
  double sigma = 5 * SIZE_WINDOW * SIZE_WINDOW;
  matrix(wp_x,wp_y).similarity = exp( - best_dist / (2.0 * sigma * sigma));
}

/**
 *
 *  Applies the mean shift algorithm to a set of input colors stored in matrix
 *  candidates. Colors are grouped together according to its similarity. The
 *  number of candidate colors associated to each group is counted. The output
 *  matrix mask is a binary mask with values 1 and 0. A value of 1 indicates
 *  that the corresponding candidate color belongs to the group with most
 *  candidate colors. Otherwise the value is 0. 
 *
 */

 void sc_search_wp5(
    int wp_x,
    int wp_y,
    CImg<unsigned char> &img,
    CImg<unsigned char> &mask_Vq,
    CImg<dataWp> &matrix)
{
  int length = (SIZE_WINDOW - 1) / 2;

  double best_dist = 1e300;
  #pragma omp parallel
  {
    #pragma omp for schedule(dynamic) nowait
  	
  // Full search through all the image
	  for(int vq_y = 0; vq_y < img.height(); vq_y++)
	  {
	    for(int vq_x = 0; vq_x < img.width(); vq_x++)
	    {
      		// We just search for Vq for the pixels labelled
      		// with 255 in the "mask_Vq" image
	      if (mask_Vq(vq_x,vq_y) == 255)
	      {
			double dist = 0.0;

        	// We compare window Wp with Vq. For that reason
			// we run through all the pixels of the window.
			for(int k = -length; k <= length; k++)
			{
			  for(int l = -length; l <= length; l++)
			  {
			    // For all the color components
			    for(int c = 0; c < img.spectrum(); c++)
			    {
			      double value = (double) img(wp_x + l, wp_y + k, c) - (double) img(vq_x + l, vq_y + k, c);
			      dist += value * value;
			    }
			  }
			}

        // This is used to compute the best candidate Vq 
			if (dist < best_dist) 
			{
			  best_dist = dist;
			  matrix(wp_x, wp_y).x = vq_x;
			  matrix(wp_x, wp_y).y = vq_y;
			}
	      }
	    }
	  }
	
  }

  // Fill data in matrix
  double sigma = 5 * SIZE_WINDOW * SIZE_WINDOW;
  matrix(wp_x,wp_y).similarity = exp( - best_dist / (2.0 * sigma * sigma));
}
 
 
void sc_mean_shift(
	CImg<unsigned char> &candidates,
  	CImg<int> &mask)
{
  // For each color in candidates we indicate here its corresponding mean shift
  CImg<float> mean_shift(SIZE_WINDOW, SIZE_WINDOW, candidates.spectrum());

  // This matrix is used to indicate, for each mean_shift pixel, its corresponding cluster
  CImg<int> correspondence(SIZE_WINDOW, SIZE_WINDOW);

  // This vector is used to store the different clusters obtained in mean_shift
  CImg<float> clusters(SIZE_WINDOW * SIZE_WINDOW, candidates.spectrum());

  int length = (SIZE_WINDOW - 1) / 2;

  double value1, norm, sum;
  float *color, *new_color;

  color     = new float[candidates.spectrum()];
  new_color = new float[candidates.spectrum()]; 

  for(int k = -length; k <= length; k++) 
  {
    for(int l = -length; l <= length; l++)
    {

      // Initialization of the mean-shift
      for(int c = 0; c < candidates.spectrum(); c++)
	new_color[c] = (float) candidates(l + length, k + length, c);

      // Proceed to compute the mean-shift
      do {

	// Initialization of the colors 
	for(int c = 0; c < candidates.spectrum(); c++)
	{
	  color[c] = new_color[c];
	  new_color[c] = 0.0;
	}

	// Compute new_color according to mean-shift
	sum = 0.0;
	for(int k_shift = -length; k_shift <= length; k_shift++) 
	{
	  for(int l_shift = -length; l_shift <= length; l_shift++)
	  {
	    norm = 0.0;

	    // Compute the color norm difference
	    for(int c = 0; c < candidates.spectrum(); c++)
	    {
	      value1 = color[c] - (float) candidates(l_shift + length, k_shift + length, c);
	      norm += value1 * value1;
	    }

	    // Weight
	    value1 = exp(- norm / (2.0 * 5.0 * 5.0));
	    sum += value1;

	    // Apply it to the color (k_shift, l_shift)
	    for(int c = 0; c < candidates.spectrum(); c++)
	      new_color[c] += value1 * (float) candidates(l_shift + length, k_shift + length, c); 
	  }
	}

	// Normalize by sum
	for(int c = 0; c < candidates.spectrum(); c++)
	  new_color[c] /= sum;

	// Compute the norm between the new color and the previous one

	norm = 0.0;
	for(int c = 0; c < candidates.spectrum(); c++)
	{
	  value1 = new_color[c] - color[c];
	  norm += value1 * value1;
	}

	// Stop if the color has not changed a lot 

      } while (norm > 0.01);

      // Store the obtained color

      for(int c = 0; c < candidates.spectrum(); c++) 
	mean_shift(l + length, k + length, c) = new_color[c];
      
    }
  }

  // Free not necessary memory
  delete [] color;
  delete [] new_color;

  // Clusters initialization
  int *nitems = new int[SIZE_WINDOW * SIZE_WINDOW];
  int nclusters = 0;

  nitems[0] = 0;
  for(int c = 0; c < candidates.spectrum(); c++)  
    clusters(nclusters, c) = mean_shift(0, 0, c);
  nclusters++;

  // Scan through all the mean_shift image
  int cl;
  for(int k = -length; k <= length; k++) 
  {
    for(int l = -length; l <= length; l++)
    {
      // Search in which cluster does mean_shift(l + length, k + length, c) fit

      for(cl = 0; cl < nclusters; cl++)
      {
	norm = 0.0;
	for(int c = 0; c < candidates.spectrum(); c++)
	{
	  value1 = clusters(cl, c) - mean_shift(l + length, k + length, c);
	  norm += value1;
	}

	if (norm < 1.0)
	  break;
      }

      // Check if found
      if (cl < nclusters) {
	correspondence(l + length, k + length) = cl;
	nitems[cl]++;
      } else {
	// New cluster

	nitems[nclusters] = 1;
	correspondence(l + length, k + length) = nclusters;
	for(int c = 0; c < candidates.spectrum(); c++)  
	  clusters(nclusters, c) = mean_shift(l + length, k + length, c);

	nclusters++;
      }
    }
  }

  int pos_max = 0;
  value1 = nitems[0];

  // Search the maximum value and its position in vector nitems
  for(cl = 0; cl < nclusters; cl++)
  {
    if (nitems[cl] > value1) {
      value1 = nitems[cl];
      pos_max = cl;
    }
  }

  // Now mark all candidate pixels that belong to this cluster
  for(int k = -length; k <= length; k++) 
  {
    for(int l = -length; l <= length; l++)
    {
      if (correspondence(l + length, k + length) == pos_max)
	mask(l + length, k + length) = 1;
      else
	mask(l + length, k + length) = 0;
    }
  }

  delete [] nitems;
}

/**
 *
 * This function computes for a given pixel (that belongs to the hole) its
 * corresponding color value. The algorithm uses the windows Vq found for the
 * neighboring pixels as described in the paper.
 *
 */

void sc_get_color(
    int wp_x,
    int wp_y,
    CImg<unsigned char> &img,
    CImg<dataWp> &matrix,
    CImg<float> &dist)
{
  int length = (SIZE_WINDOW - 1) / 2;

  CImg<unsigned char> candidates(SIZE_WINDOW, SIZE_WINDOW, 1, img.spectrum());

  // We look for the neighboring pixels of (wp_x, wp_y)
  for(int k = -length; k <= length; k++) 
  {
    for(int l = -length; l <= length; l++)
    {
      // For the neighboring pixel we look at the 
      // corresponding Vq
      dataWp &data = matrix(wp_x + l, wp_y + k);

      // We store the color of each color in candidates
      for(int c = 0; c < img.spectrum(); c++)
	candidates(l + length, k + length, c) = img(data.x - l, data.y - k, c); 
    }
  }

  // This matrix is used to indicate, for each mean_shift pixel, its corresponding cluster
  CImg<int> mask(SIZE_WINDOW, SIZE_WINDOW);

  // Compute mean shift
  sc_mean_shift(candidates, mask);

  // We now compute the new color 
  double *color, sum_wi;

  // Color vector
  color = new double[img.spectrum()];

  sum_wi = 0.0;
  for(int c = 0; c < img.spectrum(); c++)
    color[c] = 0.0;

  // We look for the neighboring pixels of (wp_x, wp_y)
  for(int k = -length; k <= length; k++) 
  {
    for(int l = -length; l <= length; l++)
    {
      if (mask(l + length, k + length)) {
	// For the neighboring pixel we look at the 
	// corresponding Vq
	dataWp &data = matrix(wp_x + l, wp_y + k);

	// We compute the value of wi as explained in the
	// paper
	double alpha = pow(1.3, - dist(wp_x + l, wp_y + k));
	double wi = data.similarity * alpha; 

	// Here we compute the new color for pixel (wp_x, wp_y)
	sum_wi += wi;
	for(int c = 0; c < img.spectrum(); c++)
	  color[c] += wi * img(data.x - l, data.y - k, c); 
      }
    }
  }

  // Normalize color as indicated in the paper
  for(int c = 0; c < img.spectrum(); c++)
  {
    color[c] /= sum_wi;
    img(wp_x, wp_y, c) = (int) color[c];
  }

  delete [] color;
}

/**
 *
 * This is the function that allows to fill a hole for a given image.  Indeed,
 * the algorithm iteratively computes 1) the windows Vq that correspond to each
 * window Wq, 2) computes the color for each point of the hole using the
 * windows Vq.
 *
 */

void sc_fill_hole(
    CImg<unsigned char> &img,
    CImg<unsigned char> &mask)
	
{
	int iter;
	
  // Compute distance function within the mask
  CImg<float> dist = mask;
  dist.distance(0, 2);

  // Compute points at which we can perform the search
  CImg<unsigned char> mask_Vq(mask.width(), mask.height());
  CImg<unsigned char> mask_Wp(mask.width(), mask.height());

  sc_get_mask_Vq_Wp(mask, mask_Vq, mask_Wp);

  // Image at which search results are stored
  CImg<dataWp> matrix(img.width(), img.height()); 

  // This process is repeated several times. At each iteration
  // the color inside the hole is improved.
  iter = 0;
  while (iter < MAX_ITER)
  {
    cout << "  Iteration number " << iter+1 << endl;

    // Make a copy of the current image
    CImg<unsigned char> copy = img;

    // Search for all pixel holes
    for(int y = 0; y < img.height(); y++)
      for(int x = 0; x < img.width(); x++)
	if (mask_Wp(x,y) == 255)
	  //sc_search_wp(x, y, img, mask_Vq, matrix);
	  sc_search_wp5(x, y, img, mask_Vq, matrix);  //< --------------

    // Once we have searched for all windows Wp in "mask_Vq", 
    // we fill the hole with the new color
    for(int y = 0; y < img.height(); y++)
      for(int x = 0; x < img.width(); x++)
	if (mask(x,y) == 255)
	  sc_get_color(x, y, img, matrix, dist);

    // Check if any pixel has changed. If no pixel
    // has changed, we can exit the iterations.
    if (copy == img)
      break; 

    iter++;
  }
}


/**
 *  Class is needed at sc_fill_hole_2 to save coordinates
 */
	
class Coordinates
{
	public:
		int x,y;
};
/**
 *
 * Same method as sc_fill_hole but parallel computing improvement
 *
 */

void sc_fill_hole2(
    CImg<unsigned char> &img,
    CImg<unsigned char> &mask)
	
{
	int iter;
	vector<Coordinates> vector_xy ;

  // Compute distance function within the mask
  CImg<float> dist = mask;
  dist.distance(0, 2);

  // Compute points at which we can perform the search
  CImg<unsigned char> mask_Vq(mask.width(), mask.height());
  CImg<unsigned char> mask_Wp(mask.width(), mask.height());

  sc_get_mask_Vq_Wp(mask, mask_Vq, mask_Wp);

  // Image at which search results are stored
  CImg<dataWp> matrix(img.width(), img.height()); 

  // This process is repeated several times. At each iteration
  // the color inside the hole is improved.
	
  iter = 0;
	
    // <----------   Added  P1  ----------->
    // save hole pixels coordinates

    for(int y = 0; y < img.height(); y++)		
    {
        for(int x = 0; x < img.width(); x++)
        {
            if (mask_Wp(x,y) == 255)
            {
                Coordinates c;
                c.x = x;
                c.y = y;
                vector_xy.push_back(c);
            }
        }
    }

   // <---------------- end -------------->
  while (iter < 10)
  {
    cout << "  Iteration number " << iter+1 << endl;

    // Make a copy of the current image
    CImg<unsigned char> copy = img;

    /*<--------  P1 Added section --------->*/
    
    #pragma omp parallel num_threads(6)
    {
        #pragma omp for nowait
        for (unsigned int j =0; j<vector_xy.size(); j++)
        {     
            sc_search_wp(vector_xy[j].x, vector_xy[j].y, img, mask_Vq, matrix);
        }
    }
	
    /*------------- End section ----------------*/
    // Once we have searched for all windows Wp in "mask_Vq", 
    // we fill the hole with the new color

    for(int y = 0; y < img.height(); y++)
      for(int x = 0; x < img.width(); x++)
        if (mask(x,y) == 255)
            sc_get_color(x, y, img, matrix, dist);



		

    // Check if any pixel has changed. If no pixel
    // has changed, we can exit the iterations.
    if (copy == img)
      break; 
    iter++;
  }
	
}
void sc_fill_hole3(
    CImg<unsigned char> &img,
    CImg<unsigned char> &mask)
	
{
	int iter;
	
  // Compute distance function within the mask
  CImg<float> dist = mask;
  dist.distance(0, 2);

  // Compute points at which we can perform the search
  CImg<unsigned char> mask_Vq(mask.width(), mask.height());
  CImg<unsigned char> mask_Wp(mask.width(), mask.height());

  sc_get_mask_Vq_Wp(mask, mask_Vq, mask_Wp);

  // Image at which search results are stored
  CImg<dataWp> matrix(img.width(), img.height()); 

  // This process is repeated several times. At each iteration
  // the color inside the hole is improved.
  iter = 0;
  while (iter < MAX_ITER)
  {
    cout << "  Iteration number " << iter+1 << endl;

    // Make a copy of the current image
    CImg<unsigned char> copy = img;
   /*<--------  P1 Added section --------->*/                              
    #pragma omp parallel  num_threads(6)
   {
      // Search for all pixel holes
      
      #pragma omp single
      {
        for(int y = 0; y < img.height(); y++)
        {
          for(int x = 0; x < img.width(); x++)
          {
            if (mask_Wp(x,y) == 255)
            {
              #pragma omp task
              {
                sc_search_wp(x, y, img, mask_Vq, matrix);
              }
            }
          }
        }
      }
    }
    // <---------------- end -------------->            

    // Once we have searched for all windows Wp in "mask_Vq", 
    // we fill the hole with the new color
    for(int y = 0; y < img.height(); y++)
      for(int x = 0; x < img.width(); x++)
	if (mask(x,y) == 255)
	  sc_get_color(x, y, img, matrix, dist);

    // Check if any pixel has changed. If no pixel
    // has changed, we can exit the iterations.
    if (copy == img)
      break; 

    iter++;
  }
}


void sc_fill_hole4(
    CImg<unsigned char> &img,
    CImg<unsigned char> &mask)
	
{
	int iter;
	vector<Coordinates> vector_xy ;

  // Compute distance function within the mask
  CImg<float> dist = mask;
  dist.distance(0, 2);

  // Compute points at which we can perform the search
  CImg<unsigned char> mask_Vq(mask.width(), mask.height());
  CImg<unsigned char> mask_Wp(mask.width(), mask.height());

  sc_get_mask_Vq_Wp(mask, mask_Vq, mask_Wp);

  // Image at which search results are stored
  CImg<dataWp> matrix(img.width(), img.height()); 

  // This process is repeated several times. At each iteration
  // the color inside the hole is improved.
	
  iter = 0;
	
    // <----------   Added  P1  ----------->
    // save hole pixels coordinates

    for(int y = 0; y < img.height(); y++)		
    {
        for(int x = 0; x < img.width(); x++)
        {
            if (mask_Wp(x,y) == 255)
            {
                Coordinates c;
                c.x = x;
                c.y = y;
                vector_xy.push_back(c);
            }
        }
    }

   // <---------------- end -------------->
  while (iter < 10)
  {
    cout << "  Iteration number " << iter+1 << endl;

    // Make a copy of the current image
    CImg<unsigned char> copy = img;

    /*<--------  P1 Added section --------->*/
    
    #pragma omp parallel num_threads(6)
    {
        #pragma omp single
        for (unsigned int j =0; j<vector_xy.size(); j++)
        {   
	    #pragma omp task shared(img,matrix)
            sc_search_wp(vector_xy[j].x, vector_xy[j].y, img, mask_Vq, matrix);
        }
    }
	
    /*------------- End section ----------------*/
    // Once we have searched for all windows Wp in "mask_Vq", 
    // we fill the hole with the new color

    for(int y = 0; y < img.height(); y++)
      for(int x = 0; x < img.width(); x++)
        if (mask(x,y) == 255)
            sc_get_color(x, y, img, matrix, dist);



		

    // Check if any pixel has changed. If no pixel
    // has changed, we can exit the iterations.
    if (copy == img)
      break; 
    iter++;
  }
	
}

/**
 *
 *  Binarizes the mask setting its values to 0 or 255. This function is used
 *  when constructing the multiscale images of the mask.
 *
 */

void sc_binarize_mask(
    CImg<unsigned char> &mask)
{
  for(int y = 0; y < mask.height(); y++)
  {
    for(int x = 0; x < mask.width(); x++)
    {
      if (mask(x,y) < 128)
	mask(x,y) = 0;
      else
	mask(x,y) = 255;
    }
  }
}

/**
 *
 *  Given a mask this function constructs a reduced version of the image. For
 *  this issue the mask image is first filtered (blurred) out and then one of
 *  each four pixels are taken to reduce the image in size.
 *
 */

void sc_downscale_mask(
    CImg<unsigned char> &mask,
    CImg<double> &kernel,
    CImg<unsigned char> &mask_scaled)
{
  mask.convolve(kernel);

  for(int y = 0; y < mask_scaled.height(); y++)
    for(int x = 0; x < mask_scaled.width(); x++)
      mask_scaled(x,y) = mask(2*x, 2*y);
}

/**
 *
 *  Given an imge this function constructs a reduced version of the image. For
 *  this issue the image is first filtered (blurrred) out and then one of each
 *  four pixels are taken to reduce the image size. It should be pointed out
 *  that filtering is performed such that pixel that fall inside the whole are
 *  not taken into account.
 *
 */

void sc_downscale_img(
    CImg<unsigned char> &img,
    CImg<double> &kernel,
    CImg<unsigned char> &mask,
    CImg<unsigned char> &mask_scaled,
    CImg<unsigned char> &img_scaled)
{
  // We implement the filtering for the image: we only take those pixels that
  // do no belong to the mask.

  float *value = new float[img.spectrum()];

  // We run through all the pixels of the image
  for(int y = 0; y < img.height(); y++)
  {
    for(int x = 0; x < img.width(); x++)
    {
      // Check the value of the mask
      if (mask(x, y) == 255) {

	for(int c = 0; c < img.spectrum(); c++)
	  img(x, y, c) = 0;

      } else {

	// If the mask value belongs to a known pixel, we use
	// the neighboring known pixels to filter

	int m = 0;
	float sum_w = 0.0;

	for(int c = 0; c < img.spectrum(); c++) 
	  value[c] = 0;

	// Look for neighboring pixels
	for(int k = -1; k <= 1; k++)
	{
	  for(int l = -1; l <= 1; l++, m++)
	  {
	    int pos_y = y + k;
	    int pos_x = x + l;

	    // Check if the pixel falls inside the image support and if it is a known pixel
	    if ((pos_y >= 0) && (pos_y < mask.height()) && (pos_x >= 0) && (pos_x < mask.width())) 
	      if (mask(pos_x, pos_y) == 0)
	      {
		sum_w += kernel[m];

		for(int c = 0; c < img.spectrum(); c++)
		  value[c] += kernel[m] * img(pos_x, pos_y, c);
	      }
	  }
	}

	// Final value
	for(int c = 0; c < img.spectrum(); c++)
	  img(x, y, c) = value[c] / sum_w;
      }
    }
  }

  delete [] value;

  // Scale the image 
  for(int y = 0; y < mask_scaled.height(); y++)
    for(int x = 0; x < mask_scaled.width(); x++)
    {
      if (mask_scaled(x,y) == 0)
      {
	// If mask is not hole take the pixel value

	for(int c = 0; c < img.spectrum(); c++)
	  img_scaled(x, y, c) = img(2*x, 2*y, c);
      }
      else
      {
	// If mask is hole put image value to 0

	for(int c = 0; c < img.spectrum(); c++)
	  img_scaled(x, y, c) = 0.0;
      }
    }
}

/**
 *
 *  Reduce the size of the image and mask by two.
 *
 */

int sc_scale_images(
    CImg<unsigned char> &img,
    CImg<unsigned char> &mask,
    CImg<unsigned char> &img_scaled,
    CImg<unsigned char> &mask_scaled)
{
  double filter1d[3] = {0.25, 0.5, 0.25};

  int m = 0;
  CImg<double> kernel(3,3);
  for(int i = 0; i < 3; i++)
    for(int j = 0; j < 3; j++, m++)
      kernel[m] = filter1d[i] * filter1d[j];

  // Reduce the size of the mask by 2
  sc_downscale_mask(mask, kernel, mask_scaled);

  // Binarize the mask (only 0 and 255 values)
  sc_binarize_mask(mask_scaled);

  // If the mask has no hole pixel or the hole pixel is too near to the border
  // return a 0 so that the multiscale construction is stopped
  if ((mask_scaled == 0) || (sc_check_mask(mask_scaled)))
    return 0;

  // If the mask contains at least one hole pixel scale the image
  sc_downscale_img(img, kernel, mask, mask_scaled, img_scaled);

  // Return 1 to indicate that the mask pixel contains at least one hole pixel.
  return 1;
}

/**
 *
 *  Construct a multiscale representation of the image. That is, construct
 *  reduced versions of the image and mask image to make processing easier. The
 *  size of the image is reduced until the mask has no pixels. Each of the
 *  images of the multiscale representation is stored in a list and afterwards
 *  used for the hole filling.
 *
 */

void sc_multiscale(
    CImg<unsigned char> img,
    CImg<unsigned char> mask,
    CImgList<unsigned char> &multiscale_img,
    CImgList<unsigned char> &multiscale_mask)
{
  int nrow, ncol, new_nrow, new_ncol;

  int pos = 0;
  bool end = false;

  multiscale_img.insert(img, pos);
  multiscale_mask.insert(mask, pos);

  ncol = img.width();
  nrow = img.height();

  // Repeat the reduction of the images and mask until
  // the mask does not contain any hole pixel.
  while (!end)
  {
    // Compute the new sizes of the images

    if (EVEN(ncol))
      new_ncol = ncol >> 1;
    else
      new_ncol = (ncol + 1) >> 1;

    if (EVEN(nrow))
      new_nrow = nrow >> 1;
    else
      new_nrow = (nrow + 1) >> 1;

    CImg<unsigned char> img_scaled(new_ncol, new_nrow, 1, img.spectrum()); 
    CImg<unsigned char> mask_scaled(new_ncol, new_nrow);

    // Scale the image and mask

    if (sc_scale_images(img, mask, img_scaled, mask_scaled))
    {
      // If the mask image contains at least a hole pixel
      // insert the image and the mask image into the 
      // list.

      pos++;

      multiscale_img.insert(img_scaled, pos);
      multiscale_mask.insert(mask_scaled, pos);

      ncol = new_ncol;
      nrow = new_nrow;

      img  = img_scaled;
      mask = mask_scaled;
    }
    else
    {
      end = true;
    }
  }
}

/**
 *
 *  Assume an image with a hole and its correspoding reduced version on which
 *  the hole has been filled. This function initializes the image using the
 *  information of the reduced image.
 *
 */

void sc_upscale_img(
    CImg<unsigned char> &img_scaled,
    CImg<unsigned char> &mask,
    CImg<unsigned char> &img)
{
  // Take the reduced image and scale it to the current
  // size of the image
  CImg<unsigned char> img_new = img_scaled;
  img_new.resize(img.width(), img.height(), 1, img.spectrum(), 3);

  // For the current image just initialize those pixels
  // that are inside the hole.
  for(int y = 0; y < img.height(); y++)
    for(int x = 0; x < img.width(); x++)
      if (mask(x,y) == 255)
      {
	for(int c = 0; c < img.spectrum(); c++)
	  img(x,y,c) = img_new(x,y,c); 
      }
}

/**
 *
 *  Main function. Reads the image and the mask from disc and fills in the hole
 *  using a multiscale approach.
 *
 */

int main(int argc, char **argv)
{
	
  if (argc != 3) {
    cout << argv[0] << " image mask" << endl;
    exit(1);
  }

  // Read input images from disc
  CImg<unsigned char> ori(argv[1]);
  CImg<unsigned char> mask(argv[2]);

  // Check that everthing if ok 
  if ((ori.width() != mask.width()) || (ori.height() != mask.height())) 
    sc_error("Input image size differ."); 

  if (sc_check_mask(mask))
    sc_error("Input mask should be a one channel image with 0 or 255 values.");

  // Set pixels of the image that are within the mask to zero.
  // Thus, we are sure we do not cheat. Instead of working with
  // ori, we will work with img.
  CImg<unsigned char> img;
  sc_merge_ori_mask(ori, mask, img);

  // List of images to store the multiscale images
  CImgList<unsigned char> multiscale_img;
  CImgList<unsigned char> multiscale_mask;

  // Construct the multiscale images and mask

  sc_multiscale(img, mask, multiscale_img, multiscale_mask);

  // We go through all the multiscale images starting
  // at the most reduced image.

  for(int i = multiscale_img.size() - 1; i >= 0; i--)
  {
    cout << "Level " << i << endl;

    // Fill the hole of the image
    start = clock();// <------------------------------------ time control, START -----------------------------> START
    sc_fill_hole2(multiscale_img(i), multiscale_mask(i));
    end = clock();// <------------------------------------ time control, END -----------------------------> END
    // Upscale the image to the next level
    if (i > 0)
      sc_upscale_img(multiscale_img(i), multiscale_mask(i-1), multiscale_img(i-1));
  }
  
  total += end - start;// <------------------------------------ time control, ACCUMULATE TOTAL -----------------------------> ACCUMULATE TOTAL
  printf("Time : %f\n",(double)total/CLOCKS_PER_SEC );
  
  // Result...
  CImg<unsigned char> output = multiscale_img(0);
  output.save("output.png");
		
  return 0;
	
	
	
}
