
#include <iostream>
#include "CImg.h"

using namespace std;
using namespace cimg_library;

CImg<float> pattern_matching(CImg<unsigned char> &img, CImg<unsigned char> &pat);

CImg<unsigned char> setBorder(CImg<unsigned char> &gray)
{
  CImg<unsigned char> gray_border(gray.width()+15,gray.height()+15,1,1,0);

  unsigned char *grayPtr = gray._data;
  unsigned char *bordPtr = gray_border._data;

  for(int i = 0; i < gray.height(); i++)
    for(int j = 0; j < gray.width(); j++)
    {
      int offset_gray = i * gray.width() + j;
      int offset_bord = i * gray_border.width() + j;

      bordPtr[offset_bord] = grayPtr[offset_gray];
    }

  return gray_border;
}

int main( int argc, char** argv )
{
  if( argc != 4)
  {
    cout <<" Usage: " << argv[0] << " <image_in> <pattern_in> <image_out>" << endl;
    return -1;
  }

  CImg<unsigned char> image(argv[1]);   
  CImg<unsigned char> pattern(argv[2]);

  if ((pattern.width() != 16) || (pattern.height() != 16)) {
    cout << "Pattern has not size 16 x 16" << endl;
    exit(1);
  }

  CImg<unsigned char> gray;

  int color = image.spectrum(); 
  switch(color)
  {
    case 1:
      gray = image;
      break;
    case 3: 
      gray = image.get_RGBtoYCbCr().get_channel(0);
      break;
    default:
      cout << "Error: image is not graylevel nor color image" << endl;
  }

  CImg<unsigned char> gray_border = setBorder(gray);
  CImg<float> out = pattern_matching(gray_border, pattern);
  CImg<unsigned char> norm_out= out.get_normalize(0, 255);

  norm_out.save(argv[3]);

  return 0;
}

void search()
{
  for (int I_x  = 0; I_x<image.width; I_x++)
  {
    for (int I_y  = 0; I_y<image.height; I_y++)
    {
        for(int i = 0; i<16; i++)
        {
          for(int j = 0; j<16; j++)
          {
             float acumulat += pow((img(I_x + j,I_y + i) - pattern(j,i)),2);
          }
        }
        E(I_x,I_y)= (1/256)*acumulat;
        acumulat = 0;
    }
  }

}
