// Kernel code
include 
#define BLOCK_SIZE 32



/* cambiar width y heigth
*/
unsigned char getValue(__global unsigned char *img, int cols, int i, int j)
{
  float val = img[i * cols + j]; 
  return val;
}

void setValue(__global float *out, int cols, int i, int j, float value)
{
  out[i * cols + j] = value; 
}

__kernel void pattern_matching(
    __global unsigned char *img,
    __global unsigned char *pat,
    int rows,
    int cols,
    __global float *out)
{

/*
  // pattern length = pattern width * heigth = 16*16
  __local float pattern_local[256];
  
  // copy from global to local memory
  for (int i = 0 ; i < 16 ; i++)
  {
    for (int j = 0); j < 16 ; j++)
    {
      pattern_local[i * cols + j] = pat[i][j];
    }
  }
*/
  float val= 0.0;  
  int out_rows = rows - 15;
  int out_cols = cols - 15;
  int row = get_global_id(1);
  int col = get_global_id(0);
  int paso = img.width()/BLOCK_SIZE;
  for(int k = 0; k < paso ; k++)
  {
      for (int e = 0; e < pat.width() ; e++)
      {
          val +=  img[row*paso + e] - pat[e * pat.width() + col] * img[row*paso + e] - pat[e * pat.width() + col];
      }
  } 
    out[row*out.width + col] = val ;
   //val = getValue(img, cols, i, j);
   //setValue(out, out_cols, i, j, val);
   val = 0.0; 
}

/*

  float acumulat = 0.0;
  for (int I_y  = 0; I_y<(image.height()-15); I_y++)
  {
    for (int I_x  = 0; I_x<(image.width()-15); I_x++)
    {
        for(int j = 0; j<16; j++)
        {
          for(int i = 0; i<16; i++)
          {

             acumulat += pow((image(I_y + j,I_x + i) - pattern(j,i)),2);
          }
        }
        // (1.0/256.0)= 0.00390625;
        E(I_y,I_x)= (0.00390625)*acumulat;
        acumulat = 0.0;
    }
  }
  return E;
*/
