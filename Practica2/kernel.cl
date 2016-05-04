// Kernel code


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

  float val;  
  int out_rows = rows - 15;
  int out_cols = cols - 15;

  for(int i = 0; i < out_rows; i++)
    for(int j = 0; j < out_cols; j++)
    {
      val = getValue(img, cols, i, j);
      setValue(out, out_cols, i, j, val);
    }
}



