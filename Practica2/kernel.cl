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



