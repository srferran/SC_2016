
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

float calculateValue(int i_pos, int j_pos, int cols, __global unsigned char *img, __global unsigned char *pat)
{
  float val_img, val_pat;
  float result = 0;
  for (int i = 0; i < 16; i++)
  {
    for (int j = 0; j < 16 ; j++)
    {
      val_img = getValue(img, cols, i_pos + i, j_pos + j);
      val_pat = getValue(pat, 16, i , j);
      result += (val_img - val_pat)*(val_img - val_pat);     
    }
  }
  // (1.0/256.0)= 0.00390625;
  result *= 0.00390625;
  return result;
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
      //val = getValue(img, cols, i, j);
      val = calculateValue(i, j, cols, img, pat);
      setValue(out, out_cols, i, j, val);
    }
}





