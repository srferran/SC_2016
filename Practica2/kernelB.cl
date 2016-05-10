
// Kernel code
#define BLOCK_SIZE 32
unsigned char getValue(__global unsigned char *img, int cols, int i, int j)
{
  float val = img[i * cols + j]; 
  return val;
}

void setValue(__global float *out, int cols, int i, int j, float value)
{
  out[i * cols + j] = value; 
}

float calculateValue(int i_pos, int j_pos, int cols, __global unsigned char *img, __local unsigned char **pat)
{
  float val_img, val_pat;
  float result = 0;
  for (int i = 0; i < 16; i++)
  {
    for (int j = 0; j < 16 ; j++)
    {
      val_img = getValue(img, cols, i_pos + i, j_pos + j);
      val_pat =pat[i][j];
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
  __local unsigned char pat_local[16][16];
  __local unsigned char img_local[48][64];
  // copiem el pattern a memoria local 
  // es un sol work grup de 8x32 work items (16*16)
  int pos0 = get_local_id(1)/16 + get_local_id(0) * 2;
  int pos1 = get_local_id(1)%16;

  pat_local[pos0][pos1] = pat[get_local_id(0) * 32 + get_local_id(1)];
  barrier(CLK_LOCAL_MEM_FENCE);

  //copiem la imatge a la memoria local en bocins de 64x64 per 
  //processar 32x32 work items, es així perque per processar 
  // p.e el pixel (32,32) ens cal tenir el valor dels pixels
  // fins al (32+15,32+15) i ja que hi som, com els workgrups
  // son de 8*32 -> 32+32 = 64 rows i [8+8+8+8]32  [+8+8]16 de border = 48 cols

  

  int row = get_global_id(0);
  int col = get_global_id(1);
  float val; 
  int out_rows = rows - 15;
  int out_cols = cols - 15;

  for(int i = row; i < out_rows; i++)
    for(int j = col; j < out_cols; j++)
    { 
      //val = getValue(img, cols, i, j);
      val = calculateValue(i, j, cols, img, pat_local);

      setValue(out, out_cols, i, j, val);
    }
}





