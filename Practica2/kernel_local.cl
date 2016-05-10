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

float calculateValue(int i_pos, int j_pos, int cols, __local unsigned char **img, __local unsigned char **pat)
{
  float val_img, val_pat;
  float result = 0;
  for (int i = 0; i < 16; i++)
  {
    for (int j = 0; j < 16 ; j++)
    {
      val_img = img[i][j];
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

  int row = get_global_id(0);
  int col = get_global_id(1);

   
  int out_rows = rows - 15;
  int out_cols = cols - 15;

  float val;

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
  for(int i=0;i<rows/BLOCK_SIZE; i++)
  {
  	for(int j=0; j<rows/BLOCK_SIZE; j++)
    {
    	for(int k=0; k<16; k++)
    	{
    		for(int l=0; l<2; l++)
	    	{

	    		img_local[8 * k + get_local_id(0)][32 * l + get_local_id(1)] = 
	    		img[((rows/BLOCK_SIZE)*32 + get_global_id(0)) * 256 + (rows/BLOCK_SIZE)*32 + get_global_id(1)];

	    	}
    	}
      barrier(CLK_LOCAL_MEM_FENCE);
      val = calculateValue(i, j, cols, img_local, pat_local);
      barrier(CLK_LOCAL_MEM_FENCE);
      out[((i * BLOCK_SIZE + get_local_id(1))* 5) + (j * 32) + get_local_id(0) ] = val;
      //setValue(out, out_cols, i, j, val);
      barrier(CLK_LOCAL_MEM_FENCE);
    }
  }
}



