// Kernel code

#define BLOCK_SIZE_H 32
#define BLOCK_SIZE_W 32

#define pattern_H 16
#define pattern_W 16


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
  float val= 0.0;  
  int out_rows = rows - 15;
  int out_cols = cols - 15;
  int row = get_global_id(1);
  int col = get_global_id(0);
  int paso = rows/BLOCK_SIZE_W; //rows = img_width
*/
  // k y e no pueden ser 0 pues todos los hilos empezarian por cero
  //for(int k = row; k < 256 ; k++)
  //{
  int row = get_global_id(1);
  int col = get_global_id(0);

  float out_value = 0;
  for (int e = 0; e < rows; e++)
    out_value = img[row * rows + e] ;
  out[row * rows + col] = out_value;
 
      

      
  //} 
   // (1.0/256.0)= 0.00390625;
    //out[row*out_rows + col] = val ;
   //val = getValue(img, cols, i, j);
   //setValue(out, out_cols, i, j, val);
    
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
