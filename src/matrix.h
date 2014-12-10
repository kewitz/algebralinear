#ifdef __MATRIX
#define __MATRIX

void SwapMLine(int, int, int, double*);
void SwapVLine(int, int, double*);

void LUCroutDecompose(int, double*, double*, double*);
void LUCroutInplaceDecompose(int, double*);

double GaussPivot(int, int, double*, double*);
void GaussJordan(int, double*, double*, double*);

#endif
