
#ifndef INPUTBASIC
#define INPUTBASIC

void InputBasic ( int &  NEl_DRM, int & Solver, int & Dis_History, double & alpha1, double & alpha2, double & amplitude, double & omega, int& Wave_Type, int& Wave_Func, int& NInt, int& NDim,int& NNode, int& NDOF, int& NMat, int& NPM, int& NNodePWaveL, int& LoadType, double& Beta, double& Gama, double& DT, int& NStep, double& L, double& Alpha, double& P, double& A, int& OShFunc, ifstream& InputFile);
void InputArrays(double & omega, int & NNodePWaveL, int & Dis_History, int& NMat, double*& Length, double& L, double*& Loc_History, int*& Element_Layer, int*& Layer_Depth, int & NEl, int& OShFunc, int & NJ,  double **& PMat, ifstream& InputFile);

void InputMatlab ( int & NJ, int & Dis_History, int& NMat, int& NEl, double& DT, int& NStep, double& L, ofstream& OutputMatlab, int& NEqM, int& NPM, double*& Length, double**& PMat, double*& Loc_History, double**& XYZ);


#endif
