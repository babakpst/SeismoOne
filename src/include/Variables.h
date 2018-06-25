
extern  int NStep;            // number of time steps

extern  int LoadType;         // Load Type 0:prssure- 1: DRM

extern  int Wave_Type;        // Wave type (0: SV- 1: P) 
extern  int Wave_Func;        // Wave_Func (0: sine-1: Ricker)
extern  int Dis_History;      // Number of nodes for the history of displacement

extern  int NNodePWaveL;      // Number of Nodes per wavelength



// - Real variables -------------------------------------------------------------------------------------------------------------------------------
extern  double DT ;           // time step
extern  double Gama, Beta ;   // Newmark parameters
extern  double Alpha, P ;     // loading parameters
extern  double A ;            // area of the cross section of the beam
extern  double L;             // length of the beam
extern  double alpha1, alpha2;// Ricker pulse signal
extern  double amplitude ;    // Amplitude of the incident wave
extern  double omega ;        // central cyclic frequency in the Ricker pulse wave


 
extern  int * Nodal_History;  // a vector that holds the node numbers to record the history of displacement.
extern  int * Element_Layer;  // Number of elements in each layer
extern  int * Layer_Depth;    // Depth of each layer

// - Real arrays ----------------------------------------------------------------------------------------------------------------------------------
extern  double * F;           // global force vector

extern  double * Loc_History;    // a vector that holds the required locations to store the time history of displacement



