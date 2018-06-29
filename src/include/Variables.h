

extern  int LoadType;         // Load Type 0:prssure- 1: DRM


extern  int Dis_History;      // Number of nodes for the history of displacement

extern  int NNodePWaveL;      // Number of Nodes per wavelength


extern  int * Nodal_History;  // a vector that holds the node numbers to record the history of displacement.
extern  int * Element_Layer;  // Number of elements in each layer
extern  int * Layer_Depth;    // Depth of each layer

// - Real arrays ----------------------------------------------------------------------------------------------------------------------------------
extern  double * F;           // global force vector

extern  double * Loc_History;    // a vector that holds the required locations to store the time history of displacement



