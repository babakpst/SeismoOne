
#ifndef INPUTBASIC
#define INPUTBASIC

namespace main_ns
{

	namespace input_ns
	{

		class input_cls {

      private:
        string TempS;         // temporary variable for reading strings from input files

  		public:
        // variables
        string Name;          // name of the input file
        string Directory;     // Input/output directory

        string Input_Dir;            // Input directory
        string OutputMatlab_Dir;     // directory to write the input file for Matlab visualizer interface
        string Info_Dir;             // directory to write general information
        string FullFile_Dir;         // directory to write the full results in the time domain analysis
        string HistoryFile_Dir;      // directory to write the time history of displacement in the time domain analysis
        string TransferFunction_Dir; // directory to write the frequency domain results

        // functions
      	input_cls();
  			~input_cls();
		    void Address_fn();
	
		};

	}
}
#endif
