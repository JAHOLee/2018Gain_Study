#ifndef setup_config_h
#define setup_config_h

#include <string>
#include <map>
#include <fstream>
using namespace std;
const int MAXBOARDS = 94;
const int MAXSKI    = 4;
const int MAXCH     = 32;

class setup_config{
 public:
  setup_config();
  ~setup_config();

  //member
  string dirpath;
  int Module_List[MAXBOARDS];
  std::map<int,int> moduleID2BDorder;

  //function
  void Make_dir();  //Make output directories
  void Read_Module_List(string Module_Layout, int config); // Read config file (module->BD)
  
  
 private:
  
  bool DirectoryExists( const char* pzPath ); // Check if a directory exist
  
};

#endif
