int main(int argc, char **argv){
  // Optional arguments for infile, outfile: adapted from Analysis code for WCSim
  char * filename=NULL;
  char * outfilename=NULL;
  double NRG = 0;

  int verb = false;
  int c=-1;
  int nPMTtypes = 1;

  while( (c = getopt(argc,argv,"f:o:e:npmt:v")) != -1 ){
    // Input in c the argument (-f etc...) and in optarg the next argument.
    // When the above test becomes -1, it means it fails to find a new argument.
    switch(c){
    case 'f':
      filename = optarg;
      break;
    case 'o':
      outfilename = optarg;
      break;
    case 'e':
      NRG = atoi(optarg);
      break;
    case 'npmt':
      nPMTtypes = atoi(optarg);
    case 'v':
      verb = true;
      break;
    }
  }
  
  TApplication theApp("theApp", &argc, argv);

}
