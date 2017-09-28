/*-----------------------------------------------------------------*/
/*  Desc:     A main function to read and place                    */
/*                                                                 */
/*  Authors:  Hsieh, Yuan-Ting, NTUEE (b01901046@ntu.edu.tw)       */
/*            Jin Hu, IBM Corporation (jinhu@us.ibm.com)           */
/*                                                                 */
/*  Created:  06/10/2015                                           */
/*-----------------------------------------------------------------*/

#include "evaluate.h"
#define _DEBUG

int main(int argc, char** argv)
{

  
  if(argc != 4 &&argc !=5 ) // can accept four argc
  {
    cout << "Incorrect arguments. exiting .."<<endl;
    cout << "Usage : myplacer ICCAD15.parm [.iccad2015] [target_util] (optional)[outputfile name]" << endl ;
    return 0;
  }
	cout << "Command line : " << endl;
	for(int i=0 ; i<argc ; i++)
		cout << " " << argv[i];
	cout << endl;

  circuit ckt;
  
  ckt.read_parameters(argv[1]);
  ckt.read_iccad2015_file(argv[2]);

  ckt.buildTraverseOrder();

  //ckt.copy_init_to_final();

//for parameter
  //double tt=atof(argv[5]);
  //cout<<endl<<endl<<endl<<endl<<endl<<endl<<tt<<endl<<endl<<endl<<endl<<endl<<endl<<endl;
  //ckt.setA1(tt);

#ifdef DEBUG
  ckt.traverseCells();
#endif

  ckt.initWeight();
  ckt.initHPWL();

	
  //int maxIter=5;// really need to be 5? or 3?
  int maxIter=1;
  ckt.copy_init_to_final();
  ckt.buildSiteMapping();
  while (maxIter)
  {
     ckt.staticTimingAnalysis();
     ckt.netWeightUpdate();
     ckt.discreteLocalSearch();
     ckt.buildSiteMapping();
     ckt.buildSizeTree();
     ckt.buildDistTreesQQQ();
     ckt.removeOverlapping();
     --maxIter;
  }
#ifdef DEBUG
  ckt.print();
#endif

  
  if (argc==5) 
	  {
      ofstream os(argv[4]);
      ckt.outputDEF(os);
    }
  else
    {
      ofstream os("test.def");
      ckt.outputDEF(os);
    }
	

  cout << endl;
  cout << "Evaluating the solution file ..                                                 " <<endl;
  cout << "--------------------------------------------------------------------------------" <<endl;

  cout << "Analyzing placement .. "<<endl;
  ckt.measure_displacement();
  cout << "  max displ. (um) : " << ckt.displacement << endl;
  ckt.measure_ABU(BIN_DIM, atof(argv[3]));
  cout << "  ABU penalty     : " << ckt.ABU_penalty <<endl;
  cout << "  alpha           : " << ALPHA  <<endl;
  cout << endl;

  cout << "Analyzing timing .. "<<endl;
  ckt.measure_HPWL();
  bool timer_done = ckt.measure_timing();
  cout << "  HPWL, StWL (um) : " << ckt.total_HPWL << ", " <<ckt.total_StWL << endl;
  cout << "  Scaled StWL     : " << ckt.total_StWL * (1 + ALPHA * ckt.ABU_penalty) << " ( "<< ALPHA * ckt.ABU_penalty*100 << "% )"<<endl;
  cout << "  Clock period    : " << ckt.clock_period << endl;
  if(timer_done)
  {
    cout << "  early WNS, TNS  : " << ckt.eWNS << ", " << ckt.eTNS <<endl;
    cout << "  late  WNS, TNS  : " << ckt.lWNS << ", " << ckt.lTNS <<endl;
  }
  else
    cout << "  WNS, TNS        : Timer failed. The values are not available." <<endl;

  cout << "--------------------------------------------------------------------------------" <<endl;
  cout << endl;
  cout << "Checking legality .. " <<endl;
  cout << "--------------------------------------------------------------------------------" <<endl;
  if(ckt.check_legality())
    cout << "Placement is LEGAL."<<endl;
  else
    cout << "Placement is ILLEGAL. see check_legality.log."<<endl;
  cout << "--------------------------------------------------------------------------------" <<endl;
	

  return 1;
}
