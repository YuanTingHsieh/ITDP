#include "evaluate.h"
#include <iostream>
#include <list>
#define _DEBUG

using namespace std;

void circuit::discreteLocalSearch()
{
	//below needs to be reverse topological order!!!!!
	for (list<unsigned>::iterator it = traverseList.begin(); it != traverseList.end(); it++)
	{
		int i = *it;
		if (!(cells[i].isFixed))
		   localSearch(i);
		else
			{
				cells[i].x_coord=cells[i].init_x_coord;
				cells[i].y_coord=cells[i].init_y_coord;
			}
	}
}

void circuit::localSearch(int cellid)
{
	double best_cost=numeric_limits<double>::max();
	double best_sit=cells[cellid].init_x_coord;
	//int total_windows=8;//actually only seven pos is searched
	int total_windows=33;
	int A1=1;
	//int A2=1;
	//int beta=1;
	///below are weight for late and early cost
	double alpha=5.0;
	double omega=1.0;

	for (unsigned int i=1;i<total_windows;++i)
	{
		
		///may have prob!? discrete
		int newcoord=cells[cellid].init_x_coord+i*rows[0].stepX;
		if (newcoord+cells[cellid].width>rows[0].origX+rows[0].numSites*rows[0].stepX)
			break;

		cells[cellid].x_coord=newcoord;
		
		//HPWL for output interconnect
		double netMaxX, netMinX;
        double netMaxY, netMinY;
        for (map<string, unsigned>::iterator it=cells[cellid].ports.begin() ; it!=cells[cellid].ports.end() ; it++)
	    {
	    	if (pins[(*it).second].portType==1)
	    	{
	    		net* theNet=&nets[pins[(*it).second].net];//one net can have two fan outs
      	 	   	netMaxX=netMinX=pins[ (*it).second ].x_offset+cells[cellid].x_coord;
     	    	netMaxY=netMinY=pins[ (*it).second ].y_offset+cells[cellid].y_coord;
       			for(vector<unsigned>::iterator theSink=theNet->sinks.begin() ; theSink != theNet->sinks.end() ; ++theSink)
      			{
          			netMaxX=max(netMaxX, pins[ *theSink ].x_coord);
          			netMinX=min(netMinX, pins[ *theSink ].x_coord);
          			netMaxY=max(netMaxY, pins[ *theSink ].y_coord);
          			netMinY=min(netMinY, pins[ *theSink ].y_coord);
        		}
        	}
        	else
        		netMaxX=netMinX=netMaxY=netMinY=0;
  		}
  		//cout<<"where do you crash? "<<endl;
  		double new_HPWL=(netMaxX-netMinX+netMaxY-netMinY) / static_cast<double>(DEFdist2Microns);
  		double delta_CAP=(new_HPWL-cells[cellid].HPWL)*LOCAL_WIRE_CAP_PER_MICRON*10000000;

  		//cout<<"The original HPWL is "<<cells[cellid].HPWL<<" new_HPWL is "<<new_HPWL<<endl;
  		double cost_Late=0.0;
  		double cost_Early=0.0;
  		for (map<string, unsigned>::iterator it=cells[cellid].ports.begin() ; it!=cells[cellid].ports.end() ; it++)
  		{
  			//for fan in
  			pin* sinknow=&pins[(*it).second];
  			if(sinknow->portType==0)
  			{				
		        for (unsigned int j=0;j<nets[sinknow->net].sinks.size();++j)
		        {
		   	       	if (nets[sinknow->net].sinks[j]==(sinknow)->id)
		   	  		{
		   	  			// I think the two below is the same!!
  						double delay_Early=sinknow->earlyAt-pins[nets[sinknow->net].source].earlyAt;
  						double delay_Late=sinknow->lateAt-pins[nets[sinknow->net].source].lateAt;
  						cost_Late+=(delay_Late+A1*delta_CAP)*nets[sinknow->net].lateWeights[j];
  						cost_Early+=(delay_Early+A1*delta_CAP)*nets[sinknow->net].earlyWeights[j];
  						#ifdef DEBUG
  						 // cout<<"delay_Late is "<<delay_Late<<" , A1*delta_CAP is "<<A1*delta_CAP<<endl;
  						 // cout<<"delay_Early is "<<delay_Early<<endl;
  						#endif
  					}
  				}
  			}
  			//for fan out
  			else if (sinknow->portType==1)
  			{
  				for (unsigned int j=0;j<nets[sinknow->net].sinks.size();++j)
		        {
		        	double MANHAT=0.0;
		        	if (pins[nets[sinknow->net].sinks[j]].type!=2)//not PO
		        	{
		        		MANHAT=(abs(cells[cellid].x_coord-cells[pins[nets[sinknow->net].sinks[j]].owner].x_coord)
		        	              +abs(cells[cellid].y_coord-cells[pins[nets[sinknow->net].sinks[j]].owner].y_coord))
		        	              / static_cast<double>(DEFdist2Microns);
  					}
  					else
  						MANHAT=(abs(cells[cellid].x_coord-pins[nets[sinknow->net].sinks[j]].x_coord)
		        	              +abs(cells[cellid].y_coord-pins[nets[sinknow->net].sinks[j]].y_coord))
		        	              / static_cast<double>(DEFdist2Microns);
  					double tow=pins[nets[sinknow->net].sinks[j]].lateAt-sinknow->lateAt;
  					double delta_TOW=0.69*LOCAL_WIRE_RES_PER_MICRON*MANHAT
  					*(LOCAL_WIRE_CAP_PER_MICRON*new_HPWL/2+nets[sinknow->net].lumped_wire_cap)-tow;
  					cost_Late+=(tow+delta_TOW)*nets[sinknow->net].lateWeights[j];
  					cost_Early+=(tow+delta_TOW)*nets[sinknow->net].earlyWeights[j];
				}
  			}
  		}
  		double cost_now=alpha*cost_Late+omega*cost_Early;
  		#ifdef DEBUG
  		  // cout<<"LOCAL searching... cell "<<cells[cellid].name<<" windows "<<i<<"current location is "<<cells[cellid].x_coord<<endl;
  		  // cout<<"The cost is "<<cost_now<<endl<<endl<<endl;
  		#endif
  		if (cost_now<best_cost)
  		{
  			best_cost=cost_now;
  			best_sit=cells[cellid].x_coord;
  		}
	}

	cells[cellid].x_coord=best_sit;
	cells[cellid].y_coord=cells[cellid].init_y_coord;
}

void circuit::initHPWL()
{
	update_pinlocs();
	for (unsigned int cellid=0;cellid<cells.size();++cellid)
	{
		//HPWL for output interconnect
		double netMaxX, netMinX;
        double netMaxY, netMinY;
        for (map<string, unsigned>::iterator it=cells[cellid].ports.begin() ; it!=cells[cellid].ports.end() ; it++)
	    {
	    	if (pins[(*it).second].portType==1)//may have two outputs
	    	{
	    		pin* thePin=&pins[(*it).second];
	    		net* theNet=&nets[pins[(*it).second].net];//one net can have two fan outs
      	 	   	netMaxX=netMinX=thePin->x_coord;
     	    	netMaxY=netMinY=thePin->y_coord;
     	    	//thePin->print();
       			for(vector<unsigned>::iterator theSink=theNet->sinks.begin() ; theSink != theNet->sinks.end() ; ++theSink)
      			{
          			netMaxX=max(netMaxX, pins[ *theSink ].x_coord);
          			netMinX=min(netMinX, pins[ *theSink ].x_coord);
          			netMaxY=max(netMaxY, pins[ *theSink ].y_coord);
          			netMinY=min(netMinY, pins[ *theSink ].y_coord);
          			//cout<<"Changing (netMinX, netMaxX, netMinY, netMaxY) to "<<netMinX<<" , "<<netMaxX<<" , "<<netMinY<<" , "<<netMaxY<<endl;
        		}
        	}
  		}
  		cells[cellid].HPWL=(netMaxX-netMinX+netMaxY-netMinY) / static_cast<double>(DEFdist2Microns);
  		#ifdef DEBUG
  		//cout<<"Initialing cell "<<cells[cellid].name<<" HPWL is "<<cells[cellid].HPWL<<endl;
  		#endif
	} 
}

bool circuit::staticTimingAnalysis()
{
  update_pinlocs();
  build_steiner();
  slice_longwires(MAX_WIRE_SEGMENT_IN_MICRON * static_cast<double>(DEFdist2Microns));

	// 1. write .spef, .tau2015, .timing, .ops to run a timer
	string filename = design_name;
	write_tau2015(filename+".tau2015");
	write_spef(filename+".spef");
	write_timing(filename+".timing");
  
  // 2. run UI-Timer
	char* arguments[5];
	for(int i=0 ; i<5 ; i++)
		arguments[i]=new char[80];
	strcpy(arguments[0], "ui-timer2.0");
	strcpy(arguments[1], (filename+".tau2015").c_str());
	strcpy(arguments[2], (filename+".timing").c_str());

	using namespace uit;
	Timer timer;
	timer.init_iccad2015_timer(5, arguments);
	timer.set_num_threads(8);
  
	// we assume PIs are driven by PI drivers from .sdc file during timer invocation
  for(vector<unsigned>::iterator PI=PIs.begin() ; PI!=PIs.end() ; ++PI)
	{
		if(pins[ *PI ].name == clock_port)
			continue;
		timer.disconnect_pin(pins[*PI].name);
		macro *theDriver = &macros[ pins[ *PI ].driverType ];
		string driverName = pins[ *PI ].name+"_drv";
		timer.insert_instance(driverName, theDriver->name);
		timer.insert_net(driverName);
		timer.connect_pin(pins[ *PI ].name, driverName);
		// NOTE: assume PI drivers are inverters
		timer.connect_pin(driverName+":a", driverName); 
		timer.connect_pin(driverName+":o", pins[ *PI ].name);
	}

	// 3. read out results
	eWNS=eTNS=0.0;
	lWNS=lTNS=0.0;
	for(vector<pin>::iterator thePin=pins.begin() ; thePin!=pins.end() ; ++thePin)
	{
		thePin->earlySlk=min(timer.report_slack(thePin->name, EARLY, RISE), timer.report_slack(thePin->name, EARLY, FALL)); 
		thePin->lateSlk =min(timer.report_slack(thePin->name, LATE,  RISE), timer.report_slack(thePin->name, LATE,  FALL));
		thePin->earlyAt=max(timer.report_at(thePin->name, EARLY, RISE), timer.report_at(thePin->name, EARLY, FALL)); 
		thePin->lateAt =min(timer.report_at(thePin->name, LATE,  RISE), timer.report_at(thePin->name, LATE,  FALL));
		thePin->earlyRat=max(timer.report_rat(thePin->name, EARLY, RISE), timer.report_rat(thePin->name, EARLY, FALL)); 
		thePin->lateRat =min(timer.report_rat(thePin->name, LATE,  RISE), timer.report_rat(thePin->name, LATE,  FALL)); 

		// collect slacks @ timing end points
		if(thePin->type == PO_PIN || (thePin->isFlopInput && !thePin->isFlopCkPort))
		{
#ifdef DEBUG
			cout << thePin->name << " " << thePin->earlySlk << " " << thePin->lateSlk <<endl;
//this is HSIEH modify
			cout<<"FUCKFUCK"<<endl;
			cout<< thePin->name<<" "<<thePin->delay<<endl;
#endif
			eWNS = min(eWNS, thePin->earlySlk);
			lWNS = min(lWNS, thePin->lateSlk);
			eTNS += min(0.0, thePin->earlySlk);
			lTNS += min(0.0, thePin->lateSlk);
		}
		//cout << thePin->name << " " << thePin->earlySlk << " " << thePin->lateSlk <<endl;
	}

  return true;
}

void circuit::initWeight()
{
	for (unsigned int i=0;i<nets.size();++i)
	{
		//double sourceEarlySlk=pins[nets[i].source].earlySlk;
		//double sourceLateSlk=pins[nets[i].source].lateSlk;
		for (unsigned int j=0;j<nets[i].sinks.size();++j)
		{
			nets[i].earlyWeights.push_back(1.0);
			nets[i].lateWeights.push_back(1.0);
		}
	}


	//below needs to be reverse topological order!!!!!
	for (list<unsigned>::iterator it = traverseList.begin(); it != traverseList.end(); it++)
		netWeightDistribute(*it);
}

void circuit::netWeightUpdate()
{
	//update PO
	for (unsigned int i=0;i<POs.size();++i)
	{
		pin* POPO=&pins[POs[i]];
		POPO->lateWeight=POPO->lateWeight*(POPO->lateAt)/(POPO->lateRat);
		POPO->earlyWeight=POPO->earlyWeight*(POPO->earlyRat)/(POPO->earlyAt);
	}

	//below needs to be reverse topological order!!!!!
	for (list<unsigned>::iterator it = traverseList.begin(); it != traverseList.end(); it++)
		netWeightBalance(*it);
	/////////////
	for (list<unsigned>::iterator it = traverseList.begin(); it != traverseList.end(); it++)
		netWeightDistribute(*it);
    /////////////////

}

void circuit::netWeightBalance(int cellid)
{
	if (cells[cellid].isLCB)// ignoring LCBs
	{
#ifdef DEBUG
	//cout<<"LLLLCCCCBBB"<<endl;
#endif
	}
    else
    {
	for (map<string, unsigned>::iterator it=cells[cellid].ports.begin() ; it!=cells[cellid].ports.end() ; it++)
	{	
		if ((*it).first=="a"||(*it).first=="b")
		{
		   pin* sinknow=&pins[(*it).second];
		   for (unsigned int j=0;j<nets[sinknow->net].sinks.size();++j)
		   {
		   	  if (nets[sinknow->net].sinks[j]==(sinknow)->id)
		   	  {
		         nets[sinknow->net].lateWeights[j]*=((sinknow->lateAt)/pins[cells[cellid].ports["o"]].lateAt);
		         nets[sinknow->net].earlyWeights[j]*=(pins[cells[cellid].ports["o"]].earlyAt/(sinknow->earlyAt));
		      }
		   }
		}
		else if ((*it).first=="d") //DFF
		{
		   pin* sinknow=&pins[(*it).second];
		   for (unsigned int j=0;j<nets[sinknow->net].sinks.size();++j)
		   {
		   	  if (nets[sinknow->net].sinks[j]==(sinknow)->id)
		   	  {
		         nets[sinknow->net].lateWeights[j]*=((sinknow->lateAt)/pins[cells[cellid].ports["q"]].lateAt);
		         nets[sinknow->net].earlyWeights[j]*=(pins[cells[cellid].ports["q"]].earlyAt/(sinknow->earlyAt));
		      }
		   }
		}
		
	}
    }
	
}

////for not PO, if for PO need modify?
void circuit::netWeightDistribute(int cellid)
{
    double inEarlyWeight=0.0;
    double outEarlyWeight=0.0;
    double inLateWeight=0.0;
    double outLateWeight=0.0;

    for (map<string, unsigned>::iterator it=cells[cellid].ports.begin() ; it!=cells[cellid].ports.end() ; it++)
	{		
		   pin* sinknow=&pins[(*it).second];
		   for (unsigned int j=0;j<nets[sinknow->net].sinks.size();++j)
		   {		   	  		   	  
		   	  	if ((*it).first=="a"||(*it).first=="b"||(*it).first=="d")
		        {
		         	if (nets[sinknow->net].sinks[j]==(sinknow)->id)
		   	        {
		               inLateWeight+=nets[sinknow->net].lateWeights[j];
		               inEarlyWeight+=nets[sinknow->net].earlyWeights[j];
		            }
		        }
		        else if ((*it).first=="o"||(*it).first=="q")
		        {
		         	outLateWeight+=nets[sinknow->net].lateWeights[j];
		            outEarlyWeight+=nets[sinknow->net].earlyWeights[j];
		        }
		      
		   }	
	}
#ifdef DEBUG
	cout<<"cell "<<cells[cellid].name<<" have early in total: "
	<<inEarlyWeight<<" and late in total "<<inLateWeight<<endl;
	cout<<"cell "<<cells[cellid].name<<" have early out total: "
	<<outEarlyWeight<<" and late out total "<<outLateWeight<<endl;
#endif
	for (map<string, unsigned>::iterator it=cells[cellid].ports.begin() ; it!=cells[cellid].ports.end() ; it++)
	{
		pin* sinknow=&pins[(*it).second];
		for (unsigned int j=0;j<nets[sinknow->net].sinks.size();++j)
		{
			if (nets[sinknow->net].sinks[j]==(sinknow)->id)
		   	  {
		   	  	 if ((*it).first=="a"||(*it).first=="b"||(*it).first=="d")
		         {
#ifdef DEBUG		    
		         	double gg=nets[sinknow->net].earlyWeights[j];
		         	double ggg=nets[sinknow->net].lateWeights[j];
		         	string ii=nets[sinknow->net].name;
		         	cout<<"changing net "<<ii<<" earlyWeight "<<gg<<" to "<<(gg*(outEarlyWeight/inEarlyWeight))<<endl;
		         	cout<<"changing net "<<ii<<" lateWeight "<<ggg<<" to "<<(ggg*(outLateWeight/inLateWeight))<<endl;
#endif
		           nets[sinknow->net].lateWeights[j]*=(outLateWeight/inLateWeight);
		           nets[sinknow->net].earlyWeights[j]*=(outEarlyWeight/inEarlyWeight);
		         }
		         
		      }
		}
	}
}