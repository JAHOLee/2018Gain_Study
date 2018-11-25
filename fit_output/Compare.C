void Read_Module_List(string Module_Layout , int config);
int Module_List[28];
std::map<int,int> moduleID2BDorder;
string Module_Layout = "../configs/all_config.csv";
int config = 1;

void Compare(){
  
  int MAXBD = 28;
  int MAXCHIP = 4;
  int MAXCH = 32;
  const int Calib_files = 5;
  double HLcoeff[Calib_files][MAXBD][MAXCHIP*MAXCH];
  double LTcoeff[Calib_files][MAXBD][MAXCHIP*MAXCH];
  double TOToffset[Calib_files][MAXBD][MAXCHIP*MAXCH];

  double dummy_ind[MAXCHIP*MAXCH];
  for(int i = 0 ; i < MAXCHIP*MAXCH; ++i)
    dummy_ind[i] = i;
  //  Read_Module_List(Module_Layout,config);
  
  char fileN[50],title[100];
  sprintf(fileN,"./Chia-hung_TB_Calib.txt");
  TTree *tree = new TTree() ;
  tree->ReadFile(fileN,"L_ID/I:M_ID/I:S_ID/I:C_ID/I:A2M/D:L2HT/D:L2H/D:T2LT/D:T2L/D:TOFF/D:HLTYPE/O:LTTYPE/O");
  
  int points = tree->GetEntries();

  int June_config[] = { 78, 90, 89, 88, 77, 85, 84, 32, 69, 79,
			67, 65, 76, 83, 35, 36, 70, 73, 44, 51,
			86, 87, 54, 62, 64, 55, 59, 71 };
  cout << "hu" << endl;
  //Fill result based on Oct config1 layerID
  for(int i = 0 ; i < points ; ++i){
    tree->GetEntry(i) ;
    int M_ID = tree->GetLeaf("M_ID")->GetValue(0);
    int L_ID = moduleID2BDorder.find(M_ID)->second;
    if(L_ID == 0 && M_ID != 78){ continue; }
    //cout << L_ID << " L_ID is available!" << endl;
    int S_ID = tree->GetLeaf("S_ID")->GetValue(0);
    int C_ID = tree->GetLeaf("C_ID")->GetValue(0);
    double L2H  = tree->GetLeaf("L2H" )->GetValue(0);
    double T2L  = tree->GetLeaf("T2L" )->GetValue(0);
    double TOFF  = tree->GetLeaf("TOFF" )->GetValue(0);

    HLcoeff[0][L_ID][S_ID*32+C_ID/2] = L2H;
    LTcoeff[0][L_ID][S_ID*32+C_ID/2] = T2L;
    TOToffset[0][L_ID][S_ID*32+C_ID/2] = TOFF;
  }


  
  for(int my_calib = 1 ; my_calib < Calib_files ; ++my_calib ){
    if(my_calib == 1)
      sprintf(fileN,"TPro_config1_fittingoutput.txt");
    else if(my_calib == 2)
      sprintf(fileN,"TPro_config1_v2_fittingoutput.txt");
    else if(my_calib == 3)
      sprintf(fileN,"TPro_config1_v3_fittingoutput.txt");
    else
      sprintf(fileN,"TPro_config1_v4_fittingoutput.txt");

    TTree *tree = new TTree() ;
    tree->ReadFile(fileN,"L_ID/I:M_ID/I:S_ID/I:C_ID/I:A2M/D:L2HT/D:L2H/D:T2LT/D:T2L/D:TOFF/D:HLTYPE/O:LTTYPE/O:HGLGF/I:LGTOT/I");
    points = tree->GetEntries();
    for(int i = 0 ; i < points ; ++i){
      tree->GetEntry(i) ;

      int L_ID = tree->GetLeaf("L_ID")->GetValue(0);
      int M_ID = tree->GetLeaf("M_ID")->GetValue(0);
      int S_ID = tree->GetLeaf("S_ID")->GetValue(0);
      int C_ID = tree->GetLeaf("C_ID")->GetValue(0);
      double L2H  = tree->GetLeaf("L2H" )->GetValue(0);
      double T2L  = tree->GetLeaf("T2L" )->GetValue(0);
      double TOFF  = tree->GetLeaf("TOFF" )->GetValue(0);

      HLcoeff[my_calib][L_ID][S_ID*32+C_ID/2] = L2H;
      LTcoeff[my_calib][L_ID][S_ID*32+C_ID/2] = T2L;
      TOToffset[my_calib][L_ID][S_ID*32+C_ID/2] = TOFF;
    }
  }

  
  TCanvas *c1 = new TCanvas();

  TGraph *gr[Calib_files];
  TGraph *gr2[Calib_files];
  
  for(int BD = 0 ; BD < MAXBD ; ++BD){
    TMultiGraph *mgr  = new TMultiGraph();
    TLegend *leg = new TLegend(0.67,0.57,0.89,0.89);
    leg->SetBorderSize(0);
    leg->SetHeader("GainFactor","C");

    for(int type = 0 ; type < Calib_files ; ++type){
      gr[type] = new TGraph(MAXCHIP*MAXCH,dummy_ind,HLcoeff[type][BD]);
      gr[type]->SetMarkerStyle(21);
      gr[type]->SetMarkerSize(1);
      gr[type]->SetMarkerColor(1+type);
      if(type == 4)
	gr[type]->SetMarkerColor(6);
      gr[type]->SetLineWidth(2.5);
      if(type == 0)
	sprintf(title,"LG2HG JuneTB(v11)");
      else if(type == 1)
	sprintf(title,"LG2HG OctTB(v1)");
      else if(type == 2)
	sprintf(title,"LG2HG OctTB(v2)");
      else if(type == 3)
	sprintf(title,"LG2HG OctTB(v3)");
      else
	sprintf(title,"LG2HG OctTB(v4)");
      leg->AddEntry(gr[type],title,"P");

      mgr->Add(gr[type],"AP");}
    
    for(int type = 0 ; type < Calib_files ; ++type){
      gr2[type] = new TGraph(MAXCHIP*MAXCH,dummy_ind,LTcoeff[type][BD]);
      gr2[type]->SetMarkerStyle(22);
      gr2[type]->SetMarkerSize(1);
      gr2[type]->SetMarkerColor(1+type);
      if(type == 4)
	gr2[type]->SetMarkerColor(6);
      gr2[type]->SetLineWidth(2.5);

      if(type == 0)
	sprintf(title,"TOT2LG JuneTB(v11)");
      else if(type == 1)
	sprintf(title,"TOT2LG OctTB(v1)");
      else if(type == 2)
	sprintf(title,"TOT2LG OctTB(v2)");
      else if(type == 3)
	sprintf(title,"TOT2LG OctTB(v3)");
      else
	sprintf(title,"TOT2LG OctTB(v4)");
      
      leg->AddEntry(gr2[type],title,"P");
  
      mgr->Add(gr2[type],"AP");}
    
    mgr->Draw("AP");
    sprintf(title,"Layer%d(module%d)_GainFactor",BD,Module_List[BD]);
    mgr->SetTitle(title);
    mgr->SetMaximum(15);
    mgr->SetMinimum(0);
    mgr->GetXaxis()->SetTitle("chip*32+chID/2");
    mgr->GetYaxis()->SetTitle("GainFactor");
    leg->Draw("same");
    c1->Update();
    //c1->WaitPrimitive();
    char p_name[100];
    sprintf(p_name,"plots/%s.png",title);
    c1->SaveAs(p_name);
  }

  /*
  for(int BD = 0 ; BD < MAXBD ; ++BD){
    TMultiGraph *mgr  = new TMultiGraph();
    TLegend *leg = new TLegend(0.62,0.13,0.89,0.35);
    leg->SetBorderSize(0);
    leg->SetHeader("TOT offset","C");

    for(int type = 0 ; type < Calib_files ; ++type){
      gr[type] = new TGraph(MAXCHIP*MAXCH,dummy_ind,TOToffset[type][BD]);
      gr[type]->SetMarkerStyle(21);
      gr[type]->SetMarkerSize(1);
      gr[type]->SetMarkerColor(1+type);
      gr[type]->SetLineWidth(2.5);
      if(type == 0)
	sprintf(title,"TOT2LG Inj(mip dac)");
      else if(type == 1)
	sprintf(title,"TOT2LG JuneTB(before)");
      else if(type == 2)
	sprintf(title,"TOT2LG Inj(v11)");
      else if(type == 3)
	sprintf(title,"TOT2LG JuneTB(v11)");
      else
	sprintf(title,"TOT2LG OctTB(v1)");
      leg->AddEntry(gr[type],title,"P");

      mgr->Add(gr[type],"AP");}

    mgr->Draw("AP");
    sprintf(title,"Layer%d(module%d)_TOTOffset",BD,layer_to_moduleID[BD]);
    mgr->SetTitle(title);
    mgr->SetMinimum(0);
    mgr->GetXaxis()->SetTitle("chip*32+chID/2");
    mgr->GetYaxis()->SetTitle("Offset(ADC)");
    leg->Draw("same");
    c1->Write(title);
    sprintf(title,"%s.png",title);
    c1->SaveAs(title);
    
  }
  */
  
}
void Read_Module_List(string Module_Layout , int config){  
  ifstream infile(Module_Layout.c_str());
  string line;
  int line_count = 0;
  int members = 6;
  string line_contents[members];
  
  // Get the headers
  getline(infile,line); 
  getline(infile,line);
  int limit = sizeof(Module_List)/sizeof(int);
  while(true){
    getline(infile,line);
    if( infile.eof() || line_count == limit) {break;};
    std::istringstream iss(line);
    for(int i = 0 ; i < members ; ++i){
      getline(iss,line_contents[i], ',' );}
    if(line_contents[(config-1)*2] == ""){ continue; } 
    int ModuleID = std::stoi( line_contents[(config-1)*2] );
    Module_List[line_count] = ModuleID;
    moduleID2BDorder.insert( std::pair<int,int>(ModuleID,line_count) );
    //cout << "Module "<< ModuleID << " correspond to BD " << line_count << endl;
    line_count++;
  }
  infile.close();
  
}
