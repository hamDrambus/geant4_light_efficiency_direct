#include "B1DetectorConstruction.hh"
#include "B1ActionInitialization.hh"

#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "CustomRunManager.hh"
#endif

#include "G4UImanager.hh"
//#include "QBBC.hh"
#include "OpticsPhysicsList.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif

#include "Randomize.hh"
#include "G4VSteppingVerbose.hh"

#include <time.h>
#include <windows.h>

int main(int argc,char** argv)
{
  // Choose the Random engine
  //
  G4Random::setTheEngine(new CLHEP::RanecuEngine(41));
  
  // Construct the default run manager
  //
#ifdef G4MULTITHREADED
  G4MTRunManager* runManager = new G4MTRunManager;
#else
  CustomRunManager* runManager = new CustomRunManager;
#endif

  // Set mandatory initialization classes:
  // Detector construction
  runManager->SetVerboseLevel(0);
  runManager->SetUserInitialization(new B1DetectorConstruction());
  // Physics list
  G4VUserPhysicsList* physicsList = new OpticPhysicsList();
  //physicsList->SetVerboseLevel(1);
  runManager->SetUserInitialization(physicsList);
    
  // User action initialization
  runManager->SetUserInitialization(new B1ActionInitialization());

  // Initialize G4 kernel
  //
  runManager->Initialize();
  //runManager->BeamOn(runManager->x_num*runManager->y_num);
#ifdef AVERAGE_QE_
  std::list<G4double> wavelengths;
  std::list<G4double> probs;
  G4double min_wl, max_wl;
  std::ifstream fl(WLS_SPECTRUM_FILE);
  while (!fl.eof())
  {
	  G4double x, y;
	  fl >> x;
	  if (fl.eof())
		  break;
	  fl >> y;
	  wavelengths.push_back(x);
	  probs.push_back(y);
}
  fl.close();
  min_wl = wavelengths.front(); //sorted input assumed
  max_wl = wavelengths.back();
  G4double* temp_wavelengths = new G4double[wavelengths.size()];
  G4double* temp_weights = new G4double[wavelengths.size()];
  int fgg = 0;
  auto j = probs.begin();
  for (auto u = wavelengths.begin(); u != wavelengths.end(); fgg++, j++, u++)
  {
	  temp_wavelengths[fgg] = *u;
	  temp_weights[fgg] = *j;
  }
  G4MaterialPropertyVector wls_file_vals = G4MaterialPropertyVector(temp_wavelengths, temp_weights, wavelengths.size());
  delete[] temp_wavelengths;
  delete[] temp_weights;
  wavelengths.erase(wavelengths.begin(),wavelengths.end());
  probs.erase(probs.begin(), probs.end());
  //-----------------------
  fl.open(PMT_QE_FILE);
  while (!fl.eof())
  {
	  G4double x, y;
	  fl >> x;
	  if (fl.eof())
		  break;
	  fl >> y;
	  wavelengths.push_back(x);
	  probs.push_back(y);
  }
  fl.close();
  temp_wavelengths = new G4double[wavelengths.size()];
  temp_weights = new G4double[wavelengths.size()];
  fgg = 0;
  auto s = probs.begin();
  for (auto u = wavelengths.begin(); u != wavelengths.end(); fgg++, s++, u++)
  {
	  temp_wavelengths[fgg] = *u;
	  temp_weights[fgg] = *s;
  }
  G4MaterialPropertyVector qe_file_vals = G4MaterialPropertyVector(temp_wavelengths, temp_weights, wavelengths.size());
  delete[] temp_wavelengths;
  delete[] temp_weights;
  wavelengths.erase(wavelengths.begin(), wavelengths.end());
  probs.erase(probs.begin(), probs.end());

  G4double integral_wls = 0, convolution=0;
  for (int h = 0; h < 10000; h++)
  {
	  G4double wl = min_wl + (max_wl - min_wl)*h / (10000 - 1);
	  G4double qe = qe_file_vals.Value(wl);
	  G4double w_wl = wls_file_vals.Value(wl);
	  integral_wls += w_wl;
	  convolution += w_wl*qe;
  }
  G4cout << "Average QE: " << convolution/integral_wls<< G4endl;
#endif
  time_t timer_start, timer_end;
  time(&timer_start);
#if !defined(SPATIAL_ANGLE_)
#ifdef AR_SPEC_TEST
	  runManager->BeamOn(700000);
#else
  runManager->BeamOn(50000000);
#endif
#else
  runManager->BeamOn(5000000);
#endif
  //DONE: TODO: check reemiss/noreemiss code
  //G4double pb_total, pb_no_reemiss, pb_reemissed;
  //runManager->get_total_detetion_eff(&pb_no_reemiss, &pb_reemissed, &pb_total);
  time(&timer_end);
  G4double seconds = difftime(timer_end, timer_start);
  G4cout<<*(runManager->sim_results);
  G4cout << *(runManager->ev_history);
  //DONE - no, no problem: TODO figure out: messed up with the names for some reason
  //G4cout << "total events proceded: " << (runManager->sim_results->Num_of_events())+(runManager->extra_run_id) << G4endl;
  G4cout << "Time elapsed: "<< ((int)seconds/3600)<<"h "<<((int)seconds % 3600)/60<<"m "<<(int)seconds % 60<<"s"<<G4endl;
#if !defined(SPATIAL_ANGLE_)
#ifdef AR_SPEC_TEST
  runManager->get_detected_spectrum();
#endif
#ifdef TOP_MESH_TEST
  runManager->export_to_bmp(&runManager->top_hits_xs, &runManager->top_hits_ys, &runManager->top_hits_probs, "top_mesh_test.bmp");
  runManager->export_to_bmp(&runManager->bot_hits_xs, &runManager->bot_hits_ys, &runManager->bot_hits_probs, "bot_mesh_test.bmp");
#endif
#endif
  Beep(1500, 300);
  Sleep(150);
  Beep(1500, 300);
  Sleep(150);
  Beep(1500, 300);
#ifdef G4VIS_USE
  // Initialize visualization
  G4VisManager* visManager = new G4VisExecutive;
  // G4VisExecutive can take a verbosity argument - see /vis/verbose guidance.
  // G4VisManager* visManager = new G4VisExecutive("Quiet");
  visManager->Initialize();
#endif

  // Get the pointer to the User Interface manager
  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  if (argc!=1) {
    // batch mode
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UImanager->ApplyCommand(command+fileName);
  }
  else {
    // interactive mode : define UI session
#ifdef G4UI_USE
    G4UIExecutive* ui = new G4UIExecutive(argc, argv,"qt");
	//G4UIExecutive* ui = new G4UIExecutive(0, 0);
#ifdef G4VIS_USE
    UImanager->ApplyCommand("/control/execute init_vis.mac"); 
#else
    UImanager->ApplyCommand("/control/execute init.mac"); 
#endif
    ui->SessionStart();
    delete ui;
#endif
  }

#ifdef G4VIS_USE
  delete visManager;
#endif
  delete runManager;

  return 0;
}
