#include "B1SteppingAction.hh"
#include "B1EventAction.hh"
#include "B1DetectorConstruction.hh"

#include "G4Step.hh"
#include "G4Event.hh"
#include "CustomRunManager.hh"
#include "G4LogicalVolume.hh"


B1SteppingAction::B1SteppingAction(B1EventAction* eventAction)
: G4UserSteppingAction(),
  fEventAction(eventAction),
  fScoringVolume(0),
  stopPhotonVolume(0)
{}


B1SteppingAction::~B1SteppingAction()
{}

void B1SteppingAction::UserSteppingAction(const G4Step* step)
{

  // get volume of the current step
	G4VPhysicalVolume* tremor = step->GetPostStepPoint()->GetTouchableHandle()->GetVolume();
	if (NULL==tremor)
	{// kind of error
		CustomRunManager* manman = (CustomRunManager*)(G4RunManager::GetRunManager());
		manman->SetHit(-1,1,step);  //was manman->SetHit(0);
		//step->GetTrack()->SetTrackStatus(fStopAndKill);
		//manman->next_event(step);
		return;
	}
  G4LogicalVolume* post_log_volume = tremor->GetLogicalVolume();
  CustomRunManager* manman = (CustomRunManager*)(G4RunManager::GetRunManager());
#ifdef NO_BACK_SCATTERING
  if (!manman->select_photon_BP(step, G4ThreeVector(0, 0, 0), G4ThreeVector(0, 0, 0)))
  {
	  step->GetTrack()->SetTrackStatus(fStopAndKill); //kills if photon goes inside interior volume
  }
#endif
  manman->sim_results->SteppingProc(step);
  if (step->GetPostStepPoint()->GetLocalTime() > 1000)
  {
	  G4cout << "Stuck photon, seting 0 hit"<<std::endl;
	  step->GetTrack()->SetTrackStatus(fStopAndKill);
	  manman->SetHit(-2);
	  return;
  }
  if (step->GetTrack()->GetTrackStatus() == fStopAndKill) //TODO: in case some process killed the track, I need to handle events properly
  {
	  manman->SetHit(0);
	  /*manman->next_event(step);*/
	  return;
  }
  if (manman->is_internal_reflection&&manman->is_no_absorbtion)
  {
	  step->GetTrack()->SetTrackStatus(fStopAndKill);
	  manman->SetHit(0,1,step);
	  manman->is_internal_reflection = 0;
	  manman->is_no_absorbtion = 0;
	  return;
  }
  manman->is_internal_reflection = 0;
  manman->is_no_absorbtion = 0;
  if (step->GetTrack()->GetTrackStatus() == fKillTrackAndSecondaries) //in case some process fails
  {
	  manman->SetHit(-1,1,step);
	  return;
  }
  //if (manman->get_curr_event_probab()< MIN_ALLOWED_STEPPING)
  //{
	 // manman->SetHit(0);
	 // step->GetTrack()->SetTrackStatus(fStopAndKill);
	 // manman->next_event(step);
	 // return;
  //}
  //WARNING! changed casting here
  B1DetectorConstruction* detectorConstruction = (B1DetectorConstruction*) (manman->GetUserDetectorConstruction());
  G4double detect_prob = detectorConstruction->GetHitProbability(step->GetPostStepPoint());
  if (0 == detect_prob)
  {
	  manman->SetHit(0,1,step);
	  step->GetTrack()->SetTrackStatus(fStopAndKill);
	  //manman->next_event(step);
	  return;
  }
  if (1 == detect_prob)
  {
	  manman->SetHit(1,1,step);
	  step->GetTrack()->SetTrackStatus(fStopAndKill);
#ifdef TOP_MESH_TEST
	  if (step->GetPostStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume() == detectorConstruction->top_mesh_test_detector)
		  manman->on_hit_proc(step->GetPostStepPoint()->GetPosition(), detect_prob,1);
	  if (step->GetPostStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume() == detectorConstruction->LAr_layer)
		  manman->on_hit_proc(step->GetPostStepPoint()->GetPosition(), detect_prob, 0);
#endif
	  //manman->next_event(step);
	  return;
  }
  if (-1 == detect_prob)
  {
	  //manman->next_event(step);
	  return;
  }
  //if detection probability is not a 0 or 1, then it's treated as process
#ifdef TOP_MESH_TEST
  if (step->GetPostStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume() == detectorConstruction->top_mesh_test_detector)
	  manman->on_hit_proc(step->GetPostStepPoint()->GetPosition(), detect_prob,1);
  if (step->GetPostStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume() == detectorConstruction->LAr_layer)
	  manman->on_hit_proc(step->GetPostStepPoint()->GetPosition(), detect_prob, 0);
#endif
  //manman->SetPhEvProb(detect_prob);
  //manman->SetPhEvType(RM_PHOTON_DETECTION);
  //manman->process_end(step);
  //after "detection process" hit as usual
  manman->SetHit(1, detect_prob,step);
  step->GetTrack()->SetTrackStatus(fStopAndKill);
  //manman->next_event(step);
}

