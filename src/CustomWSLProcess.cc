#include "G4ios.hh"
#include "G4PhysicalConstants.hh"
#include "G4OpProcessSubType.hh"

#include "CustomWLSProcess.hh"
#include "G4GeometryTolerance.hh"

#include "G4VSensitiveDetector.hh"
#include "G4ParallelWorldProcess.hh"

#include "G4SystemOfUnits.hh"
#include "globals.hh"

//Process is fully developed by geffdroid

CustomWLSProcess::CustomWLSProcess(const G4String& processName,
	G4ProcessType type)
	: G4VDiscreteProcess(processName, type)
{
	SetVerboseLevel(0);
	if (verboseLevel > 0) {
		G4cout << GetProcessName() << " is created " << G4endl;
	}

	SetProcessSubType(fOpAbsorption);
	theStatus = Defalt;
}

CustomWLSProcess::~CustomWLSProcess(){}


// PostStepDoIt
// ------------
//

G4VParticleChange*
CustomWLSProcess::PostStepDoIt(const G4Track& aTrack, const G4Step& aStep)
{
#ifdef DEBUG_MC_NODES
	G4cout << "CWLSP: wls PostStepDoIt()" << G4endl;
#endif
	theStatus = Defalt;

	aParticleChange.Initialize(aTrack);
	aParticleChange.ProposeVelocity(aTrack.GetVelocity());
	const G4DynamicParticle* aParticle = aTrack.GetDynamicParticle();
	aParticleChange.ProposeMomentumDirection(aParticle->GetMomentumDirection());
	aParticleChange.ProposePolarization(aParticle->GetPolarization());

	// Get hyperStep from  G4ParallelWorldProcess
	//  NOTE: PostSetpDoIt of this process should be
	//        invoked after G4ParallelWorldProcess!

	const G4Step* pStep = &aStep;
	const G4Step* hStep = G4ParallelWorldProcess::GetHyperStep();

	if (hStep) pStep = hStep;

	//G4bool isOnBoundary =
	//	(pStep->GetPostStepPoint()->GetStepStatus() == fGeomBoundary);
	G4Material* WLS_mat = pStep->GetPreStepPoint()->GetMaterial();
	
	if (verboseLevel > 0) {
		G4cout << "WSL absorbtion simulation" << G4endl;
	}

	if (aTrack.GetStepLength() <= G4GeometryTolerance::GetInstance()
		->GetSurfaceTolerance()/ 2){
		theStatus = StepTooSmallet;
		return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
	}

	G4double photon_energy = aParticle->GetTotalEnergy();

	G4ThreeVector EndPoint = pStep->GetPostStepPoint()->GetPosition();
	G4ThreeVector StartPoint = pStep->GetPreStepPoint()->GetPosition();
#ifdef DEBUG_MC_NODES
	G4cout << "CWLSP: start point: "<<StartPoint << G4endl;
	G4cout << "CWLSP: end point: " << EndPoint << G4endl;
#endif

	G4MaterialPropertiesTable* aMaterialPropertiesTable;
	G4MaterialPropertyVector* abs_len;
	G4MaterialPropertyVector* wls_efficiency;
	G4double absorb_l, wls_eff;

	aMaterialPropertiesTable = WLS_mat->GetMaterialPropertiesTable();
	if (aMaterialPropertiesTable) {
		abs_len = aMaterialPropertiesTable->GetProperty("ABSORBTION_LENGTH");
		wls_efficiency = aMaterialPropertiesTable->GetProperty("WLS_EFFICIENCY");
	}
	else {
		return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
	}
	CustomRunManager* manman = (CustomRunManager*)(G4RunManager::GetRunManager());
	manman->is_no_absorbtion = 1;
	if (abs_len) 
	{
		absorb_l = GetMaterialAbsLength(abs_len,photon_energy);
		if (absorb_l==DBL_MAX)
			return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep); //TODO: error?
	}
	else
		return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
	if (wls_efficiency)
		wls_eff = wls_efficiency->Value(photon_energy);
	else 
		wls_eff = 0;
	manman->is_no_absorbtion = 0;
	G4double path_l=(EndPoint- StartPoint).getR();
	G4double reemiss_prob = wls_eff*(1.0-exp(-path_l / absorb_l));
	G4double continue_prob = exp(-path_l / absorb_l);
	if (G4UniformRand() > continue_prob)
	{
		aParticleChange.ProposeTrackStatus(fStopAndKill);
		if (G4UniformRand() >= wls_eff) //wls_eff is supposed to be less than 1
		{
			aParticleChange.SetNumberOfSecondaries(0);
			return &aParticleChange;
		}
		else
		{
//#ifdef TEMP_CODE_
//			aParticleChange.SetNumberOfSecondaries(0);
//			return &aParticleChange;
//#endif
			aParticleChange.SetNumberOfSecondaries(1);
			G4MaterialPropertyVector* energy_spec = aMaterialPropertiesTable->GetProperty("WLS_ENERGY_SPECTRUM");
			G4double energy = energy_spec ? energy_spec->Value(G4UniformRand()) : aStep.GetPostStepPoint()->GetTotalEnergy();
			G4double phi = CLHEP::twopi*G4UniformRand();
			G4double cos_theta = 2 * (G4UniformRand()) - 1;
			G4ThreeVector momentum_direction(sin(phi)*sqrt(1 - cos_theta*cos_theta), cos(phi)*sqrt(1 - cos_theta*cos_theta), cos_theta);
			G4ThreeVector polarization(sin(phi)*cos_theta, cos(phi)*cos_theta, -sqrt(1 - cos_theta*cos_theta));
			G4ThreeVector perp = momentum_direction.cross(polarization);
			phi = twopi*G4UniformRand();
			polarization = cos(phi) * polarization + sin(phi) * perp;
			polarization = polarization.unit();
			G4double time = aStep.GetPostStepPoint()->GetGlobalTime();
			//\/Generating of new ph. position (exponential distribution)
			G4ThreeVector start_point1 = aStep.GetPreStepPoint()->GetPosition();
			G4ThreeVector start_point2 = aStep.GetPostStepPoint()->GetPosition();
			G4ThreeVector new_photon_pos;
			G4double R = (start_point1 - start_point2).getR();
			if (R == 0) new_photon_pos = start_point1;
			else
			{
				G4double rand_R = -absorb_l*log(1.0 - G4UniformRand()*(1 - exp(-R / absorb_l)));
				rand_R = rand_R / R;
				//modification because of tolerances
				//(one must not simulate emission of photons near the boundaries (that is start_point1 and start_point2)
				//because boundary processes won't trigger at very small (<G4GeometryTolerance::GetInstance()->GetSurfaceTolerance() / 2) distances
				//and particle will be trasferred directly into a new volume)
				G4double Toler = G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();
				if (Toler >= R)
					rand_R = 0.5;//just in the middle then
				else
				{
					Toler = Toler / 2;
					rand_R = Toler / R + rand_R*(R - 2 * Toler);
				}
				//end modification
				new_photon_pos = start_point1 + rand_R*(start_point2 - start_point1);
			}
			//\/assigning parameters
			G4DynamicParticle* aWLSPhoton = new G4DynamicParticle(G4OpticalPhoton::OpticalPhoton(), momentum_direction);
			aWLSPhoton->SetPolarization(polarization.x(), polarization.y(), polarization.z());
			aWLSPhoton->SetKineticEnergy(energy);
			G4Track* secondary_track = new G4Track(aWLSPhoton, time, new_photon_pos);
			secondary_track->SetTouchableHandle(aTrack.GetTouchableHandle());
			secondary_track->SetParentID(aTrack.GetTrackID());
			aParticleChange.AddSecondary(secondary_track);
			return &aParticleChange;
		}
	}
	else
		return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);//no changes
	return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
}

// GetMeanFreePath
// ---------------
//
G4double CustomWLSProcess::GetMeanFreePath(const G4Track&,
	G4double,
	G4ForceCondition* condition)
{
	*condition = Forced;
	return DBL_MAX;
}

G4bool CustomWLSProcess::IsApplicable(const G4ParticleDefinition& aParticleType)
{
	return (&aParticleType == G4OpticalPhoton::OpticalPhoton());
}

//abs_len_property is implied tobe sorted by energy;
G4double CustomWLSProcess::GetMaterialAbsLength(G4MaterialPropertyVector* abs_len_prop, G4double energy)
{
	G4int N = abs_len_prop->GetVectorLength();
	if (N <= 0)
		return DBL_MAX;
	if (N == 1)
		return (*abs_len_prop)(0);//TODO: check 
	//G4double min_en=(*abs_len_prop)[0];
	G4double min_en = (*abs_len_prop).Energy(0);
	G4double max_en = (*abs_len_prop).Energy(N-1);
	if (energy < min_en)
		return (*abs_len_prop)(0);
	if (energy>max_en)
		return (*abs_len_prop)[N - 1];
	return (abs_len_prop->Value(energy));
}
