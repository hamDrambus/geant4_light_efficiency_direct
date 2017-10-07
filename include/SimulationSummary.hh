#ifndef Simulation_Summary_hh
#define Simulation_Summary_hh

#include "G4Step.hh"
#include "CustomRunManager.hh"
#include "GlobalDefinitions.hh"
#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include <list>



class event_history
{
public:
	G4bool is_reemissed;
	G4int wls_origin_number; //0 x+, 1 y+, 2 x-, 3 y-
	G4int detection_number;	//--||--
	G4int bot_cu_refl_num;
	G4int top_cu_refl_num;

	G4int weight;
	G4double tot_probability[4];
	G4double neighbours_prob[4];
	G4double opposite_prob[4];
	G4double same_prob[4];
	G4double neighbours_not_refl_prob[4];
	G4double opposite_not_refl_prob[4];

	std::list<G4double> one_run_hits;//typically contains one or two elements
	event_history();
	void SetHit(G4int is_hit, G4double prob, const G4Step* step);
	void SteppingProc(const G4Step* step);//sets gamma 'tags'
	void RunEnd();
	friend  std::ostream& operator<<(std::ostream& str, const event_history& data);
	void clear(void);
};

//stores all simulation processed results (specta, probabilities, etc.)
class SimulationSummary
{
public:
	event_history* ev_history;
	G4double tot_probability;

	G4int curr_weight_pbs; //required for updating all probs and spectra by a new simulation sequence
	//^same as number of already cosidered simulations.
	G4double small_prob;
	G4int small_weight;
	G4int tot_num_of_events;
	G4int num_of_stuck_photons;

	std::list<G4double> one_run_hits;

	virtual G4int Num_of_events();
	virtual void OnEndSimulation();
#ifdef TEST_WLS_SPECTRA
	G4double *wl_spec;
#endif
	SimulationSummary(void);
	~SimulationSummary();
	//virtual void ExportEnSpectraToFile(G4String str1 = "total_spectrum.txt", G4String str2 = "no_reemiss_spectrum.txt",
		//G4String str3 = "reemissed_spectrum.txt");
	//virtual void ExportLSpectraToFile(G4String str1 = "l_total_spectrum.txt", G4String str2 = "l_no_reemiss_spectrum.txt",
		//G4String str3 = "l_reemissed_spectrum.txt");
	virtual void ClearAllData(void);
	virtual void SetHit(G4int is_hit, G4double prob,const G4Step* step);
	void SteppingProc(const G4Step* step);
	virtual void RunEnd(); //called when there are no more secondaries

	friend  std::ostream& operator<<(std::ostream& str, const SimulationSummary& data);
};

#endif