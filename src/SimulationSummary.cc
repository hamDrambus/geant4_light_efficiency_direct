#include "SimulationSummary.hh"

SimulationSummary::SimulationSummary(void)
{
#ifdef TEST_WLS_SPECTRA
	wl_spec = new G4double[EN_BINS_N_];
#endif
	ClearAllData();
}

SimulationSummary::~SimulationSummary()
{
#ifdef TEST_WLS_SPECTRA
	if (wl_spec) delete[] wl_spec;
#endif
}

void SimulationSummary::ClearAllData(void)
{
	tot_probability = -1;
	small_prob = -1;
	small_weight = 0;
	curr_weight_pbs = 0;
#ifdef TEST_WLS_SPECTRA
	for (int j = 0; j < EN_BINS_N_; j++)
		wl_spec[j] = 0;
#endif
}

G4int SimulationSummary::Num_of_events()
{
	return curr_weight_pbs*NUM_OF_PROBS_STORE + small_weight;
}

std::ostream& operator<<(std::ostream& str, const SimulationSummary& data)
{
	str << "total efficieency: " << data.tot_probability<< G4endl;
	str << "Number of events " << data.curr_weight_pbs*NUM_OF_PROBS_STORE + data.small_weight << G4endl;
#ifdef TEST_WLS_SPECTRA
	std::ofstream file;
	file.open(DETECTED_SPECTRUM, std::ios_base::trunc);
	G4double min_wl = 1.2398e3*eV / WLS_MAX_EN;
	G4double max_wl = 1.2398e3*eV / WLS_MIN_EN;
	for (int j = 0; j < EN_BINS_N_; j++)
		file << min_wl + j*(max_wl - min_wl) / (EN_BINS_N_ - 1) << "\t" << data.wl_spec[j] << std::endl;
	file.close();
#endif
	return str;
}

void SimulationSummary::SetHit(G4int is_hit, G4double prob,const G4Step* step)
{
	if ((is_hit == -1)||(prob<0))// do not account
		return;
	one_run_hits.push_back(is_hit*prob);//because 1g->2g, but needed to be accounted as 1gamma
#ifdef TEST_WLS_SPECTRA
	if ((is_hit) && (step))
	{
		G4double wl = step->GetPostStepPoint()->GetTotalEnergy();
		wl = 1.2398e3*eV / wl;
		G4double min_wl = 1.2398e3*eV / WLS_MAX_EN;
		G4double max_wl = 1.2398e3*eV / WLS_MIN_EN;
		G4int index = (EN_BINS_N_ - 1)*(wl - min_wl) / (max_wl - min_wl);
		if ((index<0) || (index>(EN_BINS_N_ - 1)))
			return;
		wl_spec[index]++;
	}
#endif
	return;
}

void SimulationSummary::RunEnd(void)
{
	G4double hit_probab = one_run_hits.empty() ? -1 : 0;
	for (auto j = one_run_hits.begin(); j != one_run_hits.end(); j++)
		hit_probab += *j;
	one_run_hits.erase(one_run_hits.begin(),one_run_hits.end());
	hit_probab = hit_probab > 1 ? 1 : hit_probab;
	if (hit_probab < 0)
		return;
	small_prob = (small_prob*small_weight + hit_probab) / (small_weight + 1);
	small_weight += 1;
	if (small_weight == NUM_OF_PROBS_STORE)//can't be more, error otherwise
	{
		tot_probability = (tot_probability*curr_weight_pbs + (small_prob*small_weight / NUM_OF_PROBS_STORE))
			/ ((small_weight / NUM_OF_PROBS_STORE) + curr_weight_pbs);
		curr_weight_pbs += 1;
		small_weight = 0;
		small_prob = -1;
	}
	return;
}

void SimulationSummary::OnEndSimulation()
{
	if (small_weight != 0)
	{
		tot_probability = (tot_probability*curr_weight_pbs + (small_prob*small_weight / NUM_OF_PROBS_STORE)) 
			/ ((small_weight / NUM_OF_PROBS_STORE) + curr_weight_pbs);
		//curr_weight_pbs += small_weight;
		//small_weight = 0;
		//small_prob = -1;
	}
}