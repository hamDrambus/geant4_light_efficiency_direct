#include "SimulationSummary.hh"

SimulationSummary::SimulationSummary(void)
{
	ClearAllData();
}

SimulationSummary::~SimulationSummary()
{}

void SimulationSummary::ClearAllData(void)
{
	tot_probability = -1;
	small_prob = -1;
	small_weight = 0;
	curr_weight_pbs = 0;
}

G4int SimulationSummary::Num_of_events()
{
	return curr_weight_pbs*NUM_OF_PROBS_STORE + small_weight;
}

std::ostream& operator<<(std::ostream& str, const SimulationSummary& data)
{
	str << "total efficieency: " << data.tot_probability<< G4endl;
	str << "Number of events " << data.curr_weight_pbs*NUM_OF_PROBS_STORE + data.small_weight << G4endl;
	return str;
}

void SimulationSummary::SetHit(G4int is_hit, G4double prob)
{
	if (is_hit == -1)// do not account
		return;
	small_prob = (small_prob*small_weight + is_hit*prob) / (small_weight + 1);
	small_weight+=1;
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