#include "SimulationSummary.hh"


event_history::event_history()
{
	clear();
	weight = 0;
	for (int Y = 0; Y < 4; Y++)
	{
		tot_probability[Y] = 0;
		neighbours_prob[Y] = 0;
		opposite_prob[Y] = 0;
		neighbours_not_refl_prob[Y] = 0;
		opposite_not_refl_prob[Y] = 0;
		same_prob[Y] = 0;
	}
}

void event_history::clear(void)
{
	is_reemissed = false;
	wls_origin_number = -1; //0 x+, 1 y+, 2 x-, 3 y-
	detection_number = -1;	//--||--
	bot_cu_refl_num = 0;
	top_cu_refl_num = 0;
	//#ifdef TEMP_CODE_ //for gammas from WLS study
	//	wls_origin_number=0;
	//	is_reemissed = true;
	//#endif
	one_run_hits.erase(one_run_hits.begin(), one_run_hits.end());
}

void event_history::SetHit(G4int is_hit, G4double prob, const G4Step* step)
{
	if ((prob < 0) || (is_hit < 0))
		return;
	one_run_hits.push_back(is_hit*prob);
}

void event_history::SteppingProc(const G4Step* step)//sets gamma 'tags'
{
	CustomRunManager* manman = (CustomRunManager*)(G4RunManager::GetRunManager());
	B1DetectorConstruction *det = (B1DetectorConstruction *)manman->GetUserDetectorConstruction();
	G4LogicalVolume* postV = step->GetPostStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume();
	G4LogicalVolume* preV = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume();
	if ((0 != step->GetSecondary()->size()) && ((preV == det->Xp_wls) || (preV == det->Xn_wls) || (preV == det->Yp_wls) || (preV == det->Yn_wls)))
	{
		if (preV == det->Xp_wls)
			wls_origin_number = 0;
		if (preV == det->Yp_wls)
			wls_origin_number = 1;
		if (preV == det->Xn_wls)
			wls_origin_number = 2;
		if (preV == det->Yn_wls)
			wls_origin_number = 3;
		is_reemissed = true;
	}
	if (!(0>det->GetHitProbability(step->GetPostStepPoint())))
	{
		if (postV == det->Xp_PMT)
			detection_number = 0;
		if (postV == det->Yp_PMT)
			detection_number = 1;
		if (postV == det->Xn_PMT)
			detection_number = 2;
		if (postV == det->Yn_PMT)
			detection_number = 3;
	}
	if ((postV == det->top_GEM->parent) || (postV == det->top_cu_plate))
		top_cu_refl_num++;
	if ((postV == det->bot_GEM->parent) || (postV == det->bot_cu_plate))
		bot_cu_refl_num++;
}

void event_history::RunEnd()
{
	G4double hit_probab = one_run_hits.empty() ? -1 : 0;
	for (auto j = one_run_hits.begin(); j != one_run_hits.end(); j++)
		hit_probab += *j;
	hit_probab = hit_probab > 1 ? 1 : hit_probab;
	if (hit_probab < 0)
	{
		clear();
		return;
	}
	for (int p = 0; p < 4; p++)
		tot_probability[p] = (tot_probability[p] * weight + (p == detection_number ? hit_probab : 0)) / (weight + 1);
	weight++;
	hit_probab = one_run_hits.back();
	G4double ne_pb = 0;
	G4double op_pb = 0;
	G4double ne_nr_pb = 0;
	G4double op_nr_pb = 0;
	G4double sm_pb = 0;
	if (hit_probab < 0)
		goto Out;
	if (is_reemissed)
	{
		if (wls_origin_number == detection_number)
		{
			sm_pb = hit_probab;
			goto Out;
		}
		G4int num1 = wls_origin_number;
		num1++;
		num1 = num1 > 3 ? num1 - 4 : num1;
		G4int num2 = wls_origin_number;
		num2--;
		num2 = num2 <0 ? num2 + 4 : num2;
		if ((num1 == detection_number) || (num2 == detection_number))
		{
			ne_pb = hit_probab;
			if ((top_cu_refl_num == 0) && (bot_cu_refl_num == 0))
				ne_nr_pb = hit_probab;
			goto Out;
		}
		num1 = wls_origin_number + 2;
		num1 = num1 > 3 ? num1 - 4 : num1;
		if (num1 == detection_number)
		{
			op_pb = hit_probab;
			if ((top_cu_refl_num == 0) && (bot_cu_refl_num == 0))
				op_nr_pb = hit_probab;
		}
	}
Out:
	for (int p = 0; p < 4; p++)
	{
		neighbours_prob[p] = (neighbours_prob[p] * (weight - 1) + (detection_number == p ? ne_pb : 0)) / weight;
		opposite_prob[p] = (opposite_prob[p] * (weight - 1) + (detection_number == p ? op_pb : 0)) / weight;
		neighbours_not_refl_prob[p] = (neighbours_not_refl_prob[p] * (weight - 1) + (detection_number == p ? ne_nr_pb : 0)) / weight;
		opposite_not_refl_prob[p] = (opposite_not_refl_prob[p] * (weight - 1) + (detection_number == p ? op_nr_pb : 0)) / weight;
		same_prob[p] = (same_prob[p] * (weight - 1) + (detection_number == p ? sm_pb : 0)) / weight;
	}
	clear();
}

std::ostream& operator<<(std::ostream& str, const event_history& data)
{
	str << G4endl;
	str << "____event_history____" << G4endl;
	str << "Number of events: " << data.weight << G4endl;
	G4double tot_pb = 0, neib_pb = 0, opp_pb = 0, same_pb = 0, nr_neib_pb = 0, nr_opp_pb = 0;
	for (int p = 0; p < 4; p++)
	{
		tot_pb += data.tot_probability[p];
		neib_pb += data.neighbours_prob[p];
		opp_pb += data.opposite_prob[p];
		same_pb += data.same_prob[p];
		nr_neib_pb += data.neighbours_not_refl_prob[p];
		nr_opp_pb += data.opposite_not_refl_prob[p];
	}
	str << "Total probability: " << tot_pb << G4endl;
	str << "To neighbour PMTs probability: " << neib_pb << G4endl;
	str << "To opposite PMT prob: " << opp_pb << G4endl;
	str << "To the same PMT prob: " << same_pb << G4endl;
	str << "No Cu reflection to neighbour PMTs probab: " << nr_neib_pb << G4endl;
	str << "No Cu reflection to opposite PMT probab: " << nr_opp_pb << G4endl;
	str << "----------------------------" << G4endl;
	str << G4endl;
	str << "X_Positive probabilities:" << G4endl;
	str << "Total probability: " << data.tot_probability[0] << G4endl;
	str << "To neighbour PMTs probability: " << data.neighbours_prob[0] << G4endl;
	str << "To opposite PMT prob: " << data.opposite_prob[0] << G4endl;
	str << "To the same PMT prob: " << data.same_prob[0] << G4endl;
	str << "No Cu reflection to neighbour PMTs probab: " << data.neighbours_not_refl_prob[0] << G4endl;
	str << "No Cu reflection to opposite PMT probab: " << data.opposite_not_refl_prob[0] << G4endl;
	str << G4endl;
	str << "----------------------------" << G4endl;
	str << "Y_Positive probabilities:" << G4endl;
	str << "Total probability: " << data.tot_probability[1] << G4endl;
	str << "To neighbour PMTs probability: " << data.neighbours_prob[1] << G4endl;
	str << "To opposite PMT prob: " << data.opposite_prob[1] << G4endl;
	str << "To the same PMT prob: " << data.same_prob[1] << G4endl;
	str << "No Cu reflection to neighbour PMTs probab: " << data.neighbours_not_refl_prob[1] << G4endl;
	str << "No Cu reflection to opposite PMT probab: " << data.opposite_not_refl_prob[1] << G4endl;
	str << G4endl;
	str << "----------------------------" << G4endl;
	str << "X_Negative probabilities:" << G4endl;
	str << "Total probability: " << data.tot_probability[2] << G4endl;
	str << "To neighbour PMTs probability: " << data.neighbours_prob[2] << G4endl;
	str << "To opposite PMT prob: " << data.opposite_prob[2] << G4endl;
	str << "To the same PMT prob: " << data.same_prob[2] << G4endl;
	str << "No Cu reflection to neighbour PMTs probab: " << data.neighbours_not_refl_prob[2] << G4endl;
	str << "No Cu reflection to opposite PMT probab: " << data.opposite_not_refl_prob[2] << G4endl;
	str << G4endl;
	str << "----------------------------" << G4endl;
	str << "Y_Negative probabilities:" << G4endl;
	str << "Total probability: " << data.tot_probability[3] << G4endl;
	str << "To neighbour PMTs probability: " << data.neighbours_prob[3] << G4endl;
	str << "To opposite PMT prob: " << data.opposite_prob[3] << G4endl;
	str << "To the same PMT prob: " << data.same_prob[3] << G4endl;
	str << "No Cu reflection to neighbour PMTs probab: " << data.neighbours_not_refl_prob[3] << G4endl;
	str << "No Cu reflection to opposite PMT probab: " << data.opposite_not_refl_prob[3] << G4endl;
	str << G4endl;
	return str;
}


SimulationSummary::SimulationSummary(void)
{
	ev_history = new event_history;
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
	delete ev_history;
}

void SimulationSummary::ClearAllData(void)
{
	num_of_stuck_photons = 0;
	tot_probability = -1;
	small_prob = -1;
	small_weight = 0;
	curr_weight_pbs = 0;
#ifdef TEST_WLS_SPECTRA
	for (int j = 0; j < EN_BINS_N_; j++)
		wl_spec[j] = 0;
#endif
	ev_history->clear();
}

G4int SimulationSummary::Num_of_events()
{
	return curr_weight_pbs*NUM_OF_PROBS_STORE + small_weight;
}

std::ostream& operator<<(std::ostream& str, const SimulationSummary& data)
{
	str << "total efficieency: " << data.tot_probability<< G4endl;
	str << "Number of events " << data.curr_weight_pbs*NUM_OF_PROBS_STORE + data.small_weight << G4endl;
	str << "Number stuck photons is" << data.num_of_stuck_photons<< G4endl;
#ifdef TEST_WLS_SPECTRA
	std::ofstream file;
	file.open(DETECTED_SPECTRUM, std::ios_base::trunc);
	G4double min_wl = 1.2398e3*eV / WLS_MAX_EN;
	G4double max_wl = 1.2398e3*eV / WLS_MIN_EN;
	for (int j = 0; j < EN_BINS_N_; j++)
		file << min_wl + j*(max_wl - min_wl) / (EN_BINS_N_ - 1) << "\t" << data.wl_spec[j] << std::endl;
	file.close();
#endif
	str << *(data.ev_history);
	return str;
}

void SimulationSummary::SetHit(G4int is_hit, G4double prob,const G4Step* step)
{
	ev_history->SetHit(is_hit, prob, step);
	if ((is_hit == -1)||(prob<0))// do not account
		return;
	if (is_hit == -2)
		num_of_stuck_photons++;
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
	ev_history->RunEnd();
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

void SimulationSummary::SteppingProc(const G4Step* step)
{
	ev_history->SteppingProc(step);
}