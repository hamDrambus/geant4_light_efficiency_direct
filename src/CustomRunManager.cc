#include "CustomRunManager.hh"
#include "G4UImanager.hh"

event_history::event_history()
{
	clear();
	weight = 0;
	for (int Y = 0; Y < 4; Y++)
	{
		tot_probability[Y] = 0;
		neighbours_prob[Y] = 0;
		opposite_prob[Y]= 0;
		neighbours_not_refl_prob[Y] = 0;
		opposite_not_refl_prob[Y] = 0;
		same_prob[Y] = 0;
	}
}

void event_history::clear(void)
{
	is_reemissed=false;
	wls_origin_number=-1; //0 x+, 1 y+, 2 x-, 3 y-
	detection_number=-1;	//--||--
	bot_cu_refl_num=0;
	top_cu_refl_num=0;
//#ifdef TEMP_CODE_ //for gammas from WLS study
//	wls_origin_number=0;
//	is_reemissed = true;
//#endif
	one_run_hits.erase(one_run_hits.begin(), one_run_hits.end());
}

void event_history::SetHit(G4int is_hit, G4double prob,const G4Step* step)
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
	if ((0!= step->GetSecondary()->size()) && ((preV == det->Xp_wls) || (preV == det->Xn_wls) || (preV == det->Yp_wls) || (preV == det->Yn_wls)))
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
	if ((postV == det->top_GEM->parent)||(postV==det->top_cu_plate))
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
		tot_probability[p] = (tot_probability[p]*weight + (p==detection_number? hit_probab :0)) / (weight + 1);
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
		G4int num2= wls_origin_number;
		num2--;
		num2 = num2 <0  ? num2 + 4 : num2;
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
	str << "Total probability: " << tot_pb<< G4endl;
	str << "To neighbour PMTs probability: " << neib_pb << G4endl;
	str << "To opposite PMT prob: " << opp_pb << G4endl;
	str << "To the same PMT prob: " << same_pb << G4endl;
	str << "No Cu reflection to neighbour PMTs probab: " << nr_neib_pb << G4endl;
	str << "No Cu reflection to opposite PMT probab: " << nr_opp_pb << G4endl;
	str << "----------------------------" << G4endl;
	str << G4endl;
	str << "X_Positive probabilities:" << G4endl;
	str << "Total probability: " << data.tot_probability[0]<< G4endl;
	str << "To neighbour PMTs probability: " << data.neighbours_prob[0]<< G4endl;
	str << "To opposite PMT prob: " << data.opposite_prob[0]<< G4endl;
	str << "To the same PMT prob: " << data.same_prob[0]<< G4endl;
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

CustomRunManager::~CustomRunManager()
{
	//if (primary_Monte_Carlo) delete primary_Monte_Carlo;
	if (curr_mapping_state) delete curr_mapping_state;
	if (sim_results) delete sim_results;
	if (ev_history) delete ev_history;
};

void CustomRunManager::ProcessOneEvent(G4int i_event)
{
	OnEventStartProc();
	G4RunManager::ProcessOneEvent(i_event);
}

void CustomRunManager::TerminateOneEvent()
{
	if (sim_results)
		sim_results->RunEnd();
	if (ev_history)
		ev_history->RunEnd();
	G4RunManager::TerminateOneEvent();
}

void CustomRunManager::DoEventLoop(G4int n_event,const char* macroFile,G4int n_select)
{
  InitializeEventLoop(n_event,macroFile,n_select);
// Event loop
  G4int interval = n_event / 100.0;
  G4int counter = 0;
  for(G4int i_event=0; i_event<n_event; i_event++ )
  {
    ProcessOneEvent(i_event);
	if (interval*counter <= i_event)
	{
		counter++;
		G4cout << "simulated " << i_event << "/" << n_event << " ("<< G4int(100.0*i_event/n_event)<< "%)" << G4endl;
	}
    if(runAborted) break;
	TerminateOneEvent();
  }
  g_out:
  TerminateEventLoop();
}

void CustomRunManager::SetHit(G4int is_hit, G4double probab,const G4Step* step)
{
	if (sim_results)
		sim_results->SetHit(is_hit, probab,step);
	if (ev_history)
		ev_history->SetHit(is_hit, probab,step);
}

//0 - kill process, 1 - select reflection, 2 - select defraction, 3 - select both
G4int CustomRunManager::select_photon_BP(const G4Step* step, G4ThreeVector defl_momentum, G4ThreeVector refl_momentum)
{
	//TODO: major rework required (main point - cut off total reflection? depending on absortion length for a given energy?)
#ifdef DEBUG_MC_NODES
	G4cout << "RM: select_photon_BP()" << G4endl;
#endif
	const B1DetectorConstruction* detC
		= static_cast<const B1DetectorConstruction*> (this->GetUserDetectorConstruction());
	G4LogicalVolume* pre_volume=step->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume();
	G4LogicalVolume* post_volume = step->GetPostStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume();
	//if ((post_volume == detC->Xn_PMT) || (post_volume == detC->Xp_PMT) || (post_volume == detC->Yn_PMT) || (post_volume == detC->Yp_PMT))
	//	return RM_CHOOSE_DEFL; //TODO: should be both?
//#if defined(TOP_MESH_TEST)||defined(TEST_MESH_SIDEWAYS)
//	if (((pre_volume == detC->top_cell_hole) || (pre_volume == detC->top_cell_container) || (pre_volume == detC->top_cell_hole_dielectric)))
//	{
//		return RM_CHOOSE_BOTH;
//	}
//#else
	//6 unreachable volumes (absorbe light, but for speed increase are also here):
	//TODO: figure out dielectic's (g10's) properties
	/*if (((post_volume == detC->bot_cu_plate) || (post_volume == detC->top_cu_plate)) && (pre_volume!=post_volume))
			return RM_CHOOSE_REFL;
	if ((post_volume == detC->top_cell) && (post_volume!=pre_volume))
		return RM_CHOOSE_REFL;
	if ((post_volume == detC->top_cell_hole_dielectric) && (post_volume != pre_volume))
		return RM_CHOOSE_REFL;
#endif
	if ((post_volume == detC->bot_cell) && (post_volume != pre_volume))
		return RM_CHOOSE_REFL;
	if ((post_volume == detC->bot_cell_hole_dielectric) && (post_volume != pre_volume))
		return RM_CHOOSE_REFL;*/
	//no rules for holes, containers and pseudoplates, because by geometry defs they all have the same refractive index
	if (((pre_volume == detC->Xn_wls) || (pre_volume == detC->Xp_wls) || (pre_volume == detC->Yn_wls) || (pre_volume == detC->Yp_wls))
		&& ((post_volume == detC->box_interior) || (post_volume == detC->LAr_layer) || (post_volume == detC->box) ||
		(post_volume == detC->bot_ps_plate) || (post_volume == detC->top_ps_plate)))
		return 0;
	//^the last condition means that photons, leaving wls and heading inside are supreessed
	//^TODO: this probably should be done by adding absorbtion inside the argon

	return 1;
}

G4ThreeVector CustomRunManager::FetchPosition()
{return initial_position;}

G4ThreeVector CustomRunManager::GenPosition()
{
#ifdef DEBUG_MC_NODES
	//G4cout << "" << G4endl;
#endif
	/*if (current_working_node != NULL)
		return initial_position = current_working_node->GenPostion();*/
#ifdef TOP_MESH_TEST
	if (t_counter < x_num*y_num)
	{
		G4double x = x_start+t_step*((int)(t_counter / y_num));
		G4double y = y_start + t_step*((int)(t_counter % y_num));
		x += G4RandGauss::shoot(G4Random::getTheEngine(), 0.0, t_uncert);
		y += G4RandGauss::shoot(G4Random::getTheEngine(),0.0,t_uncert);
		t_counter++;	
		return initial_position = G4ThreeVector(x, y, 0);
	}
	else
		return initial_position = G4ThreeVector(0, 0, 0);
#endif
#ifdef TEST_MESH_SIDEWAYS
	return initial_position = G4ThreeVector(-9, 0, 10.9); //CELL SIZE = 4.5
#endif
//#ifdef TEMP_CODE_
//	return initial_position = G4ThreeVector(0.55, 0.45*sqrt(3)/2, 10.9);
//#endif
#ifdef AR_SPEC_TEST
	return initial_position = G4ThreeVector(0, -0.45*sqrt(3)/2, 10.9);
#endif
//#ifdef TEMP_CODE_
//	return initial_position = G4ThreeVector((141+WLS_FILM_WIDTH)/2, 0, 0);
//#endif
	G4double phi = 2*CLHEP::pi*G4UniformRand();
	G4double D = CONVERSION_DIAM > GEM_SIZE_TOTAL ? GEM_SIZE_TOTAL : CONVERSION_DIAM;
	G4double r = 0.5*D*sqrt(G4UniformRand());
	G4double z = (10.9-(-6.9))*G4UniformRand() - 6.9;
	return initial_position = G4ThreeVector(r*cos(phi), r*sin(phi), z);
}

G4ThreeVector CustomRunManager::FetchMomentum()
{return initial_momentum_direction;}

G4ThreeVector CustomRunManager::GenMomentum()
{
#ifdef DEBUG_MC_NODES
	//G4cout << "" << G4endl;
#endif
	/*if (current_working_node != NULL)
		return initial_momentum_direction = current_working_node->GenMomentum();*/
#ifdef TOP_MESH_TEST
	return initial_momentum_direction = G4ThreeVector(0, 0, 1);
#endif
#ifdef TEST_MESH_SIDEWAYS
	//1*cell_size/(2*h+plate_width)<z<3*cell_size/(2*h+plate_width), h==z_ps_gem_boundary - z_initial 
	G4double dx = 3 * 4.5 / (0.2 + 0.6);
	return initial_momentum_direction = G4ThreeVector(dx*cos(CLHEP::pi / 3), dx*sin(CLHEP::pi / 3), 1).unit(); //CELL SIZE = 4.5
#endif
//#ifdef TEMP_CODE_
//	return initial_momentum_direction = G4ThreeVector(2, 0, 1).unit();
//#endif
#ifdef AR_SPEC_TEST
	return initial_momentum_direction = G4ThreeVector(0, 0, 1).unit();
#endif
//#ifdef TEMP_CODE_
//	return initial_momentum_direction = G4ThreeVector(0, 1, 0).unit();
//#endif
	G4double phi = CLHEP::twopi*G4UniformRand();
	G4double cos_theta = 2*(G4UniformRand())-1;
	return initial_momentum_direction = G4ThreeVector(sin(phi)*sqrt(1 - cos_theta*cos_theta), cos(phi)*sqrt(1 - cos_theta*cos_theta), cos_theta);
}

G4ThreeVector CustomRunManager::FetchPolarization()
{return initial_polarization;}

G4ThreeVector CustomRunManager::GenPolarization()
{
#ifdef DEBUG_MC_NODES
	//G4cout << "" << G4endl;
#endif
	/*if (current_working_node != NULL)
		return initial_polarization = current_working_node->GenPolarization();*/
	G4double phi = CLHEP::twopi*G4UniformRand();
	G4double cos_theta = 2 * (G4UniformRand()) - 1;
	return initial_polarization = G4ThreeVector(sin(phi)*sqrt(1 - cos_theta*cos_theta), cos(phi)*sqrt(1 - cos_theta*cos_theta), cos_theta);
}

G4double CustomRunManager::FetchEnergy()
{return initial_energy;}

G4double CustomRunManager::GenEnergy()
{
#ifdef DEBUG_MC_NODES
		//G4cout << "" << G4endl;
#endif
		/*if (current_working_node != NULL)
			return initial_energy = current_working_node->GenEnergy();*/
#ifdef AR_EMISSION_NITRO
		if (EnergySpectrum)
			return initial_energy = EnergySpectrum->Value(G4UniformRand());
		else
#endif
	{
		// 9.65*eV ==128nm
		// 3.9236*eV == 316nm
		// 3.6791*eV == 337nm
		// 3.4632*eV == 358nm
		// 3.2585*eV == 380.5nm
		// 3.0538*eV == 406nm
//#ifdef TEMP_CODE_
//		return initial_energy = 2.9*eV; //~415nm
//#endif
		return initial_energy = 9.65*eV; //primary event parameters are such for a while
	}
}

PseudoMeshData*	CustomRunManager::GenMappingState() //Position is assumed to be already updated before the call
{
	PseudoMeshData* temp=NULL;
	/*if (current_working_node != NULL)
		temp = current_working_node->GenMappingData();*/
	if (temp == NULL) //in case current_working_node!=NULL this is erroneous situation but not fatal
	{
		B1DetectorConstruction* detectorConstruction = (B1DetectorConstruction*)(GetUserDetectorConstruction());
		detectorConstruction->GetPseudoMeshByPoint(curr_mapping_state, FetchPosition(),FetchMomentum());
		if (curr_mapping_state->curr_mesh != NULL)
			curr_mapping_state->curr_mesh->GetDefaultMappingData(curr_mapping_state);
		else
			curr_mapping_state->SetDefauldInd();
	}
	else
		*curr_mapping_state = *temp;
	return curr_mapping_state;
}

#ifdef AR_EMISSION_NITRO
void CustomRunManager::GenEnergySpectrum()
{
	B1DetectorConstruction* detectorConstruction = (B1DetectorConstruction*)(GetUserDetectorConstruction());
	G4VPhysicalVolume *ph_v = detectorConstruction->GetPVolumeByPoint(FetchPosition(), FetchMomentum());
	G4Material* mat;
	G4MaterialPropertiesTable* aMaterialPropertiesTable;
	G4LogicalVolume *l_v = ph_v ? ph_v->GetLogicalVolume():NULL;
	mat=l_v? l_v->GetMaterial():NULL;
	aMaterialPropertiesTable = mat?mat->GetMaterialPropertiesTable():NULL;
	if (aMaterialPropertiesTable)
		EnergySpectrum = aMaterialPropertiesTable->GetProperty("WLS_ENERGY_SPECTRUM");
}
#endif

#ifdef TOP_MESH_TEST
//!!WARNING - horizontal positioni of test detector is assumed
void CustomRunManager::on_hit_proc(G4ThreeVector point, G4double prob, G4int bot_top) //called when test detector is hit with photon
{
	G4double net_prob = prob;
	if (net_prob < 0) return;
	if (bot_top == 0)
	{
		bot_hits_xs.push_back(point.x());
		bot_hits_ys.push_back(point.y());
		bot_hits_probs.push_back(net_prob);
	}
	if (bot_top == 1)
	{
		top_hits_xs.push_back(point.x());
		top_hits_ys.push_back(point.y());
		top_hits_probs.push_back(net_prob);
	}
}

void CustomRunManager::export_to_bmp(std::list<G4double>* hits_xs, std::list<G4double> *hits_ys, std::list<G4double> *hits_probs, G4String filename) 
//^first step is to map position of hits with arbitrary coordinates to bitmap. Then write to file 
{
	t_counter = 0; //nullifies history so "scan" can be done again
	G4double * bits = new G4double[x_num*y_num];
	for (int g = 0; g < x_num*y_num; g++)
		bits[g] = 0;
	G4double lx = t_step*x_num;
	G4double ly = t_step*y_num;
	while ((!hits_xs->empty()) && (!hits_ys->empty()) && (!hits_probs->empty()))
	{
		G4double x = hits_xs->back()-x_start;
		G4double y = hits_ys->back()-y_start;
		G4double prob = hits_probs->back();
		hits_xs->pop_back();
		hits_ys->pop_back();
		hits_probs->pop_back();
		if ((x<0) || (x>t_step*x_num))
			continue;
		if ((y<0) || (y>t_step*y_num))
			continue;
		G4int x_ind = (int)(x /t_step);
		G4int y_ind = (int)(y /t_step);
#define NUM_TO_CONSIDER 3
		//^should be odd
		G4double neighbours[NUM_TO_CONSIDER*NUM_TO_CONSIDER];
		//firstly - obtain squared lengths to the centers of bits in bmp, then exp distribution to the neighbours
		//then - renormalization and writing to global double array
		G4double sum=0;
		for (G4int h = 0; h < NUM_TO_CONSIDER*NUM_TO_CONSIDER; h++)
		{
			G4int te_x_ind = (h / NUM_TO_CONSIDER) - NUM_TO_CONSIDER/2;
			G4int te_y_ind = h % NUM_TO_CONSIDER - NUM_TO_CONSIDER / 2;
			neighbours[h] = (x - t_step*(x_ind+te_x_ind))*(x - t_step*(x_ind + te_x_ind)) 
				+ (y - t_step*(y_ind + te_y_ind))*(y - t_step*(y_ind + te_y_ind));
			//^lengths to centers of neighbour cells
			neighbours[h] = exp(-neighbours[h] / (t_step*t_step));
			if (((te_x_ind + x_ind) > 0) && ((te_x_ind + x_ind) < x_num) && ((te_y_ind + y_ind) > 0) && ((te_y_ind + y_ind) < y_num))
				sum += neighbours[h];
		}
		if (0 != sum)
		{
			sum = prob / sum;
			for (G4int h = 0; h < NUM_TO_CONSIDER*NUM_TO_CONSIDER; h++)
			{
				G4int te_x_ind = (h / NUM_TO_CONSIDER) - NUM_TO_CONSIDER / 2;
				G4int te_y_ind = h % NUM_TO_CONSIDER - NUM_TO_CONSIDER / 2;
				neighbours[h] *= sum;
				if (((te_x_ind + x_ind) >=0 ) && ((te_x_ind + x_ind) < x_num) && ((te_y_ind + y_ind) >= 0) && ((te_y_ind + y_ind) < y_num))
					bits[(te_x_ind + x_ind)*y_num + (te_y_ind + y_ind)] += neighbours[h];
			}
		}
#undef NUM_TO_CONSIDER
	}
	G4double max = -1;
	for (int g = 0; g < x_num*y_num; g++)
	{
		if (max < bits[g]) max = bits[g];
	}
	if (max <= 0)
	{
		delete[] bits;
		return;
	}
	unsigned char *_bits = new unsigned char[x_num*y_num];
	for (int g = 0; g < x_num*y_num; g++)
		_bits[g] = (bits[g] < 0 ? 0 : (int)(255.0*bits[g] / max));
	delete[] bits;
	//WORK WITH BMP
	unsigned int headers[13];
	G4int extra_bytes = (4 - ((x_num * 3 )% 4))%4;
	int paddle_size = (x_num * 3 + extra_bytes)*y_num;
	headers[0] = paddle_size + 54;
	headers[1] = 0;
	headers[2] = 54;
	headers[3] = 40;
	headers[4] = x_num;  
	headers[5] = y_num;
	headers[7] = 0;
	headers[8] = paddle_size;
	headers[9] = 0;
	headers[10] = 0;
	headers[11] = 0;
	headers[12] = 0;
	std::ofstream bmp;
	bmp.open(filename,std::ios_base::trunc|std::ios_base::binary);
	bmp << "BM";
	for (G4int h = 0; h <= 5; h++)
	{
		bmp << (char)(headers[h] & 0x000000ff);
		bmp << (char)((headers[h] & 0x0000ff00)>>8);
		bmp << (char)((headers[h] & 0x00ff0000)>>16);
		bmp << (char)((headers[h] & (unsigned int) 0xff000000)>>24);
	}
	bmp << (char)1 << (char)0 << (char)24 << (char)0;
	for (G4int h = 7; h <= 12; h++)
	{
		bmp << (char)(headers[h] & 0x000000ff);
		bmp << (char)((headers[h] & 0x0000ff00) >> 8);
		bmp << (char)((headers[h] & 0x00ff0000) >> 16);
		bmp << (char)((headers[h] & (unsigned int)0xff000000) >> 24);
	}
	for (G4int h = 0; h < y_num;h++)
	{
		for (G4int g = 0; g < x_num; g++)
			bmp << _bits[g*y_num + h] << _bits[g*y_num + h] << _bits[g*y_num + h];
		for (int g = 0; g < extra_bytes; g++)
			bmp << (char)0;
	}
	bmp.close();
	//END BMP OUTPUT
	delete[] _bits;
}
#endif

void CustomRunManager::OnEventStartProc()
{
	B1DetectorConstruction* detectorConstruction = (B1DetectorConstruction*)(GetUserDetectorConstruction());
	GenPosition();
	GenPolarization();
	GenMomentum();
	GenMappingState();
#ifdef AR_EMISSION_NITRO
	GenEnergySpectrum();//after GenPos. and GenMom.
#endif
	GenEnergy();
	detectorConstruction->OnEventStartProc(this);
}

G4ThreeVector CustomRunManager::MappingProc(const G4Track& track, const G4Step& aStep, G4TouchableHandle &fCurrentTouchableHandle)
{
	B1DetectorConstruction* detectorConstruction = (B1DetectorConstruction*)GetUserDetectorConstruction();
	return detectorConstruction->MappingProc(curr_mapping_state, track, aStep, fCurrentTouchableHandle);
}

void CustomRunManager::OnNewSimulationProc(void)
{
#ifndef AR_SPEC_TEST
	/*sim_results->UpdateProbabilities(primary_Monte_Carlo);
	sim_results->UpdateSpectra(primary_Monte_Carlo);
	primary_Monte_Carlo->events.pop_back();*/ //discard previous simulation in order to decrease memory usage
#endif
}

void CustomRunManager::OnRunEndProc()
{
	if (sim_results)
		sim_results->OnEndSimulation();
}
