package sequential

var Default_Parameters = Parameters{
	Generations:                  10,
	R_Init_Odds:                  0,
	S_Init_Odds:                  0,
	Signal_Threshold:             3,
	Cooperation_Effect_Threshold: 3,
	S_Radius:                     1,
	PG_Radius:                    1,
	Board_Size:                   10,
	D:                            1,
	Mut_Odds_R:                   1e-4,
	Mut_Odds_S:                   1e-4,
	Basal_Cost:                   100,
	Cooperation_Cost:             30,
	Signal_Cost:                  3,
	Receptor_Cost:                1,
	Public_Goods_Effect:          0.6}

var Default_Settings = Settings{
	Snapshots_sample_num:                5,
	Frequencies_sample_num:              10,
	Neighborhood_frequencies_sample_num: 15,
	Data_filename:                       "model-data.dat"}
