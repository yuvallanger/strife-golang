package sequential

var DefaultParameters = Parameters{
	Generations:                10,
	RInitOdds:                  0,
	SInitOdds:                  0,
	GInitOdds:                  0,
	SignalThreshold:            3,
	CooperationEffectThreshold: 3,
	SRadius:                    1,
	PGRadius:                   1,
	BoardSize:                  10,
	D:                          1,
	MutOddsR:                   1e-4,
	MutOddsS:                   1e-4,
	MutOddsG:                   1e-4,
	BasalCost:                  100,
	CooperationCost:            30,
	SignalCost:                 3,
	ReceptorCost:               1,
	PublicGoodsEffect:          0.6}

var DefaultSettings = Settings{
	GenerationsPerSnapshotSample:  50,
	GenerationsPerFrequencySample: 50,
	DataFilename:                  "model-data.dat"}
