package sequential

import ()

// We need only one init function for both models
func (model *Model) InitBoardSignalNum() {
	// init on the metal
	model.BoardSignalNum = make([][][]int, 2)
	for strain_i := range model.BoardSignalNum {
		model.BoardSignalNum[strain_i] = make([][]int, model.Parameters.BoardSize)
		for row_i := range model.BoardSignalNum[strain_i] {
			model.BoardSignalNum[strain_i][row_i] = make([]int, model.Parameters.BoardSize)
		}
	}

	// init in the model
	center_coord := Coordinate{}
	for center_coord.r = 0; center_coord.r < model.Parameters.BoardSize; center_coord.r++ {
		for center_coord.c = 0; center_coord.c < model.Parameters.BoardSize; center_coord.c++ {
			rad_coord := Coordinate{}
			for rad_coord.r = center_coord.r - model.Parameters.SRadius; rad_coord.r <= center_coord.r+model.Parameters.SRadius; rad_coord.r++ {
				for rad_coord.c = center_coord.c - model.Parameters.SRadius; rad_coord.c <= center_coord.c+model.Parameters.SRadius; rad_coord.c++ {
					// here we count the number of signals for each cell

					// get strain at rad_coord
					strain_at_rad := model.CellStrain(rad_coord.ToroidCoordinates(model.Parameters.BoardSize))

					// the allele of signal of the strain at rad_coord
					signal_strain_at_rad := s4strain[strain_at_rad]

					// add one signal at center_coord of the signal allele from rad_coord
					model.AddToCellSignalNum(center_coord, signal_strain_at_rad, 1)
				}
			}
		}
	}
}

// TODO this comment sucks. initialize all the data samples
func (model *Model) initDataSamples() {
	model.initDataSamplesSnapshots()
	model.initDataSamplesFrequencies()
	model.initDataSamplesNeighborhoodFrequencies()
}

// Initialize the board snapshots samples
func (model *Model) initDataSamplesSnapshots() {
	if model.Settings.SnapshotsSampleNum != 0 {
		if model.Parameters.Generations%model.Settings.SnapshotsSampleNum == 0 {
			// the case in which the last snapshot is the same as the last generation
			model.DataSamples.Snapshots = make([]Snapshot, model.Settings.SnapshotsSampleNum)
		} else {
			// the case in which the last snapshot isn't the same as the last generation
			model.DataSamples.Snapshots = make([]Snapshot, model.Settings.SnapshotsSampleNum+1)
		}
		for sample_i := range model.DataSamples.Snapshots {
			model.DataSamples.Snapshots[sample_i].Data = make([][]int, model.Parameters.BoardSize)
			for row := range model.DataSamples.Snapshots[sample_i].Data {
				model.DataSamples.Snapshots[sample_i].Data[row] = make([]int, model.Parameters.BoardSize)
			}
			//fmt.Println(sample_i, model.DataSamples.Snapshots[sample_i])
		}
		model.DataSamples.Snapshots = model.DataSamples.Snapshots[0:1]
	}
}

// Initialize the strain frequencies samples
func (model *Model) initDataSamplesFrequencies() {
	if model.Settings.FrequenciesSampleNum != 0 {
		if model.Parameters.Generations%model.Settings.FrequenciesSampleNum == 0 {
			// the case in which the last sample is the same as the last generation
			model.DataSamples.Frequencies = make([]Frequency, model.Settings.FrequenciesSampleNum)
		} else {
			// the case in which the last sample isn't the same as the last generation
			model.DataSamples.Frequencies = make([]Frequency, model.Settings.FrequenciesSampleNum+1)
		}
		for sample_i := range model.DataSamples.Frequencies {
			model.DataSamples.Frequencies[sample_i].Data = make([]int, 8)
		}
		model.DataSamples.Frequencies = model.DataSamples.Frequencies[0:1]
	}
}

// Initialize the neighbors frequencies samples
func (model *Model) initDataSamplesNeighborhoodFrequencies() {
	if model.Settings.FrequenciesSampleNum != 0 {
		if model.Parameters.Generations%model.Settings.FrequenciesSampleNum == 0 {
			// the case in which the last sample is the same as the last generation
			model.DataSamples.NeighborsFrequencies = make([]NeighborsFrequency, model.Settings.FrequenciesSampleNum)
		} else {
			// the case in which the last sample isn't the same as the last generation
			model.DataSamples.NeighborsFrequencies = make([]NeighborsFrequency, model.Settings.FrequenciesSampleNum+1)
		}
		for sample_i := range model.DataSamples.NeighborsFrequencies {
			data := make([][]int, 8)
			// for each strain we'll count how many strains are around it.
			for strain_i := range data {
				data[strain_i] = make([]int, 8)
			}
			model.DataSamples.NeighborsFrequencies[sample_i].Data = data
		}
		model.DataSamples.NeighborsFrequencies = model.DataSamples.NeighborsFrequencies[0:1]
	}
}
