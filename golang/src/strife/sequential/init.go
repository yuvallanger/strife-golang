package sequential

import (
	"math/rand"
)

func (model *Model) init_board_strain() {
	// init on the metal
	model.Board_strain = make([][]int, model.Parameters.Board_Size)
	for row_i := range model.Board_strain {
		model.Board_strain[row_i] = make([]int, model.Parameters.Board_Size)
	}

	// init in the model
	for row_i := range model.Board_strain {
		for col_i := range model.Board_strain[row_i] {
			if rand.Float64() < model.Parameters.R_Init_Odds {
				model.Board_strain[row_i][col_i] += 1
			}
			if rand.Float64() < model.Parameters.S_Init_Odds {
				model.Board_strain[row_i][col_i] += 2
			}

		}
	}
}

func (model *Model) init_board_signal_num() {
	// init on the metal
	model.Board_signal_num = make([][][]int, 2)
	for strain_i := range model.Board_signal_num {
		model.Board_signal_num[strain_i] = make([][]int, model.Parameters.Board_Size)
		for row_i := range model.Board_signal_num[strain_i] {
			model.Board_signal_num[strain_i][row_i] = make([]int, model.Parameters.Board_Size)
		}
	}

	// init in the model
	center_coord := Coordinate{}
	for center_coord.r = 0; center_coord.r < model.Parameters.Board_Size; center_coord.r++ {
		for center_coord.c = 0; center_coord.c < model.Parameters.Board_Size; center_coord.c++ {
			rad_i := Coordinate{}
			for rad_i.r = center_coord.r - model.Parameters.S_Radius; rad_i.r <= center_coord.r+model.Parameters.S_Radius; rad_i.r++ {
				for rad_i.c = center_coord.c - model.Parameters.S_Radius; rad_i.c <= center_coord.c+model.Parameters.S_Radius; rad_i.c++ {
					// here we count the number of signals for each cell

					// get strain at rad_i
					strain_at_rad := model.get_cell_strain(rad_i.get_toroid_coordinates(model.Parameters.Board_Size))

					// the allele of signal of the strain at rad_i
					signal_strain_at_rad := s4strain[strain_at_rad]

					// add one signal at center_coord of the signal allele from rad_i
					model.Add_To_Cell_Signal_Num(center_coord, signal_strain_at_rad, 1)
				}
			}
		}
	}
}

func (model *Model) init_board_prod() {
	// init on the metal
	model.Board_prod = make([][]bool, model.Parameters.Board_Size)
	for i0 := range model.Board_prod {
		model.Board_prod[i0] = make([]bool, model.Parameters.Board_Size)
	}

	// init in the model
	center_coord := Coordinate{}
	for center_coord.r = range model.Board_prod {
		for center_coord.c = range model.Board_prod[center_coord.r] {
			strain_at_center_coord := model.get_cell_strain(center_coord)
			receptor_allele_at_center_coord := r4strain[strain_at_center_coord]
			if model.get_cell_signal_num(center_coord, receptor_allele_at_center_coord) >= model.Parameters.Signal_Threshold {
				model.set_cell_prod(center_coord, 1 > 0)
			}
		}
	}
}

func (model *Model) init_board_pg_num() {
	model.Board_pg_num = make([][]int, model.Parameters.Board_Size)
	for row_i := range model.Board_pg_num {
		model.Board_pg_num[row_i] = make([]int, model.Parameters.Board_Size)
	}

	center_coord := Coordinate{}
	for center_coord.r = range model.Board_pg_num {
		for center_coord.c = range model.Board_pg_num[center_coord.r] {
			rad_i := Coordinate{}
			for rad_i.r = center_coord.r - model.Parameters.PG_Radius; rad_i.r < center_coord.r+model.Parameters.PG_Radius+1; rad_i.r++ {
				for rad_i.c = center_coord.c - model.Parameters.PG_Radius; rad_i.c < center_coord.c+model.Parameters.PG_Radius+1; rad_i.c++ {
					rad_i_t := rad_i.get_toroid_coordinates(model.Parameters.Board_Size)
					if model.get_cell_prod(rad_i_t) {
						model.Add_To_Cell_PG_Num(center_coord, 1)
					}
				}
			}
		}
	}
}

func (model *Model) init_boards() {
	model.init_board_strain()
	model.init_board_signal_num()
	model.init_board_prod()
	model.init_board_pg_num()
}

// TODO this comment sucks. initialize all the data samples
func (model *Model) init_data_samples() {
	model.init_data_samples_snapshots()
	model.init_data_samples_frequencies()
	model.init_data_samples_neighborhood_frequencies()
}

// Initialize the board snapshots samples
func (model *Model) init_data_samples_snapshots() {
	if model.Settings.Snapshots_sample_num != 0 {
		if model.Parameters.Generations%model.Settings.Snapshots_sample_num == 0 {
			// the case in which the last snapshot is the same as the last generation
			model.Data_samples.Snapshots = make([]Snapshot, 1, model.Settings.Snapshots_sample_num)
		} else {
			// the case in which the last snapshot isn't the same as the last generation
			model.Data_samples.Snapshots = make([]Snapshot, 1, model.Settings.Snapshots_sample_num+1)
		}
		for sample_i := range model.Data_samples.Snapshots {
			data := make([][]int, model.Parameters.Board_Size)
			for row := range data {
				data[row] = make([]int, model.Parameters.Board_Size)
			}
			model.Data_samples.Snapshots[sample_i].Data = data
		}
	}
}

// Initialize the strain frequencies samples
func (model *Model) init_data_samples_frequencies() {
	if model.Settings.Frequencies_sample_num != 0 {
		if model.Parameters.Generations%model.Settings.Frequencies_sample_num == 0 {
			// the case in which the last sample is the same as the last generation
			model.Data_samples.Frequencies = make([]Frequency, 1, model.Settings.Frequencies_sample_num)
		} else {
			// the case in which the last sample isn't the same as the last generation
			model.Data_samples.Frequencies = make([]Frequency, 1, model.Settings.Frequencies_sample_num+1)
		}
		for sample_i := range model.Data_samples.Frequencies {
			model.Data_samples.Frequencies[sample_i].Data = make([]int, 8)
		}
	}
}

// Initialize the neighbors frequencies samples
func (model *Model) init_data_samples_neighborhood_frequencies() {
	if model.Settings.Neighborhood_frequencies_sample_num != 0 {
		if model.Parameters.Generations%model.Settings.Neighborhood_frequencies_sample_num == 0 {
			// the case in which the last sample is the same as the last generation
			model.Data_samples.Neighbors_frequencies = make([]Neighbors_frequency, model.Settings.Neighborhood_frequencies_sample_num)
		} else {
			// the case in which the last sample isn't the same as the last generation
			model.Data_samples.Neighbors_frequencies = make([]Neighbors_frequency, model.Settings.Neighborhood_frequencies_sample_num+1)
		}
		for sample_i := range model.Data_samples.Neighbors_frequencies {
			data := make([][]int, 8)
			// for each strain we'll count how many strains are around it.
			for strain_i := range data {
				data[strain_i] = make([]int, 8)
			}
			model.Data_samples.Neighbors_frequencies[sample_i].Data = data
		}
	}
	model.Data_samples.Neighbors_frequencies = model.Data_samples.Neighbors_frequencies[0:1]
}
