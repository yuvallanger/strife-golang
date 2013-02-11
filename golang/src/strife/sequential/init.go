package sequential

import (
	"math/rand"
	"miscow"
)

func init_board_strain(model *Model) {
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

func init_board_signal_num(model *Model) {
	// init on the metal
	model.Board_signal_num = make([][][]int, 2)
	for strain_i := range model.Board_signal_num {
		model.Board_signal_num[strain_i] = make([][]int, model.Parameters.Board_Size)
		for row_i := range model.Board_signal_num[strain_i] {
			model.Board_signal_num[strain_i][row_i] = make([]int, model.Parameters.Board_Size)
		}
	}

	// init in the model
	center_i := Coordinate{}
	for center_i.r = 0; center_i.r < model.Parameters.Board_Size; center_i.r++ {
		for center_i.c = 0; center_i.c < model.Parameters.Board_Size; center_i.c++ {
			rad_i := Coordinate{}
			for rad_i.r = center_i.r - model.Parameters.S_Radius; rad_i.r <= center_i.r+model.Parameters.S_Radius; rad_i.r++ {
				for rad_i.c = center_i.c - model.Parameters.S_Radius; rad_i.c <= center_i.c+model.Parameters.S_Radius; rad_i.c++ {
					// here we count the number of signals for each cell

					// get strain at rad_i
					strain_at_rad := model.Get_Cell_Strain(rad_i.Get_Toroid_Coordinates(model.Parameters.Board_Size))

					// the allele of signal of the strain at rad_i
					signal_strain_at_rad := s4strain[strain_at_rad]

					// add one signal at center_i of the signal allele from rad_i
					model.Add_To_Cell_Signal_Num(center_i, signal_strain_at_rad, 1)
				}
			}
		}
	}
}

func init_board_prod(model *Model) Board_prod {
	// init on the metal
	var board_prod Board_prod
	board_prod = make([][]bool, model.Parameters.Board_Size)
	for i0 := range board_prod {
		board_prod[i0] = make([]bool, model.Parameters.Board_Size)
	}

	// init in the model
	for center_row_i := range model.Board_prod {
		for center_col_i := range (model.Board_prod)[center_row_i] {
			current_strain := (model.Board_strain)[center_row_i][center_col_i]
			current_receptor_strain := r4strain[current_strain]
			if (model.Board_signal_num)[current_receptor_strain][center_row_i][center_col_i] >= model.Parameters.Signal_Threshold {
				board_prod[center_row_i][center_col_i] = true
			}
		}
	}
	return board_prod
}

func init_board_pg_num(model *Model) Board_pg_num {
	var rad_row_i_t, rad_col_i_t int
	var board_pg_num Board_pg_num
	board_pg_num = make([][]int, model.Parameters.Board_Size)
	for row_i := range board_pg_num {
		board_pg_num[row_i] = make([]int, model.Parameters.Board_Size)
	}

	for center_row_i := range board_pg_num {
		for center_col_i := range board_pg_num[center_row_i] {
			for rad_row_i := center_row_i - model.Parameters.PG_Radius; rad_row_i < center_row_i+model.Parameters.PG_Radius+1; rad_row_i++ {
				rad_row_i_t = miscow.MyMod(rad_row_i, model.Parameters.Board_Size)
				for rad_col_i := center_col_i - model.Parameters.PG_Radius; rad_col_i < center_col_i+model.Parameters.PG_Radius+1; rad_col_i++ {
					rad_col_i_t = miscow.MyMod(rad_col_i, model.Parameters.Board_Size)
					if (model.Board_prod)[rad_row_i_t][rad_col_i_t] {
						board_pg_num[center_row_i][center_col_i]++
					}
				}
			}
		}
	}
	return board_pg_num
}

func init_boards(model *Model) {
	init_board_strain(model)
	init_board_signal_num(model)
	init_board_prod(model)
	init_board_pg_num(model)
}
