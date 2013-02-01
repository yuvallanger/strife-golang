package main

import (
	"math/rand"
	"miscow"
)

func init_board_strain(model *Model) *Board_strain {
	miscow.Trace("init_board_strain")
	defer miscow.Untrace("init_board_strain")
	// init on the metal
	var board_strain Board_strain
	board_strain = make([][]int, model.Parameters.BoardSize)
	for row_i := range board_strain {
		board_strain[row_i] = make([]int, model.Parameters.BoardSize)
	}

	// init in the model
	for row_i := range board_strain {
		for col_i := range board_strain[row_i] {
			if rand.Float64() < model.Parameters.RInitOdds {
				board_strain[row_i][col_i] += 1
			}
			if rand.Float64() < model.Parameters.SInitOdds {
				board_strain[row_i][col_i] += 2
			}

		}
	}
	return &board_strain
}

func init_board_signal_num(model *Model) *Board_signal_num {
	miscow.Trace("init_board_signal_num")
	defer miscow.Untrace("init_board_signal_num")
	// init on the metal
	var board_signal_num Board_signal_num
	board_signal_num = make([][][]int, 2)
	for strain_i := range board_signal_num {
		board_signal_num[strain_i] = make([][]int, model.Parameters.BoardSize)
		for row_i := range board_signal_num[strain_i] {
			board_signal_num[strain_i][row_i] = make([]int, model.Parameters.BoardSize)
		}
	}

	// init in the model
	for center_row_i := range *model.Board_strain {
		for center_col_i := range (*model.Board_strain)[center_row_i] {
			for rad_row_i := center_row_i - model.Parameters.SRad; rad_row_i < center_row_i+model.Parameters.SRad+1; rad_row_i++ {
				for rad_col_i := center_col_i - model.Parameters.SRad; rad_col_i < center_col_i+model.Parameters.SRad+1; rad_col_i++ {
					current_strain := (*model.Board_strain)[myMod(rad_row_i, model.Parameters.BoardSize)][myMod(rad_col_i, model.Parameters.BoardSize)]
					current_signal_strain := s4strain[current_strain]
					board_signal_num[current_signal_strain][center_row_i][center_col_i] = board_signal_num[current_signal_strain][center_row_i][center_col_i] + 1
					//fmt.Println(center_row_i, center_col_i, rad_row_i, rad_col_i, current_strain, current_signal_strain, board_signal_num[current_signal_strain][center_row_i][center_col_i])
				}
			}
		}
	}
	return &board_signal_num
}

func init_board_prod(model *Model) *Board_prod {
	miscow.Trace("init_board_prod")
	defer miscow.Untrace("init_board_prod")
	// init on the metal
	var board_prod Board_prod
	board_prod = make([][]int, model.Parameters.BoardSize)
	for i0 := range board_prod {
		board_prod[i0] = make([]int, model.Parameters.BoardSize)
	}

	// init in the model
	for center_row_i := range *model.Board_signal_num {
		for center_col_i := range (*model.Board_signal_num)[center_row_i] {
			current_strain := (*model.Board_strain)[center_row_i][center_col_i]
			current_receptor_strain := r4strain[current_strain]
			if (*model.Board_signal_num)[current_receptor_strain][center_row_i][center_col_i] >= model.Parameters.STh {
				board_prod[center_row_i][center_col_i] = 1
			}
		}
	}
	return &board_prod
}

func init_board_pg_num(model *Model) *Board_pg_num {
	miscow.Trace("init_board_pg_num")
	defer miscow.Untrace("init_board_pg_num")
	var board_pg_num Board_pg_num
	board_pg_num = make([][]int, model.Parameters.BoardSize)
	for row_i := range board_pg_num {
		board_pg_num[row_i] = make([]int, model.Parameters.BoardSize)
	}

	for center_row_i := range board_pg_num {
		for center_col_i := range board_pg_num[center_row_i] {
			for rad_row_i := center_row_i - model.Parameters.CRad; rad_row_i < center_row_i+model.Parameters.CRad+1; rad_row_i++ {
				for rad_col_i := center_col_i - model.Parameters.CRad; rad_col_i < center_col_i+model.Parameters.CRad+1; rad_col_i++ {
					if (*model.Board_prod)[myMod(rad_row_i, model.Parameters.BoardSize)][myMod(rad_col_i, model.Parameters.BoardSize)] >= model.Parameters.CTh {
						board_pg_num[center_row_i][center_col_i]++
					}
				}
			}
		}
	}
	return &board_pg_num
}

func init_boards(model *Model) {
	miscow.Trace("init_boards")
	defer miscow.Untrace("init_boards")
	model.Board_strain = init_board_strain(model)
	model.Board_signal_num = init_board_signal_num(model)
	model.Board_prod = init_board_prod(model)
	model.Board_pg_num = init_board_pg_num(model)
}
