package main

import (
	"code.google.com/p/gcfg"
	"fmt"
	"math/rand"
	"miscow"
)

var s4strain = []int{0, 0, 1, 1}
var r4strain = []int{0, 1, 0, 1}

func strain_spec(r, s int) int {
	return r + s*2
}

func myMod(a, b int) int {
	return ((a % b) + b) % b
}

func load_config() (parameters *Parameters_T, settings *Settings_T) {
	var cfg Config
	err := gcfg.ReadFileInto(&cfg, "strife.conf")
	if err != nil {
		fmt.Println(err)
		parameters = &Parameters
		settings = &Settings
	} else {
		parameters = &cfg.Parameters
		settings = &cfg.Settings
	}
	return
}

func score(row, col int, p float64) float64 {
	return 1
}

func mutate(model *Model, row, col int) int {
	var r, s int
	var strain int = (*model.Board_strain)[row][col]

	r = r4strain[strain]
	if rand.Float64() < model.Parameters.MutOddsR {
		r = 1 - r
	}

	s = s4strain[strain]
	if rand.Float64() < model.Parameters.MutOddsS {
		s = 1 - s
	}

	return (r + 2*s)
}

func update_arrays(model *Model, row, col, newstrain, oldstrain int) {
	var (
		s_row_i, s_col_i, s_row_i_t, s_col_i_t int
		c_row_i, c_col_i, c_row_i_t, c_col_i_t int
		s_rad, c_rad, s_th, board_size         int
		oldprod                                int
	)
	s_rad = model.Parameters.SRad
	c_rad = model.Parameters.CRad
	s_th = model.Parameters.STh
	board_size = model.Parameters.BoardSize

	for s_row_i = row - s_rad; s_row_i <= col+s_rad; s_row_i++ {
		for s_col_i = col - s_rad; s_col_i <= col+s_rad; s_col_i++ {
			s_row_i_t = myMod(s_row_i, board_size)
			s_col_i_t = myMod(s_col_i, board_size)
			// update signal level in sig range
			(*model.Board_signal_num)[s4strain[oldstrain]][s_row_i_t][s_col_i_t]--
			(*model.Board_signal_num)[s4strain[newstrain]][s_row_i_t][s_col_i_t]++
			// update producer status in sig range
			oldprod = (*model.Board_prod)[s_row_i_t][s_col_i_t]
			neighbor_strain := (*model.Board_strain)[s_row_i_t][s_col_i_t]
			if s_th <= (*model.Board_signal_num)[r4strain[neighbor_strain]][s_row_i_t][s_col_i_t] {
				(*model.Board_prod)[s_row_i_t][s_col_i_t] = 1
			} else {
				(*model.Board_prod)[s_row_i_t][s_col_i_t] = 0
			}
			// update pg level in sig+pg range
			if oldprod != (*model.Board_prod)[s_row_i_t][s_col_i_t] {
				for c_row_i = c_row_i_t - c_rad; c_row_i <= c_row_i_t+c_rad; c_row_i++ {
					for c_col_i = c_col_i_t - c_rad; c_col_i <= c_col_i_t+c_rad; c_col_i++ {
						c_row_i_t = (c_row_i + board_size) % board_size
						c_col_i_t = (c_col_i + board_size) % board_size
						(*model.Board_pg_num)[c_row_i_t][c_col_i_t] += ((*model.Board_prod)[s_row_i_t][s_col_i_t] - oldprod)
					}
				}
			}
		}
	}
}

func endgame(model *Model, wr, wc, lr, lc int) {
	/*
	   wr - winner's row
	   wc - winner's column
	   lr - loser's row
	   lc - loser's column
	*/
	newstrain := mutate(model, wr, wc)
	if newstrain != (*model.Board_strain)[lr][lc] {
		oldstrain := (*model.Board_strain)[lr][lc]
		(*model.Board_strain)[lr][lc] = newstrain
		update_arrays(model, lr, lc, newstrain, oldstrain)
	}
}

func competition(model *Model) {
	var c1r, c1c, c2r, c2c int
	switch direction := rand.Intn(4); direction {
	case 0:
		c2r = c1r
		c2c = myMod(c1c-1, model.Parameters.BoardSize)
	case 1:
		c2r = myMod(c1r-1, model.Parameters.BoardSize)
		c2c = c2r
	case 2:
		c2r = c1r
		c2c = myMod(c1c+1, model.Parameters.BoardSize)
	case 3:
		c2r = myMod(c1r+1, model.Parameters.BoardSize)
		c2c = c2r
	}
	score_1 := score(c1r, c1c, rand.Float64())
	score_2 := score(c2r, c2c, rand.Float64())
	if score_1 > score_2 {
		endgame(model, c1r, c1c, c2r, c2c)
	} else {
		endgame(model, c1r, c1c, c2r, c2c)
	}
}

func diffuse(model *Model) {
	// at each whirl, 4 cells are moved
	var r0, c0, r1, c1 int
	var before, after [2][2]int
	diffusion_num := model.Parameters.D * model.Parameters.BoardSize / 4
	for diffusion_i := 0; diffusion_i < diffusion_num; diffusion_i++ {
		// We get the coordinates for the cell tetrade
		r0 = rand.Intn(model.Parameters.BoardSize)
		c0 = rand.Intn(model.Parameters.BoardSize)
		r1 = myMod(r0, model.Parameters.BoardSize)
		c1 = myMod(c0, model.Parameters.BoardSize)

		// Save the tetrades
		before[0][0] = (*model.Board_strain)[r0][c0]
		before[0][1] = (*model.Board_strain)[r0][c1]
		before[1][0] = (*model.Board_strain)[r1][c0]
		before[1][1] = (*model.Board_strain)[r1][c1]

		switch rand.Intn(2) {
		case 0:
			// 0 is anticlockwise
			after[0][0] = before[0][1]
			after[0][1] = before[1][1]
			after[1][0] = before[0][0]
			after[1][1] = before[1][0]
		case 1:
			// 1 is clockwise
			after[0][0] = before[1][0]
			after[0][1] = before[0][0]
			after[1][0] = before[1][1]
			after[1][1] = before[0][1]
		}

		// Assign the rotated cells
		(*model.Board_strain)[r0][c0] = after[0][0]
		(*model.Board_strain)[r0][c1] = after[0][1]
		(*model.Board_strain)[r1][c0] = after[1][0]
		(*model.Board_strain)[r1][c1] = after[1][1]

		update_arrays(model, r0, c0, after[0][0], before[0][0])
		update_arrays(model, r0, c1, after[0][1], before[0][1])
		update_arrays(model, r0, c0, after[1][0], before[1][0])
		update_arrays(model, r1, c1, after[1][1], before[1][1])
	}
}

func run(model *Model) {
	miscow.Trace("strife.go: run()")
	defer miscow.Trace("strife.go: run()")
	for ; model.Generation_i < model.Parameters.Generations; model.Generation_i++ {
		fmt.Println(i)
		diffuse(model)
	}
	return
}
