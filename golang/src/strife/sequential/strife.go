package sequential

import (
	"fmt"
	"github.com/yuvallanger/game-of-strife/golang/src/miscow"
	"math/rand"
	"time"
)

var s4strain = []int{0, 0, 1, 1}
var r4strain = []int{0, 1, 0, 1}

func strain_spec(r, s int) int {
	miscow.Trace("strain_spec")
	defer miscow.Untrace("strain_spec")
	return r + s*2
}

func fitness(model *Model, coord Coordinate) float64 {
	miscow.Trace("fitness")
	defer miscow.Untrace("fitness")
	// computes the fitness of a given cell
	var cost float64 = model.Parameters.Basal_Cost

	if model.Board_prod.get_cell(coord) == true {
		cost += model.Parameters.Cooperation_Cost
	}

	if model.Parameters.Cooperation_Effect_Threshold <= model.Board_pg_num.get_cell(coord) {
		cost *= model.Parameters.Public_Goods_Effect
	}

	return model.Parameters.Basal_Cost / cost

}

func mutate(model *Model, coord Coordinate) int {
	miscow.Trace("mutate")
	defer miscow.Untrace("mutate")
	var r, s int
	var strain int = model.Board_strain.get_cell(coord)

	r = r4strain[strain]
	if rand.Float64() < model.Parameters.Mut_Odds_R {
		r = 1 - r
	}

	s = s4strain[strain]
	if rand.Float64() < model.Parameters.Mut_Odds_S {
		s = 1 - s
	}

	return (r + 2*s)
}

func update_arrays(model *Model, coord Coordinate, newstrain, oldstrain int) {
	miscow.Trace("update_arrays")
	defer miscow.Untrace("update_arrays")
	var (
		s_row_i, s_col_i, s_row_i_t, s_col_i_t           int
		pg_row_i, pg_col_i, pg_row_i_t, pg_col_i_t       int
		s_rad, pg_rad, s_th, board_size, neighbor_strain int
	)
	var oldprod bool
	(func(int, int, int, int, int, int, int, int, int, int, int, int, int, bool) {})(
		s_row_i, s_col_i, s_row_i_t, s_col_i_t,
		pg_row_i, pg_col_i, pg_row_i_t, pg_col_i_t,
		s_rad, pg_rad, s_th, board_size,
		neighbor_strain,
		oldprod)
	s_rad = model.Parameters.S_Radius
	pg_rad = model.Parameters.PG_Radius
	s_th = model.Parameters.Signal_Threshold
	board_size = model.Parameters.Board_Size
	// TODO convert s_rad, c_rad, s_th, board_size back to model.Parameters...

	for s_row_i = coord.r - s_rad; s_row_i <= coord.r+s_rad; s_row_i++ {
		for s_col_i = coord.c - s_rad; s_col_i <= coord.c+s_rad; s_col_i++ {
			s_row_i_t = miscow.MyMod(s_row_i, board_size)
			s_col_i_t = miscow.MyMod(s_col_i, board_size)
			// update signal level in sig range
			(*model.Board_signal_num)[s4strain[oldstrain]][s_row_i_t][s_col_i_t]--
			(*model.Board_signal_num)[s4strain[newstrain]][s_row_i_t][s_col_i_t]++
			// update producer status in sig range
			oldprod = model.GetCellProd(Coordinate{s_row_i_t, s_col_i_t})
			neighbor_strain = model.GetCellStrain(Coordinate{s_row_i_t, s_col_i_t})
			if s_th <= model.GetCellSignalNum(r4strain[neighbor_strain], Coordinate{s_row_i_t, s_col_i_t}) {
				model.SetCellProd(Coordinate{s_row_i_t, s_col_i_t}, true)
			} else {
				model.SetCellProd(Coordinate{s_row_i_t, s_col_i_t}, false)
			}
			// update pg level in sig+pg range
			if oldprod != (*model.Board_prod)[s_row_i_t][s_col_i_t] {
				for pg_row_i = s_row_i_t - pg_rad; pg_row_i <= s_row_i_t+pg_rad; pg_row_i++ {
					for pg_col_i = s_col_i_t - pg_rad; pg_col_i <= s_col_i_t+pg_rad; pg_col_i++ {
						pg_row_i_t = (pg_row_i + board_size) % board_size
						pg_col_i_t = (pg_col_i + board_size) % board_size
						switch (*model.Board_prod)[s_row_i_t][s_col_i_t] {
						case true:
							(*model.Board_pg_num)[pg_row_i_t][pg_col_i_t] += 1
						default:
							(*model.Board_pg_num)[pg_row_i_t][pg_col_i_t] -= 1

						}
					}
				}
			}
		}
	}
}

func endgame(model *Model, winner, loser Coordinate) {
	miscow.Trace("endgame")
	defer miscow.Untrace("endgame")
	/*
	   wr - winner's row
	   wc - winner's column
	   lr - loser's row
	   lc - loser's column
	*/
	newstrain := mutate(model, winner)
	if newstrain != model.Board_strain.get_cell(loser) {
		oldstrain := model.Board_strain.get_cell(loser)
		(model.Board_strain).set_cell(loser, newstrain)
		update_arrays(model, loser, newstrain, oldstrain)
	}
}

func rand_neighbor(coord Coordinate, board_size int) (coord2 Coordinate) {
	miscow.Trace("rand_neighbor")
	defer miscow.Untrace("rand_neighbor")
	switch direction := rand.Intn(4); direction {
	case 0:
		coord2.r = coord.r
		coord2.c = miscow.MyMod(coord.c-1, board_size)
	case 1:
		coord2.r = miscow.MyMod(coord.r-1, board_size)
		coord2.c = coord.r
	case 2:
		coord2.r = coord.r
		coord2.c = miscow.MyMod(coord.c+1, board_size)
	case 3:
		coord2.r = miscow.MyMod(coord.r+1, board_size)
		coord2.c = coord.r
	}
	return
}

func competition(model *Model) {
	miscow.Trace("competition")
	defer miscow.Untrace("competition")
	var c1, c2 Coordinate
	c1 = rand_coord(model.Parameters.Board_Size)
	c2 = rand_neighbor(c1, model.Parameters.Board_Size)

	fitness_1 := fitness(model, c1)
	fitness_2 := fitness(model, c2)

	// Randomize fo shizzles
	score := rand.Float64()*fitness_1 - rand.Float64()*fitness_2

	if score > 0 {
		// cell 1 wins
		endgame(model, c1, c2)
	} else {
		// cell 2 wins
		endgame(model, c2, c1)
	}
}

func diffuse(model *Model) {
	miscow.Trace("diffuse")
	defer miscow.Untrace("diffuse")
	// at each whirl, 4 cells are moved
	//var r0, c0, r1, c1 int
	var before, after [2][2]int
	// We get the coordinates for the cell tetrade
	r0 := rand.Intn(model.Parameters.Board_Size)
	c0 := rand.Intn(model.Parameters.Board_Size)
	r1 := miscow.MyMod(r0+1, model.Parameters.Board_Size)
	c1 := miscow.MyMod(c0+1, model.Parameters.Board_Size)

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

	//fmt.Println("[r0][c0]", r0, c0)
	//fmt.Println("[r0][c1]", r0, c1)
	//fmt.Println("[r1][c0]", r1, c0)
	//fmt.Println("[r1][c1]", r1, c1)
	//fmt.Println("before[0][0]", before[0][0])
	//fmt.Println("before[0][1]", before[0][1])
	//fmt.Println("before[1][0]", before[1][0])
	//fmt.Println("before[1][1]", before[1][1])
	//fmt.Println("after[0][0]", after[0][0])
	//fmt.Println("after[0][1]", after[0][1])
	//fmt.Println("after[1][0]", after[1][0])
	//fmt.Println("after[1][1]", after[1][1])
	//fmt.Println("(*model.Board_strain)[r0][c0]", (*model.Board_strain)[r0][c0])
	//fmt.Println("(*model.Board_strain)[r0][c1]", (*model.Board_strain)[r0][c1])
	//fmt.Println("(*model.Board_strain)[r1][c0]", (*model.Board_strain)[r1][c0])
	//fmt.Println("(*model.Board_strain)[r1][c1]", (*model.Board_strain)[r1][c1])

	// Assign the rotated cells
	(*model.Board_strain)[r0][c0] = after[0][0]
	(*model.Board_strain)[r0][c1] = after[0][1]
	(*model.Board_strain)[r1][c0] = after[1][0]
	(*model.Board_strain)[r1][c1] = after[1][1]

	(func(int, int, int, int, [2][2]int, [2][2]int) {})(r0, c0, r1, c1, before, after)
	//fmt.Println("(*model.Board_strain)[r0][c0]", (*model.Board_strain)[r0][c0])
	//fmt.Println("(*model.Board_strain)[r0][c1]", (*model.Board_strain)[r0][c1])
	//fmt.Println("(*model.Board_strain)[r1][c0]", (*model.Board_strain)[r1][c0])
	//fmt.Println("(*model.Board_strain)[r1][c1]", (*model.Board_strain)[r1][c1])

	update_arrays(model, Coordinate{r: r0, c: c0}, after[0][0], before[0][0])
	update_arrays(model, Coordinate{r: r0, c: c1}, after[0][1], before[0][1])
	update_arrays(model, Coordinate{r: r0, c: c0}, after[1][0], before[1][0])
	update_arrays(model, Coordinate{r: r1, c: c1}, after[1][1], before[1][1])
}

func showtiming(t_start time.Time, dt_iter time.Duration) {
	miscow.Trace("showtiming")
	defer miscow.Untrace("showtiming")
	t_elapsed := time.Now().Sub(t_start)
	dt_tot_runtime := time.Duration(dt_iter.Nanoseconds()*10000) * time.Nanosecond
	t_finish := t_start.Add(dt_tot_runtime)

	fmt.Println("Since start:", t_elapsed)
	fmt.Println("Expected total run time:", dt_tot_runtime)
	fmt.Println("Finish time:", t_finish)
}

func run(model *Model) {
	miscow.Trace("strife.go: run()")
	defer miscow.Untrace("strife.go: run()")
	var competition_i, diffusion_i int
	var t_iter_start time.Time
	var diffusion_num int = int(model.Parameters.D * float64(model.Parameters.Board_Size*model.Parameters.Board_Size) / 4)
	fmt.Println("diffusion_num: ", diffusion_num)
	fmt.Println("model.Parameters.Generations = ", model.Parameters.Generations)
	fmt.Println("model.Generation_i = ", model.Generation_i)
	tstart := time.Now()
	for model.Generation_i = 0; model.Generation_i < model.Parameters.Generations; model.Generation_i++ {
		// TODO
		/*
			if model.Generation_i%100 == 0 {
				     model.Data_Boards.Snapshots[model.Generation_i / 10] = struct {Data: model.Strain , Generation: model.Generation_i}
			}
		*/
		t_iter_start = time.Now()

		board_size := model.Parameters.Board_Size
		fmt.Println("model.Generation_i = ", model.Generation_i)
		for competition_i = 0; competition_i < board_size*board_size; competition_i++ {
			competition(model)
		}
		for diffusion_i = 0; diffusion_i < diffusion_num; diffusion_i++ {
			diffuse(model)
		}

		showtiming(tstart, time.Since(t_iter_start))

		//func (t Time) Add(d Duration) Time
		//func Since(t Time) Duration
		//func (t Time) Sub(u Time) Duration

		fmt.Println()
	}
	return
}
