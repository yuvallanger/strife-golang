package sequential

import (
	"fmt"
	"math/rand"
	"miscow"
	"time"
)

var s4strain = [4]int{0, 0, 1, 1} // index is the strain. value is the signal allele.
var r4strain = [4]int{0, 1, 0, 1} // index is the strain. value is the receptor allele.

func strain_spec(r, s int) int {
	return r + s*2
}

func fitness(model *Model, coord Coordinate) float64 {
	// computes the fitness of a given cell
	var cost float64 = model.Parameters.Basal_Cost

	if model.Get_Cell_Prod(coord) {
		cost += model.Parameters.Cooperation_Cost
	}

	if model.Parameters.Cooperation_Effect_Threshold <= model.Get_Cell_PG_Num(coord) {
		cost *= model.Parameters.Public_Goods_Effect
	}

	return model.Parameters.Basal_Cost / cost

}

func mutate(model *Model, coord Coordinate) int {
	var r, s int
	var strain int = model.Get_Cell_Strain(coord)

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
	var (
		neighbor_strain          int
		signal_i, signal_i_torus Coordinate
		pg_i, pg_i_torus         Coordinate
		oldprod                  bool
	)

	for signal_i.r = coord.r - model.Parameters.S_Radius; signal_i.r <= coord.r+model.Parameters.S_Radius; signal_i.r++ {
		for signal_i.c = coord.c - model.Parameters.S_Radius; signal_i.c <= coord.c+model.Parameters.S_Radius; signal_i.c++ {
			// Get the toroid coordinates 
			signal_i_torus = signal_i.Get_Toroid_Coordinates(model.Parameters.Board_Size)
			// update signal level in sig range
			model.Add_To_Cell_Signal_Num(signal_i_torus, s4strain[oldstrain], -1)
			model.Add_To_Cell_Signal_Num(signal_i_torus, s4strain[newstrain], 1)
			// update producer status in sig range
			oldprod = model.Get_Cell_Prod(signal_i_torus)
			neighbor_strain = model.Get_Cell_Strain(signal_i_torus)
			if model.Parameters.Signal_Threshold <= model.Get_Cell_Signal_Num(signal_i_torus, r4strain[neighbor_strain]) {
				model.Set_Cell_Prod(signal_i_torus, true)
			} else {
				model.Set_Cell_Prod(signal_i_torus, false)
			}
			// update pg level in sig+pg range
			if oldprod != model.Get_Cell_Prod(signal_i_torus) {
				for pg_i.r = signal_i.r - model.Parameters.PG_Radius; pg_i.r <= signal_i.r+model.Parameters.PG_Radius; pg_i.r++ {
					for pg_i.c = signal_i.c - model.Parameters.PG_Radius; pg_i.c <= signal_i.c+model.Parameters.PG_Radius; pg_i.c++ {
						pg_i_torus = Coordinate{r: (pg_i.r + model.Parameters.Board_Size) % model.Parameters.Board_Size,
							c: (pg_i.c + model.Parameters.Board_Size) % model.Parameters.Board_Size}
						switch model.Get_Cell_Prod(signal_i_torus) {
						case true:
							model.Add_To_Cell_PG_Num(pg_i_torus, 1)
						default:
							model.Add_To_Cell_PG_Num(pg_i_torus, -1)
						}
					}
				}
			}
		}
	}
}

func endgame(model *Model, winner_coord, loser_coord Coordinate) {
	/*
	   wr - winner_coord's row
	   wc - winner_coord's column
	   lr - loser_coord's row
	   lc - loser_coord's column
	*/
	newstrain := mutate(model, winner_coord)
	if newstrain != model.Get_Cell_Strain(loser_coord) {
		oldstrain := model.Get_Cell_Strain(loser_coord)
		model.Set_Cell_Strain(loser_coord, newstrain)
		update_arrays(model, loser_coord, newstrain, oldstrain)
	}
}

func rand_neighbor(coord Coordinate, board_size int) (coord2 Coordinate) {
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

/*
Takes a model and turns four cells 90 degrees
*/
func diffuse(model *Model, coord00 Coordinate, direction bool) {
	// at each whirl, 4 cells are moved
	//var r0, c0, r1, c1 int
	var before, after [2][2]int

	// We get the coordinates for the four cells we'll turn.'
	coord11 := Coordinate{r: miscow.MyMod(coord00.r+1, model.Parameters.Board_Size),
		c: miscow.MyMod(coord00.c+1, model.Parameters.Board_Size)}
	coord01 := coord00
	coord01.c = coord11.c
	coord10 := coord00
	coord10.c = coord11.r

	// Save the tetrades
	before[0][0] = model.Get_Cell_Strain(coord00)
	before[0][1] = model.Get_Cell_Strain(coord01)
	before[1][0] = model.Get_Cell_Strain(coord10)
	before[1][1] = model.Get_Cell_Strain(coord11)

	if direction {
		// true is clockwise
		after[0][0] = before[1][0]
		after[0][1] = before[0][0]
		after[1][0] = before[1][1]
		after[1][1] = before[0][1]
	} else {
		// false is anticlockwise
		after[0][0] = before[0][1]
		after[0][1] = before[1][1]
		after[1][0] = before[0][0]
		after[1][1] = before[1][0]
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
	model.Set_Cell_Strain(coord00, after[0][0])
	model.Set_Cell_Strain(coord01, after[0][1])
	model.Set_Cell_Strain(coord10, after[1][0])
	model.Set_Cell_Strain(coord11, after[1][1])

	//fmt.Println("(*model.Board_strain)[r0][c0]", (*model.Board_strain)[r0][c0])
	//fmt.Println("(*model.Board_strain)[r0][c1]", (*model.Board_strain)[r0][c1])
	//fmt.Println("(*model.Board_strain)[r1][c0]", (*model.Board_strain)[r1][c0])
	//fmt.Println("(*model.Board_strain)[r1][c1]", (*model.Board_strain)[r1][c1])

	update_arrays(model, coord00, after[0][0], before[0][0])
	update_arrays(model, coord01, after[0][1], before[0][1])
	update_arrays(model, coord10, after[1][0], before[1][0])
	update_arrays(model, coord11, after[1][1], before[1][1])
}

func showtiming(t_start time.Time, dt_iter time.Duration) {
	t_elapsed := time.Now().Sub(t_start)
	dt_tot_runtime := time.Duration(dt_iter.Nanoseconds()*10000) * time.Nanosecond
	t_finish := t_start.Add(dt_tot_runtime)

	fmt.Println("Since start:", t_elapsed)
	fmt.Println("Expected total run time:", dt_tot_runtime)
	fmt.Println("Finish time:", t_finish)
}
