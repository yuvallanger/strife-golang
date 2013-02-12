package sequential

import (
	"math/rand"
	"miscow"
)

var s4strain = [4]int{0, 0, 1, 1} // index is the strain. value is the signal allele.
var r4strain = [4]int{0, 1, 0, 1} // index is the strain. value is the receptor allele.

func strain_spec(r, s int) int {
	return r + s*2
}

func (model *Model) fitness(coord Coordinate) float64 {
	// computes the fitness of a given cell
	var cost float64 = model.Parameters.Basal_Cost

	coord_toroid := coord.get_toroid_coordinates(model.Parameters.Board_Size)

	if model.get_cell_prod(coord_toroid) {
		cost += model.Parameters.Cooperation_Cost
	}

	if model.Parameters.Cooperation_Effect_Threshold <= model.get_cell_pg_num(coord_toroid) {
		cost *= model.Parameters.Public_Goods_Effect
	}

	return model.Parameters.Basal_Cost / cost

}

func mutate(model *Model, coord Coordinate) int {
	var r, s int
	var strain int = model.get_cell_strain(coord)

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
		neighbor_strain                  int
		signal_coord, signal_coord_torus Coordinate
		pg_coord, pg_coord_torus         Coordinate
		oldprod                          bool
	)

	for signal_coord.r = coord.r - model.Parameters.S_Radius; signal_coord.r <= coord.r+model.Parameters.S_Radius; signal_coord.r++ {
		for signal_coord.c = coord.c - model.Parameters.S_Radius; signal_coord.c <= coord.c+model.Parameters.S_Radius; signal_coord.c++ {
			signal_coord_torus = signal_coord.get_toroid_coordinates(model.Parameters.Board_Size)

			// update signal level in sig range
			model.Add_To_Cell_Signal_Num(signal_coord_torus, s4strain[oldstrain], -1)
			model.Add_To_Cell_Signal_Num(signal_coord_torus, s4strain[newstrain], 1)

			// update producer status at signal_coord_torus
			oldprod = model.get_cell_prod(signal_coord_torus)
			neighbor_strain = model.get_cell_strain(signal_coord_torus)
			if model.Parameters.Signal_Threshold <= model.get_cell_signal_num(signal_coord_torus, r4strain[neighbor_strain]) {
				model.set_cell_prod(signal_coord_torus, true)
			} else {
				model.set_cell_prod(signal_coord_torus, false)
			}

			// update PG levels around signal_coord
			if oldprod != model.get_cell_prod(signal_coord_torus) {
				for pg_coord.r = signal_coord.r - model.Parameters.PG_Radius; pg_coord.r <= signal_coord.r+model.Parameters.PG_Radius; pg_coord.r++ {
					for pg_coord.c = signal_coord.c - model.Parameters.PG_Radius; pg_coord.c <= signal_coord.c+model.Parameters.PG_Radius; pg_coord.c++ {
						pg_coord_torus = pg_coord.get_toroid_coordinates(model.Parameters.Board_Size)
						if model.get_cell_prod(signal_coord_torus) {
							model.Add_To_Cell_PG_Num(pg_coord_torus, 1)
						} else {
							model.Add_To_Cell_PG_Num(pg_coord_torus, -1)
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
	if newstrain != model.get_cell_strain(loser_coord) {
		oldstrain := model.get_cell_strain(loser_coord)
		model.set_cell_strain(loser_coord, newstrain)
		update_arrays(model, loser_coord, newstrain, oldstrain)
	}
}

func rand_neighbor(coord Coordinate, board_size int) (coord2 Coordinate) {
	switch direction := rand.Intn(4); direction {
	case 0:
		return Coordinate{r: coord.r,
			c: miscow.MyMod(coord.c-1, board_size)}
	case 1:
		return Coordinate{r: miscow.MyMod(coord.r-1, board_size),
			c: coord.c}
	case 2:
		return Coordinate{r: coord.r,
			c: miscow.MyMod(coord.c+1, board_size)}
	default:
		return Coordinate{r: miscow.MyMod(coord.r+1, board_size),
			c: coord.c}
	}
	return
}

func (model *Model) competition() {

	var c1, c2 Coordinate
	c1 = rand_coord(model.Parameters.Board_Size)
	c2 = rand_neighbor(c1, model.Parameters.Board_Size)

	fitness_1 := model.fitness(c1)
	fitness_2 := model.fitness(c2)

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
func (model *Model) diffuse(coord00 Coordinate, direction bool) {
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
	before[0][0] = model.get_cell_strain(coord00)
	before[0][1] = model.get_cell_strain(coord01)
	before[1][0] = model.get_cell_strain(coord10)
	before[1][1] = model.get_cell_strain(coord11)

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
	model.set_cell_strain(coord00, after[0][0])
	model.set_cell_strain(coord01, after[0][1])
	model.set_cell_strain(coord10, after[1][0])
	model.set_cell_strain(coord11, after[1][1])

	//fmt.Println("(*model.Board_strain)[r0][c0]", (*model.Board_strain)[r0][c0])
	//fmt.Println("(*model.Board_strain)[r0][c1]", (*model.Board_strain)[r0][c1])
	//fmt.Println("(*model.Board_strain)[r1][c0]", (*model.Board_strain)[r1][c0])
	//fmt.Println("(*model.Board_strain)[r1][c1]", (*model.Board_strain)[r1][c1])

	update_arrays(model, coord00, after[0][0], before[0][0])
	update_arrays(model, coord01, after[0][1], before[0][1])
	update_arrays(model, coord10, after[1][0], before[1][0])
	update_arrays(model, coord11, after[1][1], before[1][1])
}
