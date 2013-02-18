package sequential

import (
	"math/rand"
)

// Compute the Czaran fitness of the cell at the specified Coordinate
func (model *CzaranModel) Fitness(coord Coordinate) float64 {
	var cost float64 = model.Parameters.BasalCost

	coordToroid := coord.ToroidCoordinates(model.Parameters.BoardSize)

	strain := model.CellStrain(coordToroid)

	// If the cell produces receptor, add cost of receptor production
	if 1 == r4strain[strain] {
		cost += model.Parameters.ReceptorCost
	}

	// If the cell produces signal, add cost of signal production
	if 1 == s4strain[strain] {
		cost += model.Parameters.SignalCost
	}

	// If there's public goods production, add cost of cooperation
	if model.CellProd(coordToroid) {
		cost += model.Parameters.CooperationCost
	}

	if model.Parameters.CooperationEffectThreshold <= model.CellPGNum(coordToroid) {
		cost *= model.Parameters.PublicGoodsEffect
	}

	return model.Parameters.BasalCost / cost
}

// Mutate the provided strain
func (model *CzaranModel) Mutate(strain int) int {
	r := r4strain[strain]
	if rand.Float64() < model.Parameters.MutOddsR {
		r = 1 - r
	}

	s := s4strain[strain]
	if rand.Float64() < model.Parameters.MutOddsS {
		s = 1 - s
	}

	g := g4strain[strain]
	if rand.Float64() < model.Parameters.MutOddsG {
		g = 1 - g
	}

	return StrainSpec(r, s, g)
}

func (model *CzaranModel) UpdateArrays(coord Coordinate, newstrain, oldstrain int) {
	var (
		neighborStrain                int
		signalCoord, signalCoordTorus Coordinate
		PGCoord, PGCoordTorus         Coordinate
		oldProd                       bool
	)

	for signalCoord.r = coord.r - model.Parameters.SRadius; signalCoord.r <= coord.r+model.Parameters.SRadius; signalCoord.r++ {
		for signalCoord.c = coord.c - model.Parameters.SRadius; signalCoord.c <= coord.c+model.Parameters.SRadius; signalCoord.c++ {
			signalCoordTorus = signalCoord.ToroidCoordinates(model.Parameters.BoardSize)

			// update signal level in sig range
			model.AddToCellSignalNum(signalCoordTorus, s4strain[oldstrain], -1)
			model.AddToCellSignalNum(signalCoordTorus, s4strain[newstrain], 1)

			oldProd = model.CellProd(signalCoordTorus)
			neighborStrain = model.CellStrain(signalCoordTorus)
			// update producer status at signalCoordTorus
			// Production is active if cell has a working PG gene and:
			//     The receptor allele is malfunctioned. OR
			//     The receptor allele is working and the signal level is above threshold.
			if g4strain[neighborStrain] == 1 {
				if r4strain[neighborStrain] == 1 {
					if model.Parameters.SignalThreshold <= model.CellSignalNum(signalCoordTorus, 1) {
						model.SetCellProd(signalCoordTorus, true)
					} else {
						model.SetCellProd(signalCoordTorus, false)
					}
				} else {
					model.SetCellProd(signalCoordTorus, true)
				}
			} else {
				model.SetCellProd(signalCoordTorus, false)
			}

			// update PG levels around signalCoord
			if oldProd != model.CellProd(signalCoordTorus) {
				for PGCoord.r = signalCoord.r - model.Parameters.PGRadius; PGCoord.r <= signalCoord.r+model.Parameters.PGRadius; PGCoord.r++ {
					for PGCoord.c = signalCoord.c - model.Parameters.PGRadius; PGCoord.c <= signalCoord.c+model.Parameters.PGRadius; PGCoord.c++ {
						PGCoordTorus = PGCoord.ToroidCoordinates(model.Parameters.BoardSize)

						// We change PG level at signalCoord's neighbors
						if model.CellProd(signalCoordTorus) {
							model.AddToCellPGNum(PGCoordTorus, 1)
						} else {
							model.AddToCellPGNum(PGCoordTorus, -1)
						}
					}
				}
			}
		}
	}
}

/*
Endgame() copies the winning cell into the losing one's position.
Before the winning cell is copied, we mutate it using Mutate().
*/
func (model *CzaranModel) Endgame(winnerCoord, loserCoord Coordinate) {
	newstrain := model.Mutate(model.CellStrain(winnerCoord))
	oldstrain := model.CellStrain(loserCoord)
	if newstrain != oldstrain {
		model.SetCellStrain(loserCoord, newstrain)
		model.UpdateArrays(loserCoord, newstrain, oldstrain)
	}
}

/*
Competition() takes two adjacent cells, decides who wins, copies the winner,
    mutates it and replaces the losing cell with the mutated copy.
*/
func (model *CzaranModel) Competition() {

	var c1, c2 Coordinate
	c1 = RandCoord(model.Parameters.BoardSize)
	c2 = RandNeighbor(c1, model.Parameters.BoardSize)

	fitness_1 := model.Fitness(c1)
	fitness_2 := model.Fitness(c2)

	// Randomize fo shizzles
	score := rand.Float64()*fitness_1 - rand.Float64()*fitness_2

	if score > 0 {
		// cell 1 wins
		model.Endgame(c1, c2)
	} else {
		// cell 2 wins
		model.Endgame(c2, c1)
	}
}

/*
Rotate four (2x2) cells 90 degrees.

Cells rotate clockwise if "direction" is true and anticlockwise if "direction" is false.
*/
func (model *CzaranModel) Diffuse(coord00 Coordinate, direction bool) {
	var before, after [2][2]int

	// We get the coordinates for the four cells we'll rotate.
	coord11 := Coordinate{r: coord00.r + 1,
		c: coord00.c + 1}
	coord11 = coord11.ToroidCoordinates(model.Parameters.BoardSize)
	coord01 := Coordinate{r: coord00.r,
		c: coord11.c}
	coord10 := Coordinate{r: coord11.r,
		c: coord00.c}

	// Save the tetrade of cells
	before[0][0] = model.CellStrain(coord00)
	before[0][1] = model.CellStrain(coord01)
	before[1][1] = model.CellStrain(coord11)
	before[1][0] = model.CellStrain(coord10)

	if direction {
		// true is clockwise
		after[0][0] = before[1][0]
		after[0][1] = before[0][0]
		after[1][1] = before[0][1]
		after[1][0] = before[1][1]
	} else {
		// false is anticlockwise
		after[0][0] = before[0][1]
		after[0][1] = before[1][1]
		after[1][1] = before[1][0]
		after[1][0] = before[0][0]
	}

	// Assign the rotated cells
	model.SetCellStrain(coord00, after[0][0])
	model.SetCellStrain(coord01, after[0][1])
	model.SetCellStrain(coord11, after[1][1])
	model.SetCellStrain(coord10, after[1][0])

	model.UpdateArrays(coord00, after[0][0], before[0][0])
	model.UpdateArrays(coord01, after[0][1], before[0][1])
	model.UpdateArrays(coord11, after[1][1], before[1][1])
	model.UpdateArrays(coord10, after[1][0], before[1][0])
}
