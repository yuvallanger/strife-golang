package sequential

import (
	"math/rand"
	"miscow"
)

func (model *CzaranModel) fitness(coord Coordinate) float64 {
	// computes the fitness of a given cell
	var cost float64 = model.Parameters.BasalCost

	coordToroid := coord.toroidCoordinates(model.Parameters.BoardSize)

	if model.cellProd(coordToroid) {
		cost += model.Parameters.CooperationCost
	}

	if model.Parameters.CooperationEffectThreshold <= model.cellPGNum(coordToroid) {
		cost *= model.Parameters.PublicGoodsEffect
	}

	return model.Parameters.BasalCost / cost

}

func (model *CzaranModel) mutate(coord Coordinate) int {
	var r, s int
	var strain int = model.getCellStrain(coord)

	r = r4strain[strain]
	if rand.Float64() < model.Parameters.MutOddsR {
		r = 1 - r
	}

	s = s4strain[strain]
	if rand.Float64() < model.Parameters.MutOddsS {
		s = 1 - s
	}

	return strainSpec(r, s, 1) // in Avigdor's model, the public goods allele is always functional (but not always expressed).
}

func (model *CzaranModel) updateArrays(coord Coordinate, newstrain, oldstrain int) {
	var (
		neighborStrain                int
		signalCoord, signalCoordTorus Coordinate
		PGCoord, PGCoordTorus         Coordinate
		oldProd                       bool
	)

	for signalCoord.r = coord.r - model.Parameters.SRadius; signalCoord.r <= coord.r+model.Parameters.SRadius; signalCoord.r++ {
		for signalCoord.c = coord.c - model.Parameters.SRadius; signalCoord.c <= coord.c+model.Parameters.SRadius; signalCoord.c++ {
			signalCoordTorus = signalCoord.toroidCoordinates(model.Parameters.BoardSize)

			// update signal level in sig range
			model.addToCellSignalNum(signalCoordTorus, s4strain[oldstrain], -1)
			model.addToCellSignalNum(signalCoordTorus, s4strain[newstrain], 1)

			oldProd = model.cellProd(signalCoordTorus)
			neighborStrain = model.getCellStrain(signalCoordTorus)
			// update producer status at signalCoordTorus
			// Production is active if cell has a working PG gene and:
			//     The receptor allele is malfunctioned. OR
			//     The receptor allele is working and the signal level is above threshold.
			if g4strain[neighborStrain] == 1 {
				if r4strain[neighborStrain] == 1 {
					if model.Parameters.SignalThreshold <= model.cellSignalNum(signalCoordTorus, 1) {
						model.setCellProd(signalCoordTorus, true)
					} else {
						model.setCellProd(signalCoordTorus, false)
					}
				} else {
					model.setCellProd(signalCoordTorus, true)
				}
			} else {
				model.setCellProd(signalCoordTorus, false)
			}

			// update PG levels around signalCoord
			if oldProd != model.cellProd(signalCoordTorus) {
				for PGCoord.r = signalCoord.r - model.Parameters.PGRadius; PGCoord.r <= signalCoord.r+model.Parameters.PGRadius; PGCoord.r++ {
					for PGCoord.c = signalCoord.c - model.Parameters.PGRadius; PGCoord.c <= signalCoord.c+model.Parameters.PGRadius; PGCoord.c++ {
						PGCoordTorus = PGCoord.toroidCoordinates(model.Parameters.BoardSize)

						// We change PG level at signalCoord's neighbors
						if model.cellProd(signalCoordTorus) {
							model.addToCellPGNum(PGCoordTorus, 1)
						} else {
							model.addToCellPGNum(PGCoordTorus, -1)
						}
					}
				}
			}
		}
	}
}

/*
Endgame copies the winning cell into the losing one's position.
Before the cell is copied, we pass it through mutate().
*/
func (model *CzaranModel) endgame(winnerCoord, loserCoord Coordinate) {
	newstrain := model.mutate(winnerCoord)
	oldstrain := model.getCellStrain(loserCoord)
	if newstrain != oldstrain {
		model.setCellStrain(loserCoord, newstrain)
		model.updateArrays(loserCoord, newstrain, oldstrain)
	}
}

func (model *CzaranModel) competition() {

	var c1, c2 Coordinate
	c1 = randCoord(model.Parameters.BoardSize)
	c2 = randNeighbor(c1, model.Parameters.BoardSize)

	fitness_1 := model.fitness(c1)
	fitness_2 := model.fitness(c2)

	// Randomize fo shizzles
	score := rand.Float64()*fitness_1 - rand.Float64()*fitness_2

	if score > 0 {
		// cell 1 wins
		model.endgame(c1, c2)
	} else {
		// cell 2 wins
		model.endgame(c2, c1)
	}
}

/*
Takes a model and turns four cells 90 degrees
*/
func (model *CzaranModel) diffuse(coord00 Coordinate, direction bool) {
	// at each whirl, 4 cells are moved
	//var r0, c0, r1, c1 int
	var before, after [2][2]int

	// We get the coordinates for the four cells we'll turn.'
	coord11 := Coordinate{r: miscow.MyMod(coord00.r+1, model.Parameters.BoardSize),
		c: miscow.MyMod(coord00.c+1, model.Parameters.BoardSize)}
	coord01 := coord00
	coord01.c = coord11.c
	coord10 := coord00
	coord10.c = coord11.r

	// Save the tetrades
	before[0][0] = model.getCellStrain(coord00)
	before[0][1] = model.getCellStrain(coord01)
	before[1][0] = model.getCellStrain(coord10)
	before[1][1] = model.getCellStrain(coord11)

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

	// Assign the rotated cells
	model.setCellStrain(coord00, after[0][0])
	model.setCellStrain(coord01, after[0][1])
	model.setCellStrain(coord10, after[1][0])
	model.setCellStrain(coord11, after[1][1])

	model.updateArrays(coord00, after[0][0], before[0][0])
	model.updateArrays(coord01, after[0][1], before[0][1])
	model.updateArrays(coord10, after[1][0], before[1][0])
	model.updateArrays(coord11, after[1][1], before[1][1])
}
