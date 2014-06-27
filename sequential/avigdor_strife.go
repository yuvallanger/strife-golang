package sequential

import (
	"math/rand"
)

// Computes the Avigdor fitness of a cell at a given Coordinate
func (model *Model) Fitness(coord Coordinate) float64 {
	var cost float64 = model.Parameters.BasalCost

	coordToroid := coord.ToroidCoordinates(model.Parameters.BoardSize)

	if model.CellProd(coordToroid) {
		cost += model.Parameters.CooperationCost
	}

	cost += model.Parameters.SignalCost
	cost += model.Parameters.ReceptorCost

	if model.Parameters.CooperationEffectThreshold <= model.CellPGNum(coordToroid) {
		cost *= (1 - model.Parameters.PublicGoodsEffect)
	}

	return model.Parameters.BasalCost / cost
}

// Mutate the provided strain
func (model *Model) Mutate(strain int) int {
	r := r4strain[strain]
	if rand.Float64() < model.Parameters.MutOddsR {
		r = 1 - r
	}

	s := s4strain[strain]
	if rand.Float64() < model.Parameters.MutOddsS {
		s = 1 - s
	}

	// in Avigdor's model, the public goods allele is always functional (but not always expressed).
	return StrainSpec(r, s, 1)
}

func (model *Model) UpdateArrays(coord Coordinate, newstrain, oldstrain int) {
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
			// Production is active if the signal level at signal_coord (that is also compatible with signal_coord's receptor) is above the quorum sensing threshold.
			if model.Parameters.SignalThreshold <= model.CellSignalNum(signalCoordTorus, r4strain[neighborStrain]) {
				model.SetCellProd(signalCoordTorus, true)
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
func (model *Model) Endgame(winnerCoord, loserCoord Coordinate) {
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
func (model *Model) Competition() {

	var c1, c2 Coordinate
	c1 = RandCoord(model.Parameters.BoardSize)
	c2 = RandNeighbor(c1, model.Parameters.BoardSize)

	if model.CellStrain(c1) == model.CellStrain(c2) {
		// doesn't matter who wins if both are of the same strain, so we pick c1 as winner and c2 as loser.
		model.Endgame(c1, c2)
	} else {
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
}

/*
Rotate four (2x2) cells 90 degrees.

Cells rotate clockwise if "direction" is true and anticlockwise if "direction" is false.
*/
func (model *Model) Diffuse(coord00 Coordinate, direction bool) {
	coordinates, before, after := model.Rotate90(coord00, direction)

	model.UpdateArrays(coordinates[0][0], after[0][0], before[0][0])
	model.UpdateArrays(coordinates[0][1], after[0][1], before[0][1])
	model.UpdateArrays(coordinates[1][1], after[1][1], before[1][1])
	model.UpdateArrays(coordinates[1][0], after[1][0], before[1][0])
}
