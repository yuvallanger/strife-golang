package sequential

import (
	"math/rand"
)

func (model *Model) initBoardStrain() {
	// init on the metal
	model.BoardStrain = make([][]int, model.Parameters.BoardSize)
	for row_i := range model.BoardStrain {
		model.BoardStrain[row_i] = make([]int, model.Parameters.BoardSize)
	}

	// init in the model
	for row_i := range model.BoardStrain {
		for col_i := range model.BoardStrain[row_i] {
			if rand.Float64() < model.Parameters.RInitOdds {
				model.BoardStrain[row_i][col_i] += 1
			}
			if rand.Float64() < model.Parameters.SInitOdds {
				model.BoardStrain[row_i][col_i] += 2
			}
			// In Avigdor's model, all public goods are wt.
			model.BoardStrain[row_i][col_i] += 4
		}
	}
}

func (model *Model) initBoardProd() {
	// init on the metal
	model.BoardProd = make([][]bool, model.Parameters.BoardSize)
	for i0 := range model.BoardProd {
		model.BoardProd[i0] = make([]bool, model.Parameters.BoardSize)
	}

	// init in the model
	centerCoord := Coordinate{}
	for centerCoord.r = range model.BoardProd {
		for centerCoord.c = range model.BoardProd[centerCoord.r] {
			centerStrain := model.CellStrain(centerCoord)
			centerReceptor := r4strain[centerStrain]
			if model.CellSignalNum(centerCoord, centerReceptor) >= model.Parameters.SignalThreshold {
				model.SetCellProd(centerCoord, true)
			}
		}
	}
}

func (model *Model) initBoardPGNum() {
	model.BoardPGNum = make([][]int, model.Parameters.BoardSize)
	for row_i := range model.BoardPGNum {
		model.BoardPGNum[row_i] = make([]int, model.Parameters.BoardSize)
	}

	center_coord := Coordinate{}
	for center_coord.r = range model.BoardPGNum {
		for center_coord.c = range model.BoardPGNum[center_coord.r] {
			rad_coord := Coordinate{}
			for rad_coord.r = center_coord.r - model.Parameters.PGRadius; rad_coord.r < center_coord.r+model.Parameters.PGRadius+1; rad_coord.r++ {
				for rad_coord.c = center_coord.c - model.Parameters.PGRadius; rad_coord.c < center_coord.c+model.Parameters.PGRadius+1; rad_coord.c++ {
					rad_coord_t := rad_coord.ToroidCoordinates(model.Parameters.BoardSize)
					if model.CellProd(rad_coord_t) {
						model.AddToCellPGNum(center_coord, 1)
					}
				}
			}
		}
	}
}

func (model *Model) initBoards() {
	model.initBoardStrain()
	model.InitBoardSignalNum()
	model.initBoardProd()
	model.initBoardPGNum()
}
