package sequential

import (
	"math/rand"
	"gitlab.com/yuvallanger/miscow"
)

func StrainSpec(r, s, g int) int {
	return r + s*2 + g*4
}

func RandNeighbor(coord Coordinate, board_size int) (coord2 Coordinate) {
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

// Rotate 2x2 cells, but do not update the boards. It is left for each model's Diffuse() function
func (model *Model) Rotate90(coord00 Coordinate, direction bool) (coordinates [2][2]Coordinate, before, after [2][2]int) {
	// We'll set the coordinates for the four cells we're rotating.
	coordinates[0][0] = coord00
	coordinates[1][1] = Coordinate{r: coordinates[0][0].r + 1, c: coordinates[0][0].c + 1}
	coordinates[1][1] = coordinates[1][1].ToroidCoordinates(model.Parameters.BoardSize)
	coordinates[0][1] = Coordinate{r: coordinates[0][0].r, c: coordinates[1][1].c}
	coordinates[1][0] = Coordinate{r: coordinates[1][1].r, c: coordinates[0][0].c}

	// Save the tetrades
	before[0][0] = model.CellStrain(coordinates[0][0])
	before[0][1] = model.CellStrain(coordinates[0][1])
	before[1][1] = model.CellStrain(coordinates[1][1])
	before[1][0] = model.CellStrain(coordinates[1][0])

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
	model.SetCellStrain(coordinates[0][0], after[0][0])
	model.SetCellStrain(coordinates[0][1], after[0][1])
	model.SetCellStrain(coordinates[1][0], after[1][0])
	model.SetCellStrain(coordinates[1][1], after[1][1])

    return
}
