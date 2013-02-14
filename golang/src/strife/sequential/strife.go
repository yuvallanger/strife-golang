package sequential

import (
	"math/rand"
	"miscow"
)

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
