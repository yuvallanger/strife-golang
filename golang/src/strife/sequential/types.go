package sequential

import (
	"fmt"
	"math/rand"
	"miscow"
)

type Board_strain [][]int       // [rows][columns]. possible values: {0,1,2,3}. s0r0 - 0; s0r1 - 1; s1r0 - 2; s1r1 - 3
type Board_signal_num [][][]int // [signal types][rows][columns]
type Board_prod [][]bool        // [rows][columns]
type Board_pg_num [][]int       // [rows][columns]

type Model struct {
	Parameters       Parameters
	Settings         Settings
	Board_strain     *Board_strain
	Board_signal_num *Board_signal_num
	Board_prod       *Board_prod
	Board_pg_num     *Board_pg_num
	Data_Boards      struct {
		Snapshots
		Frequencies
		Neighborhood_Frequencies
	}
	Generation_i int
	//RandomState      rand.Rand
}

type Snapshots struct {
	Sequence []struct {
		Generation int
		Data       [][]int
	}
}
type Frequencies struct {
	Sequence []struct {
		Generation int
		Data       []int
	}
}
type Neighborhood_Frequencies struct {
	Sequence []struct {
		Generation int
		Data       [][]int
	}
}

type Parameters struct {
	Generations                  int
	R_Init_Odds                  float64
	S_Init_Odds                  float64
	Signal_Threshold             int
	Cooperation_Effect_Threshold int
	S_Radius                     int
	PG_Radius                    int
	Board_Size                   int
	D                            float64
	Mut_Odds_R                   float64
	Mut_Odds_S                   float64
	Basal_Cost                   float64
	Cooperation_Cost             float64
	Signal_Cost                  float64
	Receptor_Cost                float64
	Public_Goods_Effect          float64
}

type Settings struct {
	Data_Filename                string
	Snapshots_num                int
	Frequencies_num              int
	Neighborhood_Frequencies_num int
}

type Config struct {
	Parameters
	Settings
}

type Coordinate struct {
	r, c int
}

func (board *Board_strain) get_cell(c Coordinate) int {
	return (*board)[c.r][c.c]
}

func (board *Board_signal_num) get_cell(signal int, c Coordinate) int {
	return (*board)[signal][c.r][c.c]
}

func (board *Board_prod) get_cell(c Coordinate) bool {
	return (*board)[c.r][c.c]
}

func (board *Board_pg_num) get_cell(c Coordinate) int {
	return (*board)[c.r][c.c]
}

func (board *Board_strain) set_cell(c Coordinate, strain int) {
	(*board)[c.r][c.c] = strain
}

func (board *Board_signal_num) set_cell(signal_coordinate int, c Coordinate, signal_num int) {
	(*board)[signal_coordinate][c.r][c.c] = signal_num
}

func (board *Board_prod) set_cell(c Coordinate, prod bool) {
	(*board)[c.r][c.c] = prod
}

func (board *Board_pg_num) set_cell(c Coordinate, pg_num int) {
	(*board)[c.r][c.c] = pg_num
}

func rand_coord(board_size int) Coordinate {
	return Coordinate{
		r: rand.Intn(board_size),
		c: rand.Intn(board_size)}
}

func (board_strain Board_strain) String() (s string) {
	miscow.Trace("(Board_strain) String()")
	defer miscow.Untrace("(Board_strain) String()")
	s = " "
	for i := range board_strain {
		s += fmt.Sprintf("%2v", i)
	}
	s += "\n"
	for i0, v0 := range board_strain {
		s += fmt.Sprintf("%v: ", i0)
		for _, v1 := range v0 {
			s += fmt.Sprintf("%v ", v1)
		}
		s += fmt.Sprintf("\n")
	}
	s += "\n"
	return
}

func (board_signal_num Board_signal_num) String() (s string) {
	miscow.Trace("(Board_signal_num) String()")
	defer miscow.Untrace("(Board_signal_num) String()")
	for i0, v0 := range board_signal_num {
		s += fmt.Sprintf("strain: %2v\n", i0)
		for i1, v1 := range v0 {
			s += fmt.Sprintf("\n%v: ", i1)
			for i2, v2 := range v1 {
				s += fmt.Sprintf("(%v,%v) ", i2, v2)
			}
		}
		s += fmt.Sprintf("\n")
	}
	s += fmt.Sprintf("\n")
	return
}

func (board_prod Board_prod) String() (s string) {
	miscow.Trace("(Board_prod) String()")
	defer miscow.Untrace("(Board_prod) String()")
	for _, row_vals := range board_prod {
		for _, prod := range row_vals {
			switch prod {
			case true:
				s += fmt.Sprintf("%v", 1)
			case false:
				s += fmt.Sprintf("%v", 0)
			}
		}
		s += fmt.Sprintf("\n")
	}
	return
}

func (p Parameters) String() (s string) {
	s += fmt.Sprintln("Parameters:")
	s += fmt.Sprintf("Generations                %2v\n", p.Generations)
	s += fmt.Sprintf("RInitOdds                  %2v\n", p.R_Init_Odds)
	s += fmt.Sprintf("SInitOdds                  %2v\n", p.S_Init_Odds)
	s += fmt.Sprintf("STh                        %2v\n", p.Signal_Threshold)
	s += fmt.Sprintf("CooperationEffectThreshold %2v\n", p.Cooperation_Effect_Threshold)
	s += fmt.Sprintf("SRad                       %2v\n", p.S_Radius)
	s += fmt.Sprintf("CRad                       %2v\n", p.PG_Radius)
	s += fmt.Sprintf("BoardSize                  %2v\n", p.Board_Size)
	s += fmt.Sprintf("D                          %2v\n", p.D)
	s += fmt.Sprintf("MutOddsR                   %2v\n", p.Mut_Odds_R)
	s += fmt.Sprintf("MutOddsS                   %2v\n", p.Mut_Odds_S)
	return
}
func (model *Model) String() (s string) {
	miscow.Trace("(Model) String()")
	defer miscow.Untrace("(Model) String()")
	s += fmt.Sprintf("model.params:\n%v\n", model.Parameters)
	s += fmt.Sprintf("%v\n", *model.Board_strain)
	return
}

type Simulation interface {
	GetCellStrain(c Coordinate) int
	GetCellProd(c Coordinate) bool
	GetCellSignalNum(c Coordinate) int
	GetCellPGNum(c Coordinate) bool
	SetCellStrain(c Coordinate, val int)
	SetCellProd(c Coordinate, val bool)
	SetCellSignalNum(c Coordinate, val int)
	SetCellPGNum(c Coordinate, val bool)
}

func (model *Model) GetCellStrain(c Coordinate) int {
	return (*model.Board_strain)[c.r][c.c]
}

func (model *Model) GetCellProd(c Coordinate) bool {
	return (*model.Board_prod)[c.r][c.c]
}

func (model *Model) GetCellSignalNum(signal_type int, c Coordinate) int {
	return (*model.Board_signal_num)[signal_type][c.r][c.c]
}

func (model *Model) GetCellPGNum(c Coordinate) int {
	return (*model.Board_pg_num)[c.r][c.c]
}

func (model *Model) SetCellStrain(c Coordinate, val int) {
	(*model.Board_strain)[c.r][c.c] = val
}

func (model *Model) SetCellProd(c Coordinate, val bool) {
	(*model.Board_prod)[c.r][c.c] = val
}

func (model *Model) SetCellSignalNum(signal_type int, c Coordinate, val int) {
	(*model.Board_signal_num)[signal_type][c.r][c.c] = val
}

func (model *Model) SetCellPGNum(c Coordinate, val int) {
	(*model.Board_pg_num)[c.r][c.c] = val
}
