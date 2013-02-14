package sequential

import (
	"fmt"
	"math/rand"
	"miscow"
)

var r4strain = [8]int{0, 1, 0, 1, 0, 1, 0, 1} // index is the strain. value is the receptor allele.
var s4strain = [8]int{0, 0, 1, 1, 0, 0, 1, 1} // index is the strain. value is the signal allele.
var g4strain = [8]int{0, 0, 0, 0, 1, 1, 1, 1} // index is the strain. value is the receptor allele.

type Board_strain [][]int       // [rows][columns]. possible values: {0,1,2,3}. s0r0 - 0; s0r1 - 1; s1r0 - 2; s1r1 - 3
type Board_signal_num [][][]int // [signal types][rows][columns]
type Board_prod [][]bool        // [rows][columns]
type Board_pg_num [][]int       // [rows][columns]

type Model struct {
	Parameters       Parameters
	Settings         Settings
	Board_strain     Board_strain
	Board_signal_num Board_signal_num
	Board_prod       Board_prod
	Board_pg_num     Board_pg_num
	Data_samples     struct {
		Snapshots             []Snapshot
		Frequencies           []Frequency
		Neighbors_frequencies []Neighbors_frequency
	}
	Generation_i int
	//RandomState      rand.Rand
	Start_Time string // start time in unix nanoseconds
}

type Snapshot struct {
	Generation int
	Data       Board_strain
}

type Snapshots []Snapshot

type Frequency struct {
	Generation int
	Data       []int
}

type Frequencies struct {
	Frequency []Frequency
}

type Neighbors_frequency struct {
	Generation int
	Data       [][]int
}

type Neighbors_frequencies struct {
	Neighbors_frequency []Neighbors_frequency
}

type Parameters struct {
	Generations                  int     // TODO
	R_Init_Odds                  float64 // TODO
	S_Init_Odds                  float64 // TODO
	Signal_Threshold             int     // TODO
	Cooperation_Effect_Threshold int     // TODO
	S_Radius                     int     // the moor radius of signal production
	PG_Radius                    int     // the moor radius of public goods production
	Board_Size                   int     // a board will be sized Board_Size^2
	D                            float64 // diffusion coefficient
	Mut_Odds_R                   float64 // TODO
	Mut_Odds_S                   float64 // TODO
	Basal_Cost                   float64 // TODO
	Cooperation_Cost             float64 // TODO
	Signal_Cost                  float64 // TODO
	Receptor_Cost                float64 // TODO
	Public_Goods_Effect          float64 // TODO
}

type Settings struct {
	Data_filename                       string // TODO
	Snapshots_sample_num                int    // TODO
	Frequencies_sample_num              int    // TODO
	Neighborhood_frequencies_sample_num int    // TODO
}

type Config struct {
	Parameters Parameters // TODO
	Settings   Settings   // TODO
}

type Coordinate struct {
	r, c int // TODO
}

/*
Returns the same coordinates in modulus the size of the board.
Board is a square, not a rectangle, so we need only one size argument.
*/
func (c *Coordinate) get_toroid_coordinates(board_size int) Coordinate {
	return Coordinate{r: miscow.MyMod(c.r, board_size),
		c: miscow.MyMod(c.c, board_size)}
}

func rand_coord(board_size int) Coordinate {
	return Coordinate{
		r: rand.Intn(board_size),
		c: rand.Intn(board_size)}
}

func (board_strain Board_strain) String() (s string) {
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
	s += fmt.Sprintln("Generations")
	s += fmt.Sprintln("    ", p.Generations)
	s += fmt.Sprintln("RInitOdds")
	s += fmt.Sprintln("    ", p.R_Init_Odds)
	s += fmt.Sprintln("SInitOdds")
	s += fmt.Sprintln("    ", p.S_Init_Odds)
	s += fmt.Sprintln("STh")
	s += fmt.Sprintln("    ", p.Signal_Threshold)
	s += fmt.Sprintln("CooperationEffectThreshold")
	s += fmt.Sprintln("    ", p.Cooperation_Effect_Threshold)
	s += fmt.Sprintln("SRad")
	s += fmt.Sprintln("    ", p.S_Radius)
	s += fmt.Sprintln("CRad")
	s += fmt.Sprintln("    ", p.PG_Radius)
	s += fmt.Sprintln("BoardSize")
	s += fmt.Sprintln("    ", p.Board_Size)
	s += fmt.Sprintln("D")
	s += fmt.Sprintln("    ", p.D)
	s += fmt.Sprintln("MutOddsR")
	s += fmt.Sprintln("    ", p.Mut_Odds_R)
	s += fmt.Sprintln("MutOddsS")
	s += fmt.Sprintln("    ", p.Mut_Odds_S)
	return
}

func (settings Settings) String() (s string) {
	s += fmt.Sprintln("Settings:")
	s += fmt.Sprintln("Data_Filename:")
	s += fmt.Sprintln("   ", settings.Data_filename)
	s += fmt.Sprintln("Snapshots_sample_num:")
	s += fmt.Sprintln("   ", settings.Snapshots_sample_num)
	s += fmt.Sprintln("Frequencies_sample_num:")
	s += fmt.Sprintln("   ", settings.Frequencies_sample_num)
	s += fmt.Sprintln("Neighborhood_Frequencies_sample_num:")
	s += fmt.Sprintln("   ", settings.Neighborhood_frequencies_sample_num)
	return
}

func (model *Model) String() (s string) {
	s += fmt.Sprintf("model.params:\n%v\n", model.Parameters)
	s += fmt.Sprintf("%v\n", model.Board_strain)
	return
}

func (model *Model) get_cell_strain(c Coordinate) int {
	return (model.Board_strain)[c.r][c.c]
}

func (model *Model) get_cell_prod(c Coordinate) bool {
	return (model.Board_prod)[c.r][c.c]
}

func (model *Model) get_cell_signal_num(c Coordinate, signal_type int) int {
	return (model.Board_signal_num)[signal_type][c.r][c.c]
}

func (model *Model) get_cell_pg_num(c Coordinate) int {
	return model.Board_pg_num[c.r][c.c]
}

func (model *Model) set_cell_strain(c Coordinate, val int) {
	(model.Board_strain)[c.r][c.c] = val
}

func (model *Model) set_cell_prod(c Coordinate, val bool) {
	(model.Board_prod)[c.r][c.c] = val
}

func (model *Model) set_cell_signal_num(c Coordinate, signal_type, val int) {
	(model.Board_signal_num)[signal_type][c.r][c.c] = val
}

func (model *Model) set_cell_pg_num(c Coordinate, val int) {
	(model.Board_pg_num)[c.r][c.c] = val
}

func (model *Model) Add_To_Cell_PG_Num(c Coordinate, val int) {
	(model.Board_pg_num)[c.r][c.c] += val
}

func (model *Model) Add_To_Cell_Signal_Num(c Coordinate, signal_type, val int) {
	model.Board_signal_num[signal_type][c.r][c.c] += val
}
