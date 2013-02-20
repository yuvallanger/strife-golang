package sequential

import (
	"fmt"
	"math/rand"
	"miscow"
	"strife/flags"
	"time"
)

var r4strain = [8]int{0, 1, 0, 1, 0, 1, 0, 1} // index is the strain. value is the receptor allele.
var s4strain = [8]int{0, 0, 1, 1, 0, 0, 1, 1} // index is the strain. value is the signal allele.
var g4strain = [8]int{0, 0, 0, 0, 1, 1, 1, 1} // index is the strain. value is the receptor allele.

type BoardStrain [][]int      // [rows][columns]. possible values: {0,1,2,3}. s0r0 - 0; s0r1 - 1; s1r0 - 2; s1r1 - 3
type BoardSignalNum [][][]int // [signal types][rows][columns]
type BoardProd [][]bool       // [rows][columns]
type BoardPGNum [][]int       // [rows][columns]

// TODO This interface isn't very helpfull
type Simulation interface {
	showtiming(time.Time, time.Duration)
	getGenerationIdx() int
	setGenerationIdx(int)
	incGenerationIdx()
	setStartTime()
	setSettings(Settings)
	setParameters(Parameters)
	getSettings() Settings
	getParameters() Parameters
	initBoards()
	initBoardStrain()
	initDataSamples()
	sample()
	save_json() error
	SaveSnapshotsAsImages()
	Fitness(Coordinate) float64
	mutate(Coordinate) int
	updateArrays(Coordinate, int, int)
	Endgame(winnerCoord, loserCoord Coordinate)
	Competition()
	Diffuse(coord00 Coordinate, direction bool)
}

func (model *Model) getSettings() Settings {
	return model.Settings
}

func (model *Model) setStartTime() {
	model.StartTime = int64(time.Now().UnixNano())
}

type AvigdorModel struct {
	Model
}

type CzaranModel struct {
	Model
}

type Model struct {
	CommandlineFlags          flags.Flags
	Parameters     Parameters
	Settings       Settings
	BoardStrain    BoardStrain
	BoardSignalNum BoardSignalNum
	BoardProd      BoardProd
	BoardPGNum     BoardPGNum
	DataSamples    struct {
		Snapshots            []Snapshot
		Frequencies          []Frequency
		NeighborsFrequencies []NeighborsFrequency
	}
	GenerationIdx int
	//RandomState      rand.Rand
	StartTime int64 // start time in unix nanoseconds
}

type Snapshot struct {
	Generation int
	Data       BoardStrain
}

type Snapshots []Snapshot

type Frequency struct {
	Generation int
	Data       []int
}

type Frequencies struct {
	Frequency []Frequency
}

type NeighborsFrequency struct {
	Generation int
	Data       [][]int
}

type NeighborsFrequencies struct {
	NeighborsFrequency []NeighborsFrequency
}

type Parameters struct {
	Generations                int     // TODO
	RInitOdds                  float64 // TODO
	SInitOdds                  float64 // TODO
	GInitOdds                  float64 // TODO
	SignalThreshold            int     // TODO
	CooperationEffectThreshold int     // TODO
	SRadius                    int     // the moor radius of signal production
	PGRadius                   int     // the moor radius of public goods production
	BoardSize                  int     // a board will be sized BoardSize^2
	D                          float64 // diffusion coefficient
	MutOddsR                   float64 // TODO
	MutOddsS                   float64 // TODO
	MutOddsG                   float64 // TODO
	BasalCost                  float64 // TODO
	CooperationCost            float64 // TODO
	SignalCost                 float64 // TODO
	ReceptorCost               float64 // TODO
	PublicGoodsEffect          float64 // TODO
}

type Settings struct {
	DataFilename                     string // TODO
	SnapshotsSampleNum               int    // TODO
	FrequenciesSampleNum             int    // TODO
	NeighborhoodFrequenciesSampleNum int    // TODO
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
func (c *Coordinate) ToroidCoordinates(boardSize int) Coordinate {
	return Coordinate{r: miscow.MyMod(c.r, boardSize),
		c: miscow.MyMod(c.c, boardSize)}
}

func RandCoord(boardSize int) Coordinate {
	return Coordinate{
		r: rand.Intn(boardSize),
		c: rand.Intn(boardSize)}
}

func (board BoardStrain) String() (s string) {
	s = " "
	for i := range board {
		s += fmt.Sprintf("%2v", i)
	}
	s += "\n"
	for i0, v0 := range board {
		s += fmt.Sprintf("%v: ", i0)
		for _, v1 := range v0 {
			s += fmt.Sprintf("%v ", v1)
		}
		s += fmt.Sprintf("\n")
	}
	s += "\n"
	return
}

func (board BoardSignalNum) String() (s string) {
	for i0, v0 := range board {
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

func (board BoardProd) String() (s string) {
	for _, rowVals := range board {
		for _, prod := range rowVals {
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

// fmt uses this function to generatie a string out of Parameters typed variables.
// TODO: must convert all names into underscore notation.
func (p Parameters) String() (s string) {
	s += fmt.Sprintln("Parameters:")
	s += fmt.Sprintln("Generations")
	s += fmt.Sprintln("    ", p.Generations)
	s += fmt.Sprintln("RInitOdds")
	s += fmt.Sprintln("    ", p.RInitOdds)
	s += fmt.Sprintln("SInitOdds")
	s += fmt.Sprintln("    ", p.SInitOdds)
	s += fmt.Sprintln("GInitOdds")
	s += fmt.Sprintln("    ", p.GInitOdds)
	s += fmt.Sprintln("SignalThreshold")
	s += fmt.Sprintln("    ", p.SignalThreshold)
	s += fmt.Sprintln("CooperationEffectThreshold")
	s += fmt.Sprintln("    ", p.CooperationEffectThreshold)
	s += fmt.Sprintln("SRad")
	s += fmt.Sprintln("    ", p.SRadius)
	s += fmt.Sprintln("CRad")
	s += fmt.Sprintln("    ", p.PGRadius)
	s += fmt.Sprintln("BoardSize")
	s += fmt.Sprintln("    ", p.BoardSize)
	s += fmt.Sprintln("D")
	s += fmt.Sprintln("    ", p.D)
	s += fmt.Sprintln("MutOddsR")
	s += fmt.Sprintln("    ", p.MutOddsR)
	s += fmt.Sprintln("MutOddsS")
	s += fmt.Sprintln("    ", p.MutOddsS)
	return
}

func (settings Settings) String() (s string) {
	s += fmt.Sprintln("Settings:")
	s += fmt.Sprintln("DataFilename:")
	s += fmt.Sprintln("   ", settings.DataFilename)
	s += fmt.Sprintln("SnapshotsSampleNum:")
	s += fmt.Sprintln("   ", settings.SnapshotsSampleNum)
	s += fmt.Sprintln("FrequenciesSampleNum:")
	s += fmt.Sprintln("   ", settings.FrequenciesSampleNum)
	s += fmt.Sprintln("NeighborhoodFrequenciesSampleNum:")
	s += fmt.Sprintln("   ", settings.NeighborhoodFrequenciesSampleNum)
	return
}

func (model *Model) String() (s string) {
	s += fmt.Sprintf("model.params:\n%v\n", model.Parameters)
	s += fmt.Sprintf("%v\n", model.BoardStrain)
	return
}

func (model *Model) CellStrain(c Coordinate) int {
	return (model.BoardStrain)[c.r][c.c]
}

func (model *Model) CellProd(c Coordinate) bool {
	return (model.BoardProd)[c.r][c.c]
}

func (model *Model) CellSignalNum(c Coordinate, signalType int) int {
	return (model.BoardSignalNum)[signalType][c.r][c.c]
}

func (model *Model) CellPGNum(c Coordinate) int {
	return model.BoardPGNum[c.r][c.c]
}

func (model *Model) SetCellStrain(c Coordinate, val int) {
	(model.BoardStrain)[c.r][c.c] = val
}

func (model *Model) SetCellProd(c Coordinate, val bool) {
	(model.BoardProd)[c.r][c.c] = val
}

func (model *Model) setCellSignalNum(c Coordinate, signalType, val int) {
	(model.BoardSignalNum)[signalType][c.r][c.c] = val
}

func (model *Model) setCellPGNum(c Coordinate, val int) {
	(model.BoardPGNum)[c.r][c.c] = val
}

func (model *Model) AddToCellPGNum(c Coordinate, val int) {
	(model.BoardPGNum)[c.r][c.c] += val
}

func (model *Model) AddToCellSignalNum(c Coordinate, signalType, val int) {
	model.BoardSignalNum[signalType][c.r][c.c] += val
}
