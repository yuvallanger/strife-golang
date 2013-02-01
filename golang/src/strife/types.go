package main

import (
	"fmt"
	"miscow"
)

type Board_strain [][]int       // [rows][columns]. possible values: {0,1,2,3}. s0r0 - 0; s0r1 - 1; s1r0 - 2; s1r1 - 3
type Board_signal_num [][][]int // [signal types][rows][columns]
type Board_prod [][]int         // [rows][columns]
type Board_pg_num [][]int       // [rows][columns]

type Model struct {
	Parameters       Parameters_T
	Settings         Settings_T
	Board_strain     *Board_strain
	Board_signal_num *Board_signal_num
	Board_prod       *Board_prod
	Board_pg_num     *Board_pg_num
	DataBoards       struct {
		Snapshots                *Snapshots_T
		Frequencies              *Frequencies_T
		Neighborhood_Frequencies *Neighborhood_Frequencies_T
	}
	Generation_i int
	//RandomState      rand.Rand
}

type Snapshots_T [][][]int
type Frequencies_T [][]int
type Neighborhood_Frequencies_T [][][]int

type Parameters_T struct {
	Generations int
	RInitOdds   float64
	SInitOdds   float64
	STh         int
	CTh         int
	SRad        int
	CRad        int
	BoardSize   int
	D           int
	MutOddsR    float64
	MutOddsS    float64
}

type Settings_T struct {
	DataFilename string
}

type Config struct {
	Parameters Parameters_T
	Settings   Settings_T
}

func (board_strain Board_strain) String() (s string) {
	miscow.Trace("(Board_strain) String()")
	defer miscow.Untrace("(Board_strain) String()")
	s = ""
	for _, gene := range [2]int{1, 2} {
		s += fmt.Sprintf("gene %v\n", gene)
		for i0, v0 := range board_strain {
			s += fmt.Sprintf("%v: ", i0)
			for i1, v1 := range v0 {
				s += fmt.Sprintf("(%v,%v) ", i1, (v1&gene)/gene)
			}
			s += fmt.Sprintf("\n")
		}
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
	for _, gene := range [2]int{1, 2} {
		s += fmt.Sprintf("gene %v\n", gene)
		for _, v0 := range board_prod {
			for _, v1 := range v0 {
				s += fmt.Sprintf("(g%v ", v1)
				s += fmt.Sprintf("%v) ", (v1&gene)/gene)
			}
			s += fmt.Sprintf("\n")
		}
	}
	return
}

func (p Parameters_T) String() (s string) {
	s += fmt.Sprintf("Generations %v\n", p.Generations)
	s += fmt.Sprintf("RInitOdds   %v\n", p.RInitOdds)
	s += fmt.Sprintf("SInitOdds   %v\n", p.SInitOdds)
	s += fmt.Sprintf("STh         %v\n", p.STh)
	s += fmt.Sprintf("CTh         %v\n", p.CTh)
	s += fmt.Sprintf("SRad        %v\n", p.SRad)
	s += fmt.Sprintf("CRad        %v\n", p.CRad)
	s += fmt.Sprintf("BoardSize   %v\n", p.BoardSize)
	s += fmt.Sprintf("D           %v\n", p.D)
	s += fmt.Sprintf("MutOddsR    %v\n", p.MutOddsR)
	s += fmt.Sprintf("MutOddsS    %v\n", p.MutOddsS)
	return
}
func (model Model) String() (s string) {
	miscow.Trace("(Model) String()")
	defer miscow.Untrace("(Model) String()")
	s += fmt.Sprintf("model.params:\n%v\n", model.Parameters)
	s += fmt.Sprintf("%v\n", *model.Board_strain)
	return
}
