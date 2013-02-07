package concurrent

import (
	"fmt"
	"math/rand"
	"miscow"
	"runtime"
	"time"
)

type Coordinate struct{ i, j int }

func Main() {
	ncpu := runtime.NumCPU()
	fmt.Printf("previous GOMAXPROCS: %+v; now: %+v\n", runtime.GOMAXPROCS(ncpu), ncpu)
	fmt.Scanln()
	var m Concur_Model = Concur_Model{
		N: 5}
	m = init_comms(m)
	fire_cells(m)
}

type CellComm struct {
	CellState CellState
	Action    struct {
		inc_signal bool
		dec_signal bool
		inc_pg     bool
		dec_pg     bool
	}
}

type Concur_Model struct {
	CommBoard [][]chan CellComm
	N         int
}

func init_comms(m Concur_Model) Concur_Model {
	m.CommBoard = make([][]chan CellComm, m.N)
	for i := range m.CommBoard {
		(m.CommBoard)[i] = make([]chan CellComm, m.N)
		for j := range (m.CommBoard)[i] {
			(m.CommBoard)[i][j] = make(chan CellComm)
		}
	}
	return m
}

type CellState struct {
	signal      bool
	receptor    bool
	cooperation bool
	signal_num  int
	pg_num      int
}

// this will fire all the cells with their initial state
func fire_cells(m Concur_Model) {
	var res chan int
	// just a placeholder for future settings
	//initstate := CellState{signal: true, receptor: true, cooperation: true}
	var neighbors_chans [8]chan CellComm
	var inbox_chan chan CellComm
	// for the toroid coordinates, we'll create:
	// ip = i plus  = i+1
	// im = i minus = i-1
	// jp = j plus  = j+1
	// jm = j minus = j-1
	var im, ip, jm, jp int
	for i := range m.CommBoard {
		for j := range (m.CommBoard)[i] {
			im = miscow.MyMod(i-1, m.N)
			ip = miscow.MyMod(i+1, m.N)
			jm = miscow.MyMod(i-1, m.N)
			jp = miscow.MyMod(i+1, m.N)
			neighbors_chans[0] = (m.CommBoard)[im][jm]
			neighbors_chans[1] = (m.CommBoard)[im][j]
			neighbors_chans[2] = (m.CommBoard)[im][jp]
			neighbors_chans[3] = (m.CommBoard)[i][jp]
			neighbors_chans[4] = (m.CommBoard)[ip][jp]
			neighbors_chans[5] = (m.CommBoard)[ip][j]
			neighbors_chans[6] = (m.CommBoard)[ip][jm]
			neighbors_chans[7] = (m.CommBoard)[i][jm]
			inbox_chan = (m.CommBoard)[i][j]
			go cell_go(rand.New(rand.NewSource(time.Now().UnixNano())),
				CellState{signal: rand.Float64() < 0.5},
				neighbors_chans,
				inbox_chan,
				Coordinate{i, j},
				res)
		}
	}
	sum := 0
	for i := 0; i < m.N*m.N; i++ {
        a_res := <-res
        fmt.Println(a_res)
        sum += a_res
	}
	fmt.Println("Done. Sum:", sum)
	fmt.Scanln()
}

func cell_go(r *rand.Rand, initstate CellState, neighbors_chans [8]chan CellComm, inbox_chan chan CellComm, id Coordinate, res chan int) {
	mystate := CellState{
		signal:      initstate.signal,
		receptor:    initstate.receptor,
		cooperation: initstate.cooperation}
	a_comm := CellComm{}
	for i := 0; i < 8; i++ {
		go func(neighbor_chan chan CellComm) {
			neighbor_chan <- CellComm{
				CellState: CellState{signal: mystate.signal}}
		}(neighbors_chans[i])
	}
    fmt.Printf("id: %v; after signal\n", id)
	for i := 0; i < 8; i++ {
		a_comm = <-inbox_chan
        fmt.Printf("id: %v; a_comm: %v\n", id, a_comm)
		switch a_comm {
		case CellComm{CellState: CellState{signal: true}}:
            fmt.Printf("id: %v; signal_num++", id)
		 	mystate.signal_num++
        case CellComm{CellState: CellState{signal:false}}:
            fmt.Printf("id: %v; signal_num=", id)
            continue
		}
	}
	res <- mystate.signal_num
}
