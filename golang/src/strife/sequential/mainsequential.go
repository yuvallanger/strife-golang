// The sequential version of Game of Strife
package sequential

import (
	"code.google.com/p/gcfg"
	"fmt"
	"math/rand"
	"miscow"
	"time"
)

func load_config() (parameters Parameters, settings Settings) {
	var cfg Config
	err := gcfg.ReadFileInto(&cfg, "strife.conf")
	if err != nil {
		fmt.Println(err)
		parameters = Default_Parameters
		settings = Default_Settings
	} else {
		parameters = cfg.Parameters
		settings = cfg.Settings
	}
	return
}

func get_sample_rate(m *Model) int {
	if m.Settings.Snapshots_num != 0 {
		return m.Parameters.Generations/m.Settings.Snapshots_num + 1
	}
	return 0
}

func init_databoard_snapshots(model *Model) {
	if model.Settings.Snapshots_num != 0 {
		if model.Parameters.Generations%model.Settings.Snapshots_num == 0 {
			model.Data_Boards.Snapshots = make([]Snapshot, model.Settings.Snapshots_num)
		} else {
			// if we have number of generations that does not divides with number of snapshots, we want another snapshot at the end.
			model.Data_Boards.Snapshots = make([]Snapshot, model.Settings.Snapshots_num+1)
		}
		for snapshot_i := range model.Data_Boards.Snapshots {
			data := miscow.Make2dIntArray(model.Parameters.Board_Size, model.Parameters.Board_Size)
			model.Data_Boards.Snapshots[snapshot_i].Data = data
		}
	}
}

func init_databoards(model *Model) {
	init_databoard_snapshots(model)
}

func Main() {
	// Reading configuration file
	model := new(Model)
	params, settings := load_config()
	model.Start_Time = fmt.Sprintf("%v", time.Now().UnixNano())
	model.Parameters = params
	model.Settings = settings
	fmt.Printf("Parameters:\n%+v\n", model.Parameters)
	fmt.Printf("Settings:\n%+v\n", model.Settings)

	init_boards(model)
	init_databoards(model)

	err := save_json(model)
	if err != nil {
		panic(fmt.Sprintf("save_json() error: %v", err))
	}

	//fmt.Println("Board strain:\n", model.Board_strain)

	run(model)

	save_json(model)
	//fmt.Println("Board strain:\n", model.Board_strain)
	//fmt.Println("model.params:\n", model.Parameters)
	//fmt.Println("model.settings:\n", model.Settings)
	fmt.Print()
}

func run(model *Model) {
	var sample_rate int = get_sample_rate(model)
	diffusion_num := int(model.Parameters.D * float64(model.Parameters.Board_Size*model.Parameters.Board_Size) / 4)
	fmt.Println("diffusion_num: ", diffusion_num)
	fmt.Println("model.Parameters.Generations = ", model.Parameters.Generations)
	fmt.Println("model.Generation_i = ", model.Generation_i)
	tstart := time.Now()
	for model.Generation_i = 0; model.Generation_i < model.Parameters.Generations; model.Generation_i++ {
		time_per_iter := time.Now()

		board_size := model.Parameters.Board_Size
		//fmt.Printf("Generation: %v\n", model.Generation_i)
		for competition_i := 0; competition_i < board_size*board_size; competition_i++ {
			competition(model)
		}
		for diffusion_i := 0; diffusion_i < diffusion_num; diffusion_i++ {
			if rand.Intn(2) == 0 {
				diffuse(model,
					Coordinate{r: rand.Intn(model.Parameters.Board_Size),
						c: rand.Intn(model.Parameters.Board_Size)},
					true)
			} else {
				diffuse(model,
					Coordinate{r: rand.Intn(model.Parameters.Board_Size),
						c: rand.Intn(model.Parameters.Board_Size)},
					false)
			}
		}

		if model.Generation_i%100 == 0 {
			showtiming(tstart, time.Since(time_per_iter))
		}

		//func (t Time) Add(d Duration) Time
		//func Since(t Time) Duration
		//func (t Time) Sub(u Time) Duration
		if model.Settings.Snapshots_num != 0 {
			if model.Generation_i%sample_rate == 0 || model.Generation_i == model.Parameters.Generations-1 {
				fmt.Printf("data to json: generation: %v\n", model.Generation_i)
				fmt.Printf("data to json: model.Generation_i/sample_rate: %v\n",
					model.Generation_i/sample_rate)
				model.Data_Boards.Snapshots[model.Generation_i/sample_rate] = Snapshot{
					Data:       model.Board_strain,
					Generation: model.Generation_i}
			}
		}
	}
	fmt.Println()
	return
}
