// The sequential version of Game of Strife
package sequential

import (
	"code.google.com/p/gcfg"
	"fmt"
	"log"
	"math/rand"
	"os"
	"runtime/pprof"
	"strife/flags"
	"time"
)

func load_config(cmdln_flags flags.Flags) (parameters Parameters, settings Settings) {
	var cfg Config
	var settings_filename string
	if cmdln_flags.Settings_filename_flag != "" {
		settings_filename = cmdln_flags.Settings_filename_flag
	} else {
		settings_filename = "default.conf"
	}
	err := gcfg.ReadFileInto(&cfg, settings_filename)
	if err != nil {
		log.Printf("settings file: %+v; error: %+v\n", settings_filename, err)
		parameters = DefaultParameters
		settings = DefaultSettings
	} else {
		parameters = cfg.Parameters
		settings = cfg.Settings
	}
	return
}

func (model *Model) SnapshotsSampleRate() int {
	if model.Settings.SnapshotsSampleNum != 0 {
		if model.Parameters.Generations%model.Settings.SnapshotsSampleNum == 0 {
			return model.Parameters.Generations / model.Settings.SnapshotsSampleNum
		} else {
			return model.Parameters.Generations/model.Settings.SnapshotsSampleNum + 1
		}
	}
	return model.Parameters.Generations
}

func (m *Model) FrequenciesSampleRate() int {
	if m.Settings.FrequenciesSampleNum != 0 {
		if m.Parameters.Generations%m.Settings.FrequenciesSampleNum == 0 {
			return m.Parameters.Generations / m.Settings.FrequenciesSampleNum
		} else {
			return m.Parameters.Generations/m.Settings.FrequenciesSampleNum + 1
		}
	}
	return m.Parameters.Generations
}

func Main(cmdln_flags flags.Flags) {
	if cmdln_flags.Cpuprofileflag != "" {
		f, err := os.Create(cmdln_flags.Cpuprofileflag)
		if err != nil {
			log.Fatal(err)
		}
		pprof.StartCPUProfile(f)
		defer pprof.StopCPUProfile()
	}

	// Reading configuration file
	params, settings := load_config(cmdln_flags)
	fmt.Printf("Parameters:\n%+v\n", params)
	fmt.Scanln()
	fmt.Printf("Settings:\n%+v\n", settings)
	fmt.Scanln()

	tstart := time.Now()

	// Which simulation model do we want to run?
	// Decided according to commandline flags.
	if cmdln_flags.Avigdorflag {
		model := new(AvigdorModel)
		model.CommandlineFlags = cmdln_flags
		model.setStartTime()
		model.Parameters = params
		model.Settings = settings
		model.initBoards()
		model.initDataSamples()
		if err := model.save_json(); err != nil {
			panic(err)
		}

		// TODO cmdln_flags should be inside model's settings
		var diffusion_num int = diffusionIterPerGeneration(model.Parameters.D, model.Parameters.BoardSize)
		fmt.Println("diffusion_num: ", diffusion_num)
		fmt.Println("model.Parameters.Generations = ", model.Parameters.Generations)
		fmt.Println("model.Generation_i = ", model.GenerationIdx)
		for model.GenerationIdx = 0; model.GenerationIdx < model.Parameters.Generations; model.GenerationIdx++ {
			time_per_iter := time.Now()
			for competition_i := 0; competition_i < model.Parameters.BoardSize*model.Parameters.BoardSize; competition_i++ {
				model.Competition()
			}

			for diffusion_i := 0; diffusion_i < diffusion_num; diffusion_i++ {
				if rand.Intn(2) == 0 {
					model.Diffuse(Coordinate{r: rand.Intn(model.Parameters.BoardSize),
						c: rand.Intn(model.Parameters.BoardSize)},
						true)
				} else {
					model.Diffuse(Coordinate{r: rand.Intn(model.Parameters.BoardSize),
						c: rand.Intn(model.Parameters.BoardSize)},
						false)
				}
			}

			model.sample()

			if model.GenerationIdx%10 == 0 {
				model.showtiming(tstart, time.Since(time_per_iter))
			}
		}
		if err := model.save_json(); err != nil {
			panic(err)
		}
		if cmdln_flags.Imagesflag {
			model.SaveSnapshotsAsImages()
		}
	}

	if cmdln_flags.Czaranflag {
		model := new(CzaranModel)
		model.CommandlineFlags = cmdln_flags
		model.setStartTime()
		model.Parameters = params
		model.Settings = settings
		model.initBoards()
		model.initDataSamples()
		if err := model.save_json(); err != nil {
			panic(err)
		}

		// TODO cmdln_flags should be inside model's settings
		var diffusion_num int = diffusionIterPerGeneration(model.Parameters.D, model.Parameters.BoardSize)
		fmt.Println("diffusion_num: ", diffusion_num)
		fmt.Println("model.Parameters.Generations = ", model.Parameters.Generations)
		fmt.Println("model.Generation_i = ", model.GenerationIdx)
		for model.GenerationIdx = 0; model.GenerationIdx < model.Parameters.Generations; model.GenerationIdx++ {
			time_per_iter := time.Now()
			for competition_i := 0; competition_i < model.Parameters.BoardSize*model.Parameters.BoardSize; competition_i++ {
				model.Competition()
			}

			for diffusion_i := 0; diffusion_i < diffusion_num; diffusion_i++ {
				if rand.Intn(2) == 0 {
					model.Diffuse(Coordinate{r: rand.Intn(model.Parameters.BoardSize),
						c: rand.Intn(model.Parameters.BoardSize)},
						true)
				} else {
					model.Diffuse(Coordinate{r: rand.Intn(model.Parameters.BoardSize),
						c: rand.Intn(model.Parameters.BoardSize)},
						false)
				}
			}

			model.sample()

			if model.GenerationIdx%10 == 0 {
				model.showtiming(tstart, time.Since(time_per_iter))
			}
		}
		if err := model.save_json(); err != nil {
			panic(err)
		}
		if cmdln_flags.Imagesflag {
			model.SaveSnapshotsAsImages()
		}
	}

}

func diffusionIterPerGeneration(d float64, boardSize int) int {
	return int(d * float64(boardSize*boardSize) / 4)
}

func (model *Model) showboards() {
	/*
		fmt.Println(model.Board_strain)
		fmt.Scanln()
		fmt.Println(model.Board_signal_num)
		fmt.Scanln()
		fmt.Println(model.Board_prod)
		fmt.Scanln()
		fmt.Println(model.Board_pg_num)
		fmt.Scanln()
	*/
}

// TODO maybe remove this when you're over?
func buga() {
	fmt.Println("buga")
	fmt.Println("buga")
	fmt.Println("buga")
	fmt.Println("buga")
	fmt.Println("buga")
	fmt.Println("buga")
}

func run(model Simulation, cmdln_flags flags.Flags) {
}

func (model *Model) sample() {
	var snapshotsSampleRate int = model.SnapshotsSampleRate()
	var frequenciesSampleRate int = model.FrequenciesSampleRate()
	// take a sample only if we were asked to.
	if model.Settings.SnapshotsSampleNum != 0 {
		// take sample only every snapshots_sample_rate generations
		// or when we're at the last generation.
		if model.GenerationIdx%snapshotsSampleRate == 0 || model.GenerationIdx == model.Parameters.Generations-1 {
			model.take_board_sample()
		}
	}

	// take a sample only if we were asked to.
	if model.Settings.FrequenciesSampleNum != 0 {
		// take sample only every frequencies_sample_rate generations
		// or when we're at the last generation.
		if model.GenerationIdx%frequenciesSampleRate == 0 || model.GenerationIdx == model.Parameters.Generations-1 {
			model.take_frequencies_sample()
			model.take_neighbors_frequencies_sample()
		}
	}
}

func (model *Model) take_board_sample() {
	for i := range model.BoardStrain {
		model.DataSamples.Snapshots[len(model.DataSamples.Snapshots)-1].Data[i] = append([]int{}, model.BoardStrain[i]...)
	}

	model.DataSamples.Snapshots[len(model.DataSamples.Snapshots)-1].Generation = model.GenerationIdx

	if len(model.DataSamples.Snapshots) < cap(model.DataSamples.Snapshots) {
		model.DataSamples.Snapshots = model.DataSamples.Snapshots[0 : len(model.DataSamples.Snapshots)+1]
	}
}

func (model *Model) take_frequencies_sample() {
	coord := Coordinate{}
	for coord.r = 0; coord.r < model.Parameters.BoardSize; coord.r++ {
		for coord.c = 0; coord.c < model.Parameters.BoardSize; coord.c++ {
			model.DataSamples.Frequencies[len(model.DataSamples.Frequencies)-1].Data[model.CellStrain(coord)]++
		}
	}

	model.DataSamples.Frequencies[len(model.DataSamples.Frequencies)-1].Generation = model.GenerationIdx

	if len(model.DataSamples.Frequencies) < cap(model.DataSamples.Frequencies) {
		model.DataSamples.Frequencies = model.DataSamples.Frequencies[0 : len(model.DataSamples.Frequencies)+1]
	}
}

func (model *Model) take_neighbors_frequencies_sample() {
	center_coord := Coordinate{}
	for center_coord.r = 0; center_coord.r < model.Parameters.BoardSize; center_coord.r++ {
		for center_coord.c = 0; center_coord.c < model.Parameters.BoardSize; center_coord.c++ {
			rad_coord := Coordinate{}
			for rad_coord.r = center_coord.r - 1; rad_coord.r < center_coord.r+1; rad_coord.r++ {
				for rad_coord.c = center_coord.r - 1; rad_coord.c < center_coord.r+1; rad_coord.c++ {
					model.DataSamples.NeighborsFrequencies[len(model.DataSamples.NeighborsFrequencies)-1].Data[model.CellStrain(center_coord)][model.CellStrain(rad_coord.ToroidCoordinates(model.Parameters.BoardSize))]++
				}
			}
		}
	}

	model.DataSamples.NeighborsFrequencies[len(model.DataSamples.NeighborsFrequencies)-1].Generation = model.GenerationIdx

	if len(model.DataSamples.NeighborsFrequencies) < cap(model.DataSamples.NeighborsFrequencies) {
		model.DataSamples.NeighborsFrequencies = model.DataSamples.NeighborsFrequencies[0 : len(model.DataSamples.NeighborsFrequencies)+1]
	}
}

func (model *Model) showtiming(t_start time.Time, dt_iter time.Duration) {
	t_elapsed := time.Now().Sub(t_start)
	dt_tot_runtime := time.Duration(dt_iter.Nanoseconds()*int64(model.Parameters.Generations)) * time.Nanosecond
	dt_tot_runtime_300_10000 := time.Duration(dt_iter.Nanoseconds()*10000*300*300/int64(model.Parameters.BoardSize*model.Parameters.BoardSize)) * time.Nanosecond
	dt_eta := time.Duration(dt_iter.Nanoseconds() * int64(model.Parameters.Generations - model.GenerationIdx)) * time.Nanosecond
	dt_eta_300_10000 := time.Duration(dt_iter.Nanoseconds() * int64(model.Parameters.Generations - model.GenerationIdx) * 10000 * 300 * 300 / int64(model.Parameters.BoardSize*model.Parameters.BoardSize)) * time.Nanosecond

	t_finish := t_start.Add(dt_tot_runtime)

	fmt.Println("Since start:", t_elapsed)
    fmt.Printf("Generation (%v tot): %v\n", model.Parameters.Generations, model.GenerationIdx)
	fmt.Println("Expected total run time:", dt_tot_runtime)
	fmt.Println("Time to finish:", dt_eta)
	fmt.Println("Expected total run time (10k gen, 300^2 board):", dt_tot_runtime_300_10000)
	fmt.Println("Time to finish (10k gen, 300^2 board):", dt_eta_300_10000)
	fmt.Println("Finish time:", t_finish)
}
