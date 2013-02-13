// The sequential version of Game of Strife
package sequential

import (
	"code.google.com/p/gcfg"
	"flag"
	"fmt"
	"log"
	"math/rand"
	"os"
	"runtime/pprof"
	"time"
)

func load_config() (parameters Parameters, settings Settings) {
	var cfg Config
	settings_filename := "strife.conf"
	err := gcfg.ReadFileInto(&cfg, settings_filename)
	if err != nil {
		log.Printf("settings file: %+v; error: %+v\n", settings_filename, err)
		parameters = Default_Parameters
		settings = Default_Settings
	} else {
		parameters = cfg.Parameters
		settings = cfg.Settings
	}
	return
}

func (m *Model) get_sample_rate() int {
	if m.Settings.Snapshots_num != 0 {
		return m.Parameters.Generations/m.Settings.Snapshots_num + 1
	}
	return 0
}

func Main(cpuprofile *string, imagesflag *bool) {
	fmt.Println(flag.Args())
	if *cpuprofile != "" {
		f, err := os.Create(*cpuprofile)
		if err != nil {
			log.Fatal(err)
		}
		pprof.StartCPUProfile(f)
		defer pprof.StopCPUProfile()
	}

	// Reading configuration file
	model := new(Model)
	params, settings := load_config()
	model.Start_Time = fmt.Sprintf("%v", time.Now().UnixNano())
	model.Parameters = params
	model.Settings = settings
	fmt.Printf("Parameters:\n%+v\n", model.Parameters)
	fmt.Printf("Settings:\n%+v\n", model.Settings)

	model.init_boards()
	model.init_data_samples()

	if err := model.save_json(); err != nil {
		panic(err)
	}

	model.run()

	if err := model.save_json(); err != nil {
		panic(err)
	}
	if *imagesflag {
		model.Save_snapshots_as_images()
	}

}

func (model *Model) get_diffusion_iter_per_generation() int {
	return int(model.Parameters.D * float64(model.Parameters.Board_Size*model.Parameters.Board_Size) / 4)
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
func (model *Model) run() {
	var snapshots_sample_rate int = model.get_sample_rate()
	var diffusion_num int = model.get_diffusion_iter_per_generation()
	fmt.Println("diffusion_num: ", diffusion_num)
	fmt.Println("model.Parameters.Generations = ", model.Parameters.Generations)
	fmt.Println("model.Generation_i = ", model.Generation_i)
	tstart := time.Now()

	for model.Generation_i = 0; model.Generation_i < model.Parameters.Generations; model.Generation_i++ {
		time_per_iter := time.Now()
		model.showboards()
		for competition_i := 0; competition_i < model.Parameters.Board_Size*model.Parameters.Board_Size; competition_i++ {
			model.competition()
		}
		model.showboards()

		for diffusion_i := 0; diffusion_i < diffusion_num; diffusion_i++ {
			if rand.Intn(2) == 0 {
				model.diffuse(Coordinate{r: rand.Intn(model.Parameters.Board_Size),
					c: rand.Intn(model.Parameters.Board_Size)},
					true)
			} else {
				model.diffuse(Coordinate{r: rand.Intn(model.Parameters.Board_Size),
					c: rand.Intn(model.Parameters.Board_Size)},
					false)
			}
		}
		model.showboards()

		model.showtiming(tstart, time.Since(time_per_iter))

		// take a snapshot only if we were asked to.
		if model.Settings.Snapshots_num != 0 {
			// take snapshot only every snapshots_sample_rate generations
			// or when we're at the last generation.
			if model.Generation_i%snapshots_sample_rate == 0 || model.Generation_i == model.Parameters.Generations-1 {
				model.take_strain_snapshot()
			}
		}
	}
}

func (model *Model) take_strain_snapshot() {
	for i := range model.Board_strain {
		model.Data_Boards.Snapshots[len(model.Data_Boards.Snapshots)-1].Data[i] = append([]int{}, model.Board_strain[i]...)
	}
	model.Data_Boards.Snapshots[len(model.Data_Boards.Snapshots)-1].Generation = model.Generation_i
	//	fmt.Println(model.Board_strain)

	/*	fmt.Println("Snapshots after assignment", model.Data_Boards.Snapshots[0],
		len(model.Data_Boards.Snapshots),
		cap(model.Data_Boards.Snapshots))
	*/

	if len(model.Data_Boards.Snapshots) < cap(model.Data_Boards.Snapshots) {
		model.Data_Boards.Snapshots = model.Data_Boards.Snapshots[0 : len(model.Data_Boards.Snapshots)+1]
	}
}

func (model *Model) showtiming(t_start time.Time, dt_iter time.Duration) {
	t_elapsed := time.Now().Sub(t_start)
	dt_tot_runtime := time.Duration(dt_iter.Nanoseconds()*int64(model.Parameters.Generations)) * time.Nanosecond
	dt_tot_runtime_300_10000 := time.Duration(dt_iter.Nanoseconds()*10000*300*300/int64(model.Parameters.Board_Size*model.Parameters.Board_Size)) * time.Nanosecond

	t_finish := t_start.Add(dt_tot_runtime)

	fmt.Println("Since start:", t_elapsed)
	fmt.Println("Expected total run time:", dt_tot_runtime)
	fmt.Println("Expected total run time (10k gen, 300^2 board):", dt_tot_runtime_300_10000)
	fmt.Println("Finish time:", t_finish)
}
