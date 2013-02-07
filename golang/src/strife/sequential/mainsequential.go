package sequential

import (
	"code.google.com/p/gcfg"
	"encoding/json"
	"fmt"
	"io/ioutil"
	"miscow"
)

func load_config() (parameters Parameters, settings Settings) {
    miscow.Trace("load_config")
	defer miscow.Untrace("load_config")
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

func get_sample_rate(m *Model, n int) (sample_rate int) {
    miscow.Trace("get_sample_rate")
    defer miscow.Untrace("get_sample_rate")
	if n == 0 {
		sample_rate = 0
	} else {
		sample_rate = m.Parameters.Board_Size/m.Settings.Snapshots_num + 1
	}
	return sample_rate
}

func init_databoards(model *Model) {
    miscow.Trace("init_databoards")
    defer miscow.Untrace("init_databoards")
	sample_rate := get_sample_rate(model, model.Settings.Snapshots_num)
	if sample_rate != 0 {
		model.Data_Boards.Snapshots.Sequence = make([]struct {
			Generation int
			Data       [][]int
		}, model.Settings.Snapshots_num)
		for sample_i := 0; sample_i < model.Settings.Snapshots_num+1; sample_i++ {
			for row_i := 0; row_i < model.Parameters.Board_Size; row_i++ {
				fmt.Print(model.Data_Boards.Snapshots.Sequence[sample_i])
				//model.Data_Boards.Snapshots.Sequence[sample_i].Data = make([][]int, model.Parameters.Board_Size)
				for col_i := 0; col_i < model.Parameters.Board_Size; col_i++ {
					if model.Settings.Snapshots_num != 0 {
						model.Data_Boards.Snapshots.Sequence[sample_i].Data[row_i] = make([]int, model.Parameters.Board_Size)
					}
				}
			}
		}
	}
}

func Main() {
    miscow.Trace("Main")
    defer miscow.Untrace("Main")
	// Reading configuration file
	model := new(Model)
	params, settings := load_config()
	model.Parameters = params
	model.Settings = settings
	fmt.Printf("Parameters:\n%+v\n", model.Parameters)
	fmt.Printf("Settings:\n%+v\n", model.Settings)

	init_boards(model)
	init_databoards(model)

	jsonModel, err := json.Marshal(model)
	switch err != nil {
	case true:
		fmt.Println("error:", err)
		return
	case false:
		f, err := ioutil.TempFile("/home/yankel/", "strife-temp")
		defer f.Close()
		switch err != nil {
		case true:
			fmt.Println("error:", err)
		case false:
			f.Write(jsonModel)
		}

		//fmt.Println("Board strain:\n", model.Board_strain)

		run(model)

		//fmt.Println("Board strain:\n", model.Board_strain)
		//fmt.Println("model.params:\n", model.Parameters)
		//fmt.Println("model.settings:\n", model.Settings)
		fmt.Print()
	}
}
