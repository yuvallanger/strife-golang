package sequential

import (
	"encoding/json"
	"fmt"
	"io/ioutil"
	"log"
)

/*
func (model *Model) save_h5() {
	var data struct {
		Model
		Board_strain
		Board_signal_num
		Board_prod
		Board_pg_num
		Parameters
		Settings
		Snapshots
		Frequencies
		Neighbors_Frequencies
	}
	data.Model = *model
	data.Board_strain = *model.Board_strain
	data.Board_signal_num = *model.Board_signal_num
	data.Board_prod = *model.Board_prod
	data.Board_pg_num = *model.Board_pg_num
	data.Parameters = model.Parameters
	data.Settings = model.Settings
	data.Snapshots.Snapshot = model.Data_Boards.Snapshots.Snapshot
	data.Frequencies = model.Data_Boards.Frequencies
	data.Neighbors_Frequencies = model.Data_Boards.Neighborhood_Frequencies

}
*/

func (model *Model) save_json() error {
	jsonModel, err := json.Marshal(model)
	if err != nil {
		log.Panicln(err)
	}
	file, err := ioutil.TempFile(".",
		fmt.Sprintf("%v-%v-%v-", model.Settings.DataFilename, model.StartTime, model.GenerationIdx))
	if err != nil {
		log.Panicln(err)
	}
	defer file.Close()
	file.Write(jsonModel)
	return err
}
