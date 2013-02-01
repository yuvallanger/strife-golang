package main

import (
	"fmt"
)

func main() {
	// Reading configuration file
	model := new(Model)
	params, settings := load_config()
	model.Parameters = *params
	model.Settings = *settings

	init_boards(model)

	fmt.Println("Board strain:\n", model.Board_strain)

	run(model)

	fmt.Println("Board strain:\n", model.Board_strain)
	fmt.Println("model.params:\n", model.Parameters)
	fmt.Println("model.settings:\n", model.Settings)
}
