package main

import (
	"flag"
	"strife/concurrent"
	"strife/sequential"
)

func main() {
	modeltype := flag.String("algo", "sequential", "The model to run.")
	flag.Parse()
	switch *modeltype {
	case "sequential":
		sequential.Main()
	case "concurrent":
		concurrent.Main()
	}
}
