package main

import (
	"flag"
	s "strife/sequential"
    c "strife/concurrent"
)

func main() {
	modeltype := flag.String("algo", "sequential", "The model to run.")
	flag.Parse()
	switch *modeltype {
		case "sequential":
			s.Main()
	case "concurrent":
		c.Main()
	}
}
