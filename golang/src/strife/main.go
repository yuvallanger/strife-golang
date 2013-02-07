package main

import (
	"flag"
	c "strife/concurrent"
	s "strife/sequential"
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
