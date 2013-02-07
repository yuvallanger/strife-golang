package main

import (
	"flag"
	c "github.com/yuvallanger/game-of-strife/golang/src/strife/concurrent"
	s "github.com/yuvallanger/game-of-strife/golang/src/strife/sequential"
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
