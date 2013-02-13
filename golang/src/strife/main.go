package main

import (
	"flag"
	"fmt"
	"log"
	"strife/concurrent"
	"strife/sequential"
)

var sequentialflag = flag.Bool("sequential", true, "Run sequential model.")
var concurrentflag = flag.Bool("concur", false, "Run concurrent model.")
var cpuprofile = flag.String("cpuprofile", "", "Write cpu profile to file.")
var imagesflag = flag.Bool("images", false, "Output snapshots as images.")

func main() {
	flag.Parse()
	fmt.Println("Just to let you know, here are some command line flags:")
	flag.Usage()

	if *concurrentflag && *sequentialflag {
		flag.Usage()
		log.Fatalf("Should only use one of the two.")
	}
	if *sequentialflag {
		sequential.Main(cpuprofile, imagesflag)
	}
	if *concurrentflag {
		concurrent.Main(cpuprofile)
	}
	flag.Usage()
}
