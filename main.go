package main

import (
	"flag"
	"fmt"
	"gitlab.com/yuvallanger/strife-golang/flags"
	"gitlab.com/yuvallanger/strife-golang/sequential"
)

func main() {
	cmdln_flags := flags.Init_flags()
	fmt.Printf("%+v\n", cmdln_flags)

	if cmdln_flags.Sequentialflag {
		sequential.Main(cmdln_flags)
	}

	if cmdln_flags.Concurrentflag {
		concurrent.Main(cmdln_flags)
	}

	flag.Usage()
}
