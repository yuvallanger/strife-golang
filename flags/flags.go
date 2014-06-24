package flags

import (
	"flag"
	"log"
)

type Flags struct {
	Sequentialflag         bool
	Concurrentflag         bool
	Cpuprofileflag         string
	Imagesflag             bool
	Settings_filename_flag string
	Czaranflag             bool
	Avigdorflag            bool
}

func Init_flags() (flags Flags) {
	const (
		defaultSequential = false
		usageSequential   = "Run sequential model."

		defaultConcurrent = false
		usageConcurrent   = "Run concurrent model."

		defaultCpuprofile = ""
		usageCpuprofile   = "Write cpu profile to file."

		defaultImages = false
		usageImages   = "Output snapshots as images."

		defaultSettings = ""
		usageSettings   = "Settings filename."

		defaultAvigdor = false
		usageAvigdor   = "Run Avigdor's model."

		defaultCzaran = false
		usageCzaran   = "Run Czaran's model."
	)
	flag.BoolVar(&flags.Sequentialflag, "sequential", defaultSequential, usageSequential)
	flag.BoolVar(&flags.Concurrentflag, "concurrent", defaultConcurrent, usageConcurrent)
	flag.StringVar(&flags.Cpuprofileflag, "cpuprofile", defaultCpuprofile, "")
	flag.BoolVar(&flags.Imagesflag, "images", defaultImages, usageImages)
	flag.StringVar(&flags.Settings_filename_flag, "settings", defaultSettings, usageSettings)
	flag.BoolVar(&flags.Avigdorflag, "avigdor", defaultAvigdor, usageAvigdor)
	flag.BoolVar(&flags.Czaranflag, "czaran", defaultCzaran, usageCzaran)

	flag.Parse()

	if (flags.Avigdorflag && flags.Czaranflag) || !(flags.Avigdorflag || flags.Czaranflag) {
		log.Fatalf("Must include one and only one of the flags -avigdor and -czaran")
		flag.Usage()
	}
	if flags.Concurrentflag && flags.Sequentialflag {
		log.Fatalf("Must include one and only one of the flags -sequential and -concurrent")
		flag.Usage()
	}
	return
}
