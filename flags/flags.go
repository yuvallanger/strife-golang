package flags

import (
	"flag"
)

type Flags struct {
	Cpuprofileflag         string
	Settings_filename_flag string
}

func Init_flags() (flags Flags) {
	const (
		defaultCpuprofile = ""
		usageCpuprofile   = "Write cpu profile to file."

		defaultSettings = ""
		usageSettings   = "Settings filename."
	)
	flag.StringVar(&flags.Cpuprofileflag, "cpuprofile", defaultCpuprofile, "")
	flag.StringVar(&flags.Settings_filename_flag, "settings", defaultSettings, usageSettings)

	flag.Parse()
	return
}
