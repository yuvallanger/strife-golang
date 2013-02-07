package miscow

import "fmt"

func MyMod(a, b int) int {
	return ((a % b) + b) % b
}

func Iterrange(lo, hi int) func() (int, bool) {
	i := lo
	return func() (int, bool) {
		var end bool
		if i < hi {
			end = true
		} else {
			end = false
		}
		return i, end
	}
}

func Trace(s string)   { fmt.Println("trace: entering:", s) }
func Untrace(s string) { fmt.Println("trace: leaving:", s) }

