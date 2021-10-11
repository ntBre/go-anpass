package anpass

// Stat is a type of stationary point
type Stat int

const (
	MAX Stat = iota
	MIN
	STAT
)

func (s Stat) String() string {
	return []string{
		"M A X I M U M",
		"M I N I M U M",
		"S T A T I O N A R Y  P O I N T",
	}[s]
}

