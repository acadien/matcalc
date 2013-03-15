package mtrlTypes

type Model struct {
	Header string
	Basis [][]float64
	Ortho bool
	Natoms int
	Types []int
	Seldyn bool
	Atoms [][]float64
}

