package main

import (
	"fmt"
	"bufio"
	"bytes"
	"io"
	"os"
	"log"
	"strconv"
	"./mtrlTypes"
	"strings"
	"math"
//	"reflect"
	)

func readLines(path string) (lines []string, err error) {
    var (
        file *os.File
        part []byte
        prefix bool
    )
    if file, err = os.Open(path); err != nil {
        return
    }
    defer file.Close()

    reader := bufio.NewReader(file)
    buffer := bytes.NewBuffer(make([]byte, 0))
    for {
        if part, prefix, err = reader.ReadLine(); err != nil {
            break
        }
        buffer.Write(part)
        if !prefix {
            lines = append(lines, buffer.String())
            buffer.Reset()
        }
    }
    if err == io.EOF {
        err = nil
    }
    return
}

func readPoscar(fname string) *mtrlTypes.Model {
	var (
		lines []string
		err error
		scale float64
	)
	m := new(mtrlTypes.Model)
	
	if lines, err = readLines(fname); err != nil {
		log.Fatal(err)
	}

	//header
	m.Header = strings.TrimSpace(lines[0])

	//scale
	if scale, err = strconv.ParseFloat(strings.TrimSpace(lines[1]),len(lines[1])); err != nil {
		log.Fatal(err)
	}

	//Basis	
	m.Basis = make([][]float64,3)
	curr := 2
	m.Ortho = true
	for i:=0;i<3;i++ {
		m.Basis[i] = make([]float64,3)
		c := 0
		for _,v := range strings.Split(lines[i+curr]," ") {
			if len(v)>0 {
				b,_ := strconv.ParseFloat(v,64)
				m.Basis[i][c] = b*scale
			
				if c!=i && m.Basis[i][c]!=0.0 {
					m.Ortho = false
				}
				
				c+=1
			}
		}
	}
	curr+=3

	//Types
	line := strings.Split(lines[curr]," ")
	c := 0
	for _,v := range line{
		if len(v)>0{
			c+=1
		}
	}
	m.Types = make([]int,c)
	m.Natoms = 0
	c = 0
	for _,v := range line{
		if len(v)>0{
			m.Types[c],_ = strconv.Atoi(v)
			m.Natoms += m.Types[c]
			c+=1
		}
	}
	curr+=1

	//Selective Dynamics
	if strings.Contains(lines[curr],"Sel") || strings.Contains(lines[curr],"sel"){
		m.Seldyn = true
		curr+=1
	}
	curr+=1	//Assume line 

	//Atoms
	m.Atoms = make([][]float64,m.Natoms)
	for i,l := range lines[curr : curr+m.Natoms] {
		m.Atoms[i] = make([]float64,3)

		c:=0
		for _,v := range strings.Split(l," ") {
			if len(v)>0 {
				m.Atoms[i][c],_ = strconv.ParseFloat(v,64)
				c+=1
			}
		}

		//Convert from fractional space to real space
		a := m.Atoms[i]
		for c:=0; c<3; c++ {
			m.Atoms[i][c] = m.Basis[0][c]*a[0] + m.Basis[1][c]*a[1] + m.Basis[2][c]*a[2]
		}
	}

	return m
}

func cross3(a []float64, b []float64) []float64{
	c := make([]float64,3)
	c[0] = a[1]*b[2] - a[2]*b[1]
	c[1] = a[2]*b[0] - a[0]*b[2]
	c[2] = a[0]*b[1] - a[1]*b[0]
	return c
}

func dot(a []float64, b []float64) float64{
	var t float64
	for i := range a {
		t += a[i]*b[i]
	}
	return t
}


func volume(vec3 [][]float64) float64{
	a,b,c := vec3[0],vec3[1],vec3[2]
	return math.Abs(dot(a,cross3(b,c)))
}

func RDFperiodic(m *mtrlTypes.Model, cutoff float64, nbins int) ([]float64, []float64) {
	rdist := make([]float64, nbins)
	rbins := make([]float64, nbins)
	
	var c,d,dr,t1,t2,t3,cmin,vol float64
	dr = cutoff/float64(nbins)
	for i:=0; i<nbins; i++ {
		rbins[i] = float64(i)*dr
	}

	b := m.Basis
	ato := m.Atoms
	for i,ai := range ato {
		for _,aj := range ato[i+1:] {

			//Find the minimum image distance of atom pairs
			//ideally should use neighbor list
			cmin = 1E6;
			/*
			for cell:=0;cell<27;cell++{
			//for t1=-1.0; t1<2; t1+=1.0 {
			//for t2=-1.0; t2<2; t2+=1.0 {
			//for t3=-1.0; t3<2; t3+=1.0 {
				t1 = float64(cell/9)    -1
				t2 = float64((cell%9)/3)-1
				t3 = float64(cell%3)    -1
			//	fmt.Println(t1,t2,t3)
				d = ai[0]-aj[0] + t1*b[0][0] + t2*b[1][0] + t3*b[2][0]
				c = d*d
				d = ai[1]-aj[1] + t1*b[0][1] + t2*b[1][1] + t3*b[2][1]
				c += d*d
				d = ai[2]-aj[2] + t1*b[0][2] + t2*b[1][2] + t3*b[2][2]
				c += d*d
				c = math.Sqrt(c)

				if c <= cmin {
					cmin = c
				}
			}//}}*/
		
			for t1=-1.0; t1<2; t1+=1.0 {
			for t2=-1.0; t2<2; t2+=1.0 {
			for t3=-1.0; t3<2; t3+=1.0 {
				d = ai[0]-aj[0] + t1*b[0][0] + t2*b[1][0] + t3*b[2][0]
				c = d*d
				d = ai[1]-aj[1] + t1*b[0][1] + t2*b[1][1] + t3*b[2][1]
				c += d*d
				d = ai[2]-aj[2] + t1*b[0][2] + t2*b[1][2] + t3*b[2][2]
				c += d*d
				c = math.Sqrt(c)

				if c <= cmin {
					cmin = c
				}
			 }}}



			if cmin < cutoff {
				rdist[int(cmin/dr)]+=2;
			}
		}
	}

	Ndensity := float64(m.Natoms*m.Natoms) / volume(m.Basis)
	for i,r := range rbins{
		if i==0 {
			vol = 4.0*math.Pi*dr*dr*dr/3.0
		} else {
			vol = 4.0*math.Pi*r*r*dr
		}
		
		if rdist[i] >= 0.0 {
			rdist[i] /= vol*Ndensity
		}
	}

	return rbins,rdist
}

func main() {
	m := readPoscar("POSCAR_large")
	r,bv := RDFperiodic(m,10.0,1000)
	
	fmt.Println("r g(r)")
	for i := range r {
		fmt.Println(r[i],bv[i])
	}
}