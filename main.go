package main

import "fmt"
import "image"
import "image/color"
import "image/png"
import "image/draw"
import _ "image/jpeg"
import "os"
import "github.com/golang/geo/r3"
import "math"

const (
	XAxis r3.Axis = iota
	YAxis
	ZAxis
)

type Luz struct {
	Posicao  r3.Vector
	Intensidade float64
}

type Materiais struct {
	IndiceDeRefracao  float64
	Albedo           []float64
	Cor    color.RGBA
	ComponenteEspecular float64
}

type Esfera struct {
	Centro   r3.Vector
	Raio   float64
	Material Materiais
}



func (s Esfera) rayIntersect(orig, dir r3.Vector, t0 *float64) bool {
	L := r3.Vector.Sub(s.Centro, orig)
	tca := r3.Vector.Dot(L, dir)
	d2 := r3.Vector.Dot(L, L) - tca*tca

	if d2 > s.Raio*s.Raio {
		return false
	}

	thc := math.Sqrt(s.Raio*s.Raio - d2)
	*t0 = tca - thc
	t1 := tca + thc

	if *t0 < 0.0 {
		*t0 = t1
	}
	if *t0 < 0.0 {
		return false
	}
	return true
}

func sceneIntersect(orig, dir r3.Vector, esferas []Esfera, hit, N *r3.Vector, material *Materiais) bool {
	esferasDist := math.MaxFloat64
	distanciaPlano := math.MaxFloat64

	for i := 0; i < len(esferas); i++ {
		distI := 0.0
		if esferas[i].rayIntersect(orig, dir, &distI) && distI < esferasDist {
			esferasDist = distI
			*hit = r3.Vector.Add(orig, r3.Vector.Mul(dir, distI))
			*N = r3.Vector.Normalize(r3.Vector.Sub(*hit, esferas[i].Centro))
			*material = esferas[i].Material
			_ = material
		}
	}
	if math.Abs(dir.Y) > 0.001 {
		d := -(orig.Y + 4.) / dir.Y // the checkerboard plane has equation y = -4
		pt := r3.Vector.Add(orig, r3.Vector.Mul(dir, d))

		if d > 0 && math.Abs(pt.X) < 10. && pt.Z < -10. && pt.Z > -30. && d < esferasDist {
			distanciaPlano = d
			*hit = pt
			*N = r3.Vector{0, 1, 0}
			if (int(.5*hit.X+1000)+int(.5*hit.Z))&1 == 1 {
				// need to re-define 'material', to a valid one, because it is local
				// to the `if` statement above, and we get a null-panic like, eg. when accessing albedo
				*material = esferas[0].Material
				material.Cor = color.RGBA{255, 255, 255, 255}
			} else {
				*material = esferas[0].Material
				material.Cor = color.RGBA{255, 177, 65, 255}
			}
		}
	}

	return min(esferasDist, distanciaPlano) < 1000

}

func castRay(orig, dir r3.Vector, esferas []Esfera, luzes []Luz, depth int, envmap image.Image) color.RGBA {
	var point, N r3.Vector
	var material Materiais
	depth += 1
	if depth > 10 || !sceneIntersect(orig, dir, esferas, &point, &N, &material) { // depth para reflexões
		// aqui normalmente teriamos um background static, substituido pelo env map
		normalized := r3.Vector.Normalize(dir)
		bounds := envmap.Bounds()

		r, g, b, _ := envmap.At(int((normalized.X/2+0.5)*float64(bounds.Max.X)), int((-normalized.Y/2+0.5)*float64(bounds.Max.Y))).RGBA()
		return color.RGBA{uint8(r >> 8), uint8(g >> 8), uint8(b >> 8), 255}
	}

	// começa reflexão e refração

	var reflexaoOrig r3.Vector
	var refracaoOrig r3.Vector

	reflexaoDirecao := r3.Vector.Normalize(Reflect(dir, N))
	refracaoDirecao := r3.Vector.Normalize(Refract(dir, N, material.IndiceDeRefracao))

	if r3.Vector.Dot(reflexaoDirecao, N) < 0 {
		reflexaoOrig = r3.Vector.Sub(point, r3.Vector.Mul(N, 0.001))
	} else {
		reflexaoOrig = r3.Vector.Add(point, r3.Vector.Mul(N, 0.001))
	}

	if r3.Vector.Dot(refracaoDirecao, N) < 0 {
		refracaoOrig = r3.Vector.Sub(point, r3.Vector.Mul(N, 0.001))
	} else {
		refracaoOrig = r3.Vector.Add(point, r3.Vector.Mul(N, 0.001))
	}

	CorDeReflaxao := castRay(reflexaoOrig, reflexaoDirecao, esferas, luzes, depth, envmap)
	CorDeRefracao := castRay(refracaoOrig, refracaoDirecao, esferas, luzes, depth, envmap)

	// termina reflexão e refração

	// iluminação difusa e especular começa

	intensidadeLuzDifusa, intensidadeLuzEspecular := 0.0, 0.0
	for i := 0; i < len(luzes); i++ {
		LuzDirecao := r3.Vector.Normalize(r3.Vector.Sub(luzes[i].Posicao, point))

		// sombra começa
		LuzDistancia := r3.Vector.Norm(r3.Vector.Sub(luzes[i].Posicao, point))

		var sombraOrig r3.Vector
		if r3.Vector.Dot(LuzDirecao, N) < 0 {
			sombraOrig = r3.Vector.Sub(point, r3.Vector.Mul(N, 0.001))
		} else {
			sombraOrig = r3.Vector.Add(point, r3.Vector.Mul(N, 0.001))
		}

		var sombraPT, sombraN r3.Vector
		var tmpmaterial Materiais

		if sceneIntersect(sombraOrig, LuzDirecao, esferas, &sombraPT, &sombraN, &tmpmaterial) && r3.Vector.Norm(r3.Vector.Sub(sombraPT, sombraOrig)) < LuzDistancia {
			continue
		}
		// sombra termina

		// produtos escalar para calculo das intensidade da luz
		intensidadeLuzDifusa += luzes[i].Intensidade * max(0, r3.Vector.Dot(LuzDirecao, N))
		// phong
		mLuzDirecao := r3.Vector.Mul(LuzDirecao, -1.)
		intensidadeLuzEspecular += math.Pow(max(0., r3.Vector.Dot(r3.Vector.Mul(Reflect(mLuzDirecao, N), -1), dir)), material.ComponenteEspecular) * luzes[i].Intensidade
	}

	// iluminação difusa e especular termina

	// apos os 4 calculos (iluminação difusa, especular, reflexão e refração) temos as cores finais
	res1x := float64(material.Cor.R) * intensidadeLuzDifusa * material.Albedo[0]
	res1y := float64(material.Cor.G) * intensidadeLuzDifusa * material.Albedo[0]
	res1z := float64(material.Cor.B) * intensidadeLuzDifusa * material.Albedo[0]

	black := color.RGBA{255, 255, 255, 255}
	res2x := float64(black.R) * intensidadeLuzEspecular * material.Albedo[1]
	res2y := float64(black.G) * intensidadeLuzEspecular * material.Albedo[1]
	res2z := float64(black.B) * intensidadeLuzEspecular * material.Albedo[1]

	res3x := float64(CorDeReflaxao.R) * material.Albedo[2]
	res3y := float64(CorDeReflaxao.G) * material.Albedo[2]
	res3z := float64(CorDeReflaxao.B) * material.Albedo[2]

	res4x := float64(CorDeRefracao.R) * material.Albedo[3]
	res4y := float64(CorDeRefracao.G) * material.Albedo[3]
	res4z := float64(CorDeRefracao.B) * material.Albedo[3]

	return AddColors(r3.Vector{res1x, res1y, res1z}, r3.Vector{res2x, res2y, res2z}, r3.Vector{res3x, res3y, res3z}, r3.Vector{res4x, res4y, res4z})
}

func render(esferas []Esfera, luzes []Luz, envmap image.Image) {

	const w = 1024.0
	const h = 768.0
	const fov = math.Pi / 2

	img := image.NewRGBA(image.Rect(0, 0, w, h))
	for j := 0; j < h; j++ {
		for i := 0; i < w; i++ {
			// trajetoria inicial de cada ray
			var x, y float64
			x = (2.0*(float64(i)+0.5)/w - 1.0) * math.Tan(fov/2.0) * w / h
			y = -(2.0*(float64(j)+0.5)/h - 1.0) * math.Tan(fov/2.0)
			dir := r3.Vector.Normalize(r3.Vector{x, y, -1.0})
			img.Set(i, j, castRay(r3.Vector{0, 0, 0}, dir, esferas, luzes, 0, envmap))
		}
	}

	f, _ := os.OpenFile("out.png", os.O_WRONLY|os.O_CREATE, 0600)
	defer f.Close()
	png.Encode(f, img)

	fmt.Println("Done!")
}

func Reflect(I, N r3.Vector) r3.Vector {
	return r3.Vector.Sub(I, r3.Vector.Mul(N, 2.0*r3.Vector.Dot(I, N)))
}

// Lei de Snell
func Refract(I, N r3.Vector, refractiveIdx float64) r3.Vector {
	cosi := -max(-1., min(1, r3.Vector.Dot(I, N)))
	etai, etat := 1., refractiveIdx
	n := N
	if cosi < 0. {
		cosi = -cosi
		etai, etat = swap(etai, etat)
		n = r3.Vector.Mul(N, -1)
	}
	eta := etai / etat
	k := 1. - eta*eta*(1.-cosi*cosi)
	if k < 0. {
		return r3.Vector{0., 0., 0.}
	} else {
		return r3.Vector.Add(r3.Vector.Mul(I, eta), r3.Vector.Mul(n, (eta*cosi-math.Sqrt(k))))
	}
}

func AddColors(i, j, k, l r3.Vector) color.RGBA {
	r, g, b := (i.X + j.X + k.X + l.X), (i.Y + j.Y + k.Y + l.Y), (i.Z + j.Z + k.Z + l.Z)
	maxc := float64(max(float64(r), max(float64(g), float64(b))))
	if maxc > 255. {
		return color.RGBA{uint8(float64(r) * 255. / maxc),
			uint8(float64(g) * 255. / maxc),
			uint8(float64(b) * 255. / maxc),
			255}
	}
	return color.RGBA{uint8(r),
		uint8(g),
		uint8(b),
		255}
}

func main() {

	marfim := Materiais{IndiceDeRefracao: 1.0, Albedo: []float64{0.3, 0.6, 0.1, 0.0}, Cor: color.RGBA{100, 100, 75, 255}, ComponenteEspecular: 50.}
	vidro := Materiais{IndiceDeRefracao: 1.5, Albedo: []float64{0.0, 0.5, 0.1, 0.8}, Cor: color.RGBA{255, 255, 255, 255}, ComponenteEspecular: 1425.}
	borracha := Materiais{IndiceDeRefracao: 1.0, Albedo: []float64{0.9, 0.1, 0.0, 0.0}, Cor: color.RGBA{76, 25, 25, 255}, ComponenteEspecular: 10.}
	espelho := Materiais{IndiceDeRefracao: 1.0, Albedo: []float64{0.0, 10.0, 0.8, 0.0}, Cor: color.RGBA{255, 255, 255, 255}, ComponenteEspecular: 1425.}

	esferas := make([]Esfera, 0)
	esferas = append(esferas, Esfera{r3.Vector{-3.0, 0.0, -16.0}, 0, marfim})
	esferas = append(esferas, Esfera{r3.Vector{-4.0, -1.5, -12.0}, 2, vidro})
	esferas = append(esferas, Esfera{r3.Vector{1.5, -0.5, -18.0}, 3, borracha})
	esferas = append(esferas, Esfera{r3.Vector{11.0, 5.0, -18.0}, 5, espelho})

	luzes := make([]Luz, 0)
	luzes = append(luzes, Luz{r3.Vector{-20, 20, 20}, 1.5})
	luzes = append(luzes, Luz{r3.Vector{30, 50, -25}, 1.8})
	luzes = append(luzes, Luz{r3.Vector{30, 20, 30}, 1.7})

	// Loading evnrioment map
	f, err := os.Open("envmap.jpg")
	if err != nil {
		panic(err)
	}
	defer f.Close()

	img, format, err := image.Decode(f)
	if err != nil {
		panic(err)
	}

	size := img.Bounds()
	skybox := image.NewRGBA(size)
	draw.Draw(skybox, size, img, size.Min, draw.Src)
	_ = img
	_ = format
	// Env map 

	render(esferas, luzes, skybox)
}

// funções auxiliares
func max(a, b float64) float64 {
    if a >= b {
        return a
    }
    return b
}
func swap(a, b float64)(float64, float64) {
    return b, a
}
func min(a, b float64) float64 {
    if a < b {
        return a
    }
    return b
}