package main

import (
	"bufio"
	"fmt"
	"math"
	"math/rand"
	"os"
	"time"
)

type vector2 struct{ x, y float64 }

func (v vector2) mag() float64                { return math.Sqrt(v.x*v.x + v.y*v.y) }
func (v vector2) multiply(w vector2) vector2  { return vector2{v.x * w.x, v.y * w.y} }
func (v vector2) divide(w vector2) vector2    { return vector2{v.x / w.x, v.y / w.y} }
func (v vector2) multiplyF(s float64) vector2 { return vector2{v.x * s, v.y * s} }
func (v vector2) divideF(s float64) vector2   { return vector2{v.x / s, v.y / s} }
func (v vector2) minus(w vector2) vector2     { return vector2{v.x - w.x, v.y - w.y} }
func (v vector2) plus(w vector2) vector2      { return vector2{v.x + w.x, v.y + w.y} }
func (v *vector2) subtract(w vector2)         { v.x -= w.x; v.y -= w.y }
func (v *vector2) add(w vector2)              { v.x += w.x; v.y += w.y }
func (v *vector2) divideSafeF(s float64) vector2 { // same as divideF but guarded
	const eps = 1e-9
	if math.Abs(s) < eps {
		return vector2{}
	}
	return v.divideF(s)
}

type particle struct {
	position, velocity, previousPosition vector2
}

type constants struct {
	g vector2

	epsilon  float64
	damping  float64
	h        float64
	k        float64
	kNear    float64
	density0 float64
	timeStep float64
	sigma    float64
	beta     float64

	wallFriction float64 // 0 → stick, 1 → fully slippery
}

// ------------------------------------------------------------------
// 2-D Poly6 kernel and Spiky gradient (monomial form)
// W_poly6(r,h)  = (4/π h⁸) (h²−r²)³           for 0 ≤ r ≤ h
// ∇W_spiky      = −(10/π h⁵) (h−r)²  r̂        for 0 ≤ r ≤ h
// ------------------------------------------------------------------

func poly6(r2, h float64) float64 {
	if r2 >= h*h {
		return 0
	}
	inv := (h*h - r2)
	return (4.0 / (math.Pi * math.Pow(h, 8))) * inv * inv * inv
}

func spikyGrad(r float64, dir vector2, h float64) vector2 {
	if r == 0 || r >= h {
		return vector2{}
	}
	f := -10.0 / (math.Pi * math.Pow(h, 5)) * (h - r) * (h - r)
	return dir.multiplyF(f)
}

type boundingBox struct{ top, left, bottom, right float64 }

func (bb *boundingBox) collide(p *particle, c constants) {
	switch {
	case p.position.y < bb.bottom: // floor
		p.position.y = bb.bottom + c.epsilon
		p.velocity.y = -p.velocity.y * (1 - c.damping)
		p.velocity.x *= c.wallFriction

	case p.position.y > bb.top: // ceiling
		p.position.y = bb.top - c.epsilon
		p.velocity.y = -p.velocity.y * (1 - c.damping)
		p.velocity.x *= c.wallFriction

	case p.position.x < bb.left: // left wall
		p.position.x = bb.left + c.epsilon
		p.velocity.x = -p.velocity.x * (1 - c.damping)
		p.velocity.y *= c.wallFriction

	case p.position.x > bb.right: // right wall
		p.position.x = bb.right - c.epsilon
		p.velocity.x = -p.velocity.x * (1 - c.damping)
		p.velocity.y *= c.wallFriction
	}

	// keep previousPosition consistent so updateVelocity sees the bounce
	p.previousPosition = vector2{
		x: p.position.x - p.velocity.x*c.timeStep,
		y: p.position.y - p.velocity.y*c.timeStep,
	}
}

func distSq(a, b vector2) float64 {
	dx, dy := a.x-b.x, a.y-b.y
	return dx*dx + dy*dy
}
func dist(a, b vector2) float64 { return math.Sqrt(distSq(a, b)) }

func unitVec(a, b vector2) vector2 {
	d := a.minus(b)
	return d.divideSafeF(d.mag())
}

func (p *particle) updateVelocity(dt float64) {
	diff := p.position.minus(p.previousPosition)
	p.velocity = diff.divideSafeF(dt)
}
func (p *particle) integrate(dt float64) {
	p.previousPosition = p.position
	p.position.add(p.velocity.multiplyF(dt))
}

func viscosity(c constants, particles []particle, neigh [][]int) {
	for i := range neigh {
		for _, j := range neigh[i] {
			if j >= i {
				break
			}
			q := dist(particles[i].position, particles[j].position) / c.h
			if q >= 1 {
				continue
			}

			relVel := particles[i].velocity.minus(particles[j].velocity)
			unit := unitVec(particles[i].position, particles[j].position)
			rel := relVel.multiply(unit)
			mag := rel.mag()
			if mag == 0 {
				continue
			}

			V := unit.multiplyF(c.timeStep * (1 - q) * (c.sigma*mag + c.beta*mag*mag))
			particles[i].velocity.add(V.multiplyF(0.5))
			particles[j].velocity.subtract(V.multiplyF(0.5))
		}
	}
}

func neighbours(particles []particle, h float64) [][]int {
	N := len(particles)
	out := make([][]int, N)
	h2 := h * h
	for i := 0; i < N; i++ {
		for j := i + 1; j < N; j++ {
			if distSq(particles[i].position, particles[j].position) < h2 {
				out[i] = append(out[i], j)
				out[j] = append(out[j], i)
			}
		}
	}
	return out
}

func sphPressureAccel(acc []vector2, p []particle, neigh [][]int, c constants, m float64) {
	N := len(p)
	rho := make([]float64, N)
	P := make([]float64, N)

	for i := 0; i < N; i++ {
		rho[i] += m * poly6(0, c.h)
		for _, j := range neigh[i] {
			r2 := distSq(p[i].position, p[j].position)
			w := m * poly6(r2, c.h)
			rho[i] += w
			rho[j] += w
		}
	}

	/* pressure from Tait EOS (simple ideal gas) */
	for i := 0; i < N; i++ {
		P[i] = c.k * (rho[i] - c.density0)
	}

	for i := 0; i < N; i++ {
		for _, j := range neigh[i] {
			if j <= i {
				continue
			}
			r := dist(p[i].position, p[j].position)
			dir := unitVec(p[i].position, p[j].position)
			grad := spikyGrad(r, dir, c.h)

			f := grad.multiplyF(m * (P[i]/(rho[i]*rho[i]) + P[j]/(rho[j]*rho[j])))
			acc[i].subtract(f.divideF(m)) // a_i -= f / m
			acc[j].add(f.divideF(m))      // a_j += f / m
		}
	}
}

func update(p []particle, bb *boundingBox, c constants) {

	neigh := neighbours(p, c.h)

	acc := make([]vector2, len(p))

	for i := range p {
		acc[i] = c.g
	}

	sphPressureAccel(acc, p, neigh, c, 1.0) // mass = 1

	viscosity(c, p, neigh)

	for i := range p {
		p[i].velocity.add(acc[i].multiplyF(c.timeStep))
		p[i].integrate(c.timeStep)
		bb.collide(&p[i], c)
		p[i].updateVelocity(c.timeStep)
	}
}

func savePositions(filename string, particles []particle) error {
	f, err := os.Create(filename)
	if err != nil {
		return err
	}
	defer f.Close()

	w := bufio.NewWriter(f)
	for i, p := range particles {
		fmt.Fprintf(w, "%d,%.6f,%.6f\n", i, p.position.x, p.position.y)
	}
	return w.Flush()
}

func main() {
	rand.Seed(time.Now().UnixNano())

	c := constants{
		g:            vector2{0, -9.81},
		epsilon:      1e-3,
		damping:      0.25,
		h:            0.6,
		density0:     4.0,
		k:            30.0,
		timeStep:     0.0005,
		sigma:        0.04,
		beta:         0.00,
		wallFriction: 0.5,
	}

	bb := boundingBox{top: 10, bottom: 0, left: 0, right: 10}

	N := 1000
	particles := make([]particle, N)
	for i := 0; i < N; i++ {
		particles[i].position = vector2{x: math.Mod(1+float64(i)*0.5, 10), y: 8}
		particles[i].velocity = vector2{x: 2 * math.Sin(float64(i)), y: -10}
		particles[i].previousPosition = particles[i].position
	}

	steps := 20000
	if err := os.MkdirAll("out", 0o755); err != nil && !os.IsExist(err) {
		panic(err)
	}
	for s := 0; s < steps; s++ {
		update(particles, &bb, c)
		if s%10 == 0 {
			if err := savePositions(fmt.Sprintf("out/positions_%06d.csv", s), particles); err != nil {
				panic(err)
			}
		}
	}
}
