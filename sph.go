// Smoothed-Particle Hydrodynamics demo (Position-Based Fluids variant)
// --------------------------------------------------------------------
// Single-file implementation with simple axis-aligned box boundaries.
// Wall-friction and stronger XSPH viscosity have been added to stop the
// bottom layer of particles from skating back and forth indefinitely.

package main

import (
	"bufio"
	"fmt"
	"math"
	"math/rand"
	"os"
	"time"
)

/*-------------------------------------------------------------------------*/
/*  Basic types                                                            */
/*-------------------------------------------------------------------------*/

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

/*-------------------------------------------------------------------------*/
/*  Simulation constants                                                   */
/*-------------------------------------------------------------------------*/

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

/*-------------------------------------------------------------------------*/
/*  Bounding box with wall-friction                                        */
/*-------------------------------------------------------------------------*/

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

/*-------------------------------------------------------------------------*/
/*  Helpers                                                                */
/*-------------------------------------------------------------------------*/

func distSq(a, b vector2) float64 {
	dx, dy := a.x-b.x, a.y-b.y
	return dx*dx + dy*dy
}
func dist(a, b vector2) float64 { return math.Sqrt(distSq(a, b)) }

func unitVec(a, b vector2) vector2 {
	d := a.minus(b)
	return d.divideSafeF(d.mag())
}

/*-------------------------------------------------------------------------*/
/*  Particle methods                                                       */
/*-------------------------------------------------------------------------*/

func (p *particle) updateVelocity(dt float64) {
	diff := p.position.minus(p.previousPosition)
	p.velocity = diff.divideSafeF(dt)
}
func (p *particle) integrate(dt float64) {
	p.previousPosition = p.position
	p.position.add(p.velocity.multiplyF(dt))
}

/*-------------------------------------------------------------------------*/
/*  External forces                                                        */
/*-------------------------------------------------------------------------*/

func applyExternalForces(p *particle, c constants) {
	p.velocity.add(c.g.multiplyF(c.timeStep))
}

/*-------------------------------------------------------------------------*/
/*  Double-density relaxation (PBF)                                        */
/*-------------------------------------------------------------------------*/

func doubleDensity(c constants, i int, neigh []int, particles []particle) {
	var density, nearDensity float64

	for _, j := range neigh {
		if i == j {
			continue
		}
		q := dist(particles[i].position, particles[j].position) / c.h
		if q < 1.0 {
			density += (1 - q) * (1 - q)
			nearDensity += (1 - q) * (1 - q) * (1 - q)
		}
	}

	pressure := c.k * (density - c.density0)
	nearPressure := c.kNear * nearDensity

	delta := vector2{}
	for _, j := range neigh {
		if i == j {
			continue
		}
		q := dist(particles[i].position, particles[j].position) / c.h
		if q >= 1.0 {
			continue
		}

		dir := unitVec(particles[i].position, particles[j].position)
		D := dir.multiplyF(c.timeStep * c.timeStep *
			(pressure*(1-q) + nearPressure*(1-q)*(1-q)))

		if m := D.mag(); m > 0.5*c.h {
			D = D.divideF(m).multiplyF(0.5 * c.h) // clamp
		}

		particles[j].position.add(D.multiplyF(0.5))
		delta.subtract(D.multiplyF(0.5))
	}
	particles[i].position.add(delta)
}

/*-------------------------------------------------------------------------*/
/*  XSPH viscosity                                                         */
/*-------------------------------------------------------------------------*/

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

/*-------------------------------------------------------------------------*/
/*  Neighbour search (naïve O(N²))                                         */
/*-------------------------------------------------------------------------*/

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

/*-------------------------------------------------------------------------*/
/*  Simulation step                                                        */
/*-------------------------------------------------------------------------*/

func update(p []particle, bb *boundingBox, c constants) {

	// 1  external forces
	for i := range p {
		applyExternalForces(&p[i], c)
	}

	// 2  viscosity
	neigh := neighbours(p, c.h)
	viscosity(c, p, neigh)

	// 3  predict positions
	for i := range p {
		p[i].integrate(c.timeStep)
	}

	// 4  density / pressure (4 passes, neighbour list refreshed)
	for pass := 0; pass < 4; pass++ {
		neigh = neighbours(p, c.h)
		for i := range p {
			doubleDensity(c, i, neigh[i], p)
		}
	}

	// 5  collisions & velocity update
	for i := range p {
		bb.collide(&p[i], c)
		p[i].updateVelocity(c.timeStep)
	}
}

/*-------------------------------------------------------------------------*/
/*  I/O helper                                                             */
/*-------------------------------------------------------------------------*/

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

/*-------------------------------------------------------------------------*/
/*  Main                                                                    */
/*-------------------------------------------------------------------------*/

func main() {
	rand.Seed(time.Now().UnixNano())

	c := constants{
		g:            vector2{0, -9.81},
		epsilon:      1e-3,
		damping:      0.25,
		h:            0.6,
		k:            0.005,
		kNear:        0.020,
		density0:     4.0,
		timeStep:     0.001,
		sigma:        0.10, // stronger linear viscosity
		beta:         0.00, // no quadratic term
		wallFriction: 0.5,  // lose 50 % tangential speed at walls
	}

	bb := boundingBox{top: 10, bottom: 0, left: 0, right: 10}

	/*  initial particle block  */
	N := 500
	particles := make([]particle, N)
	for i := 0; i < N; i++ {
		particles[i].position = vector2{x: math.Mod(1+float64(i)*0.5, 10), y: 8}
		particles[i].velocity = vector2{x: 2 * math.Sin(float64(i)), y: -10}
		particles[i].previousPosition = particles[i].position
	}

	/*  run  */
	steps := 10000
	if err := os.MkdirAll("out", 0o755); err != nil && !os.IsExist(err) {
		panic(err)
	}
	for s := 0; s < steps; s++ {
		update(particles, &bb, c)
		if s%10 == 0 {
			if err := savePositions(fmt.Sprintf("out/positions_%04d.csv", s), particles); err != nil {
				panic(err)
			}
		}
	}
}
