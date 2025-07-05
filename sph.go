package main

import (
	"bufio"
	"fmt"
	"math"
	"math/rand"
	"os"
	"time"
)

type vector2 struct {
	x float64
	y float64
}

type particle struct {
	position         vector2
	velocity         vector2
	previousPosition vector2
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
}

type boundingBox struct {
	top    float64
	left   float64
	bottom float64
	right  float64
}

func (bb *boundingBox) collide(p *particle, c constants) {

	if p.position.y < bb.bottom {

		p.position.y = bb.bottom + c.epsilon
		p.velocity.y = -p.velocity.y * (1 - c.damping)

	} else if p.position.y > bb.top {

		p.position.y = bb.top - c.epsilon
		p.velocity.y = -p.velocity.y * (1 - c.damping)

	} else if p.position.x < bb.left {

		p.position.x = bb.left + c.epsilon
		p.velocity.x = -p.velocity.x * (1 - c.damping)

	} else if p.position.x > bb.right {

		p.position.x = bb.right - c.epsilon
		p.velocity.x = -p.velocity.x * (1 - c.damping)
	}

	p.previousPosition = vector2{
		x: p.position.x - p.velocity.x*c.timeStep,
		y: p.position.y - p.velocity.y*c.timeStep,
	}

}

func (v *vector2) mag() float64 {
	return math.Sqrt(v.x*v.x + v.y*v.y)
}

func (v1 *vector2) multiply(v2 vector2) vector2 {
	return vector2{x: v1.x * v2.x, y: v1.y * v2.y}
}

func (v1 *vector2) divide(v2 vector2) vector2 {
	return vector2{x: v1.x / v2.x, y: v1.y / v2.y}
}

func (v1 *vector2) multiplyF(p float64) vector2 {
	return vector2{x: v1.x * p, y: v1.y * p}
}

func (v1 *vector2) divideF(p float64) vector2 {
	const eps = 1e-9
	if math.Abs(p) < eps { // avoid blowing up
		return vector2{x: 0, y: 0}
	}
	return vector2{x: v1.x / p, y: v1.y / p}
}

func (v1 *vector2) minus(v2 vector2) vector2 {
	return vector2{x: v1.x - v2.x, y: v1.y - v2.y}
}

func (v1 *vector2) plus(v2 vector2) vector2 {
	return vector2{x: v1.x + v2.x, y: v1.y + v2.y}
}

func (v1 *vector2) subtract(v2 vector2) {
	v1.x -= v2.x
	v1.y -= v2.y
}

func (v1 *vector2) add(v2 vector2) {
	v1.x += v2.x
	v1.y += v2.y
}

func (v1 *vector2) subtractF(s float64) {
	v1.x -= s
	v1.y -= s
}

func (v1 *vector2) addF(s float64) {
	v1.x += s
	v1.y += s
}

func distSq(p1, p2 vector2) float64 {
	dx := p1.x - p2.x
	dy := p1.y - p2.y
	return dx*dx + dy*dy
}

func dist(p1 vector2, p2 vector2) float64 {
	return math.Sqrt(distSq(p1, p2))
}

func unitVec(p1, p2 vector2) vector2 {
	d := p1.minus(p2)
	length := d.mag()
	if length < 1e-9 {
		return vector2{x: 0, y: 0}
	}
	return d.divideF(length)
}

func (p *particle) updateVelocity(timeStep float64) {
	diffPos := p.position.minus(p.previousPosition)
	p.velocity = diffPos.divideF(timeStep)
}

func (p *particle) integrate(timeStep float64) {
	p.previousPosition = p.position
	p.position.add(p.velocity.multiplyF(timeStep))
}

func applyExternalForces(p *particle, c constants) {

	// Gravity
	p.velocity.add(c.g.multiplyF(c.timeStep))
}

func doubleDensity(c constants, this int, neighbours []int, particles []particle) {
	density := 0.0
	nearDensity := 0.0

	for _, neighbour := range neighbours {
		if this == neighbour {
			continue
		}
		q := dist(particles[this].position, particles[neighbour].position) / c.h
		if q < 1.0 {
			density += (1.0 - q) * (1.0 - q)
			nearDensity += (1.0 - q) * (1.0 - q) * (1.0 - q)
		}
	}

	pressure := c.k * (density - c.density0)
	nearPressure := c.kNear * nearDensity

	deltaX := vector2{x: 0.0, y: 0.0}

	for _, neighbour := range neighbours {
		if this == neighbour {
			continue
		}
		q := dist(particles[this].position, particles[neighbour].position) / c.h
		if q < 1.0 { // was “> 1.0”
			pressureTerm := pressure * (1.0 - q)
			nearPressureTerm := nearPressure * (1.0 - q) * (1.0 - q)

			d := dist(particles[this].position, particles[neighbour].position)
			if d < 1e-12 { // 1e-12 ≈ machine epsilon for metre-scale coords
				continue // skip this interaction
			}

			D := unitVec(particles[this].position, particles[neighbour].position)
			D = D.multiplyF(c.timeStep * c.timeStep * (pressureTerm + nearPressureTerm))

			if mag := D.mag(); mag > c.h*0.5 { // ← new safety clamp
				D = D.divideF(mag)
				D = D.multiplyF(c.h * 0.5) //   max 0.5 h per step
			}

			particles[neighbour].position.add(D.multiplyF(0.5))
			deltaX.subtract(D.multiplyF(0.5))

		} // q > 1.0?
	} // End neighbour loop

	particles[this].position.add(deltaX)

}

func viscosity(c constants, particles []particle, neighboursArray [][]int) {

	for i, _ := range neighboursArray {
		for _, j := range neighboursArray[i] {

			if j >= i {
				break
			}

			q := dist(particles[i].position, particles[j].position) / c.h

			if q < 1 {

				d := dist(particles[i].position, particles[j].position)
				if d < 1e-12 { // 1e-12 ≈ machine epsilon for metre-scale coords
					continue // skip this interaction
				}

				u := particles[i].velocity.minus(particles[j].velocity)
				unit := unitVec(particles[i].position, particles[j].position)
				u = u.multiply(unit)
				magU := u.mag()

				if magU > 0 {
					// V := unitVec(particles[i].position, particles[j].position)
					V := unit.multiplyF(c.timeStep * (1.0 - q) * (c.sigma*magU + c.beta*magU*magU))

					particles[i].velocity.add(V.multiplyF(0.5))
					particles[j].velocity.subtract(V.multiplyF(0.5))

				} // magU > 0?
			} // q < 1?

		} // End loop over j
	} // End loop over i

}

func computeNeighbours(particles []particle, h float64) [][]int {
	n := len(particles)
	result := make([][]int, n)
	h2 := h * h

	for i := 0; i < n; i++ {
		for j := i + 1; j < n; j++ {
			if distSq(particles[i].position, particles[j].position) < h2 {
				result[i] = append(result[i], j)
				result[j] = append(result[j], i)
			}
		}
	}
	return result
}

func update(particles []particle, bb *boundingBox, c constants) {

	// 1. External forces (gravity etc.)
	for i := range particles {
		applyExternalForces(&particles[i], c)
	}

	// 2. Neighbour search
	neighboursArray := computeNeighbours(particles, c.h)

	// 3. Viscosity impulse
	viscosity(c, particles, neighboursArray)

	// 4. Predict positions
	for i := range particles {
		particles[i].integrate(c.timeStep)
	}

	// 5. Density / pressure correction
	for pass := 0; pass < 4; pass++ {
		neighboursArray = computeNeighbours(particles, c.h)
		for i := range particles {
			doubleDensity(c, i, neighboursArray[i], particles)
		}
	}

	// 6. Boundary collisions and velocity update
	for i := range particles {
		bb.collide(&particles[i], c)
		particles[i].updateVelocity(c.timeStep)
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

	// --- simulation parameters ------------------------------------------------
	c := constants{
		g:        vector2{x: 0, y: -9.81},
		epsilon:  1e-3,
		damping:  0.5,
		h:        0.6,
		k:        0.05,
		kNear:    0.20,
		density0: 1.0,
		timeStep: 0.005,
		sigma:    0.02,
		beta:     0.02,
	}

	// --- bounding box (y ↑, x →) ---------------------------------------------
	bb := boundingBox{
		top:    10, // ceiling
		bottom: 0,  // floor
		left:   0,
		right:  10,
	}

	// --- initial particle array ----------------------------------------------
	n := 500
	particles := make([]particle, n)
	for i := 0; i < n; i++ {
		particles[i].position = vector2{x: math.Mod(1+float64(i)*0.5, 10), y: 8}
		// particles[i].velocity = vector2{x: 0, y: 0}
		// particles[i].position = vector2{x: math.Mod(10*math.Sin(float64(i)), 10.), y: 8}
		// x0 := rand.Float64() * 10.0      // uniform in [0,10)
		// y0 := 8.0 + (rand.Float64() * 3) // uniform in [0,10)
		// particles[i].position = vector2{x: x0, y: y0}
		particles[i].velocity = vector2{x: 2 * math.Sin(float64(i)), y: -10}
		particles[i].previousPosition = particles[i].position
	}

	// --- run the simulation ---------------------------------------------------
	steps := 10000
	for step := 0; step < steps; step++ {
		update(particles, &bb, c) // advance
		if step%10 == 0 {
			fname := fmt.Sprintf("out/positions_%04d.csv", step)    // e.g. positions_0000.csv
			if err := savePositions(fname, particles); err != nil { // save snapshot
				panic(err)
			}
		}
	}
}
