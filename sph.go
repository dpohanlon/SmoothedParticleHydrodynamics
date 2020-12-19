package main

import "fmt"
import "math"

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
    g        vector2

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

func (bb * boundingBox) collide(p * particle, c constants) {

    if p.position.y < bb.bottom {

        p.position.y = bb.bottom + c.epsilon
        p.velocity.y = -p.velocity.y * (1 - c.damping)

    } else if p.position.y > bb.top {

        p.position.y = bb.top - c.epsilon
        p.velocity.y = -p.velocity.y * (1 - c.damping)

    } else if p.position.x < bb.left {

        p.position.x = bb.left + c.epsilon
        p.velocity.y = -p.velocity.x * (1 - c.damping)

    } else if p.position.x > bb.right {

        p.position.x = bb.right - c.epsilon
        p.velocity.y = -p.velocity.x * (1 - c.damping)

    }
}

func (v * vector2) mag() float64 {
    return math.Sqrt(v.x * v.x + v.y * v.y)
}

func (v1 * vector2) multiply(v2 vector2) vector2 {
    return vector2{x : v1.x * v2.x, y : v1.y * v2.y}
}

func (v1 * vector2) divide(v2 vector2) vector2 {
    return vector2{x : v1.x / v2.x, y : v1.y / v2.y}
}

func (v1 * vector2) multiplyF(p float64) vector2 {
    return vector2{x : v1.x * p, y : v1.y * p}
}

func (v1 * vector2) divideF(p float64) vector2 {
    return vector2{x : v1.x / p, y : v1.y / p}
}

func (v1 * vector2) minus(v2 vector2) vector2 {
    return vector2{x : v1.x - v2.x, y : v1.y - v2.y}
}

func (v1 * vector2) plus(v2 vector2) vector2 {
    return vector2{x : v1.x + v2.x, y : v1.y + v2.y}
}

func (v1 * vector2) subtract(v2 vector2) {
    v1.x -= v2.x
    v1.y -= v2.y
}

func (v1 * vector2) add(v2 vector2) {
    v1.x += v2.x
    v1.y += v2.y
}

func (v1 * vector2) subtractF(s float64) {
    v1.x -= s
    v1.y -= s
}

func (v1 * vector2) addF(s float64) {
    v1.x += s
    v1.y += s
}


func distSq(p vector2, p2 vector2) float64 {
    p.subtract(p2)
    return math.Pow(p.x, 2) + math.Pow(p.y, 2)
}

func dist(p1 vector2, p2 vector2) float64 {
    return math.Sqrt( distSq(p1, p2) )
}

func unitVec(p1 vector2, p2 vector2) vector2 {
    d := p1.minus(p2)
    return d.divideF(dist(p1, p2))
}

func (p * particle) updateVelocity(timeStep float64) {
    diffPos := p.position.minus(p.previousPosition)
    p.velocity = diffPos.divideF(timeStep)
}

func (p * particle) integrate(timeStep float64) {
    p.previousPosition = p.position
    p.position.add(p.velocity.multiplyF(timeStep))
}

func applyExternalForces(p * particle, c constants) {

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
        q := distSq(particles[this].position, particles[neighbour].position) / (c.h * c.h)
        if q < 1.0 {
            density += (1.0 - q) * (1.0 - q)
            nearDensity += (1.0 - q) * (1.0 - q) * (1.0 - q)
        }
    }

    pressure := c.k * (density - c.density0)
    nearPressure := c.kNear * nearDensity

    deltaX := vector2{x : 0.0, y : 0.0}

    for _, neighbour := range neighbours {
        if this == neighbour {
            continue
        }
        q := distSq(particles[this].position, particles[neighbour].position) / (c.h * c.h)
        if q > 1.0 {
            pressureTerm := pressure * (1.0 - q)
            nearPressureTerm := nearPressure * (1.0 - q) * (1.0 - q)

            D := unitVec(particles[this].position, particles[neighbour].position)
            D.multiplyF(c.timeStep * c.timeStep * (pressureTerm + nearPressureTerm))

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

            q := distSq(particles[i].position, particles[j].position) / (c.h * c.h)

            if q < 1 {
                u := particles[i].velocity.minus(particles[j].velocity)
                u = u.multiply(unitVec(particles[i].position, particles[j].position))
                magU := u.mag()

                if magU > 0 {
                    V := unitVec(particles[i].position, particles[j].position)
                    V = V.multiplyF(c.timeStep * (1.0 - q) * (c.sigma * magU + c.beta * magU * magU))

                    particles[i].velocity.add(V.multiplyF(0.5))
                    particles[j].velocity.subtract(V.multiplyF(0.5))

                } // magU > 0?
            } // q < 1?

        } // End loop over j
    } // End loop over i

}

func update(particles []particle, neighboursArray [][]int, c constants) {

    for _, p := range particles { applyExternalForces(&p, c) }

    viscosity(c, particles, neighboursArray)

    for _, p := range particles { p.integrate(c.timeStep) }

    // for _, n := range neighboursArray { doubleDensity }
}

func main() {

    point1 := vector2{x : 3.0, y : 2.0}
    point2 := vector2{x : 12.0, y : 5.0}

    fmt.Println(point1.mag())
    fmt.Println(dist(point1, point2))
    fmt.Println(point1.minus(point2))
    fmt.Println(unitVec(point1, point2))
}
